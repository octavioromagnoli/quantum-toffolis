"""Herramientas para contar la historia del notebook: de lo mínimo al resultado final.

Dueño: Santino (integración del notebook de entrega).

Nada de este módulo da energías REPORTADAS: las de la curva salen de medir con shots
(estimator.py). Acá hay diagnósticos con el estado exacto o con matriz densidad, que explican
POR QUÉ tomamos cada decisión. La consigna lo permite: los parámetros se pueden buscar sin ruido
y los métodos clásicos sirven para comparar.

    hamiltonian_structure(Hs)  clases de términos de Pauli de Ĥ y su peso a cada distancia
    top_terms(H, n)            los n términos más grandes, en notación legible
    hamiltonian_latex(H, n)    Ĥ escrito con sus n términos más grandes (para display(Math(...)))
    minimum_table(distances)   lo mínimo: UCCSD ×1 y ×2 sin ruido (error, S², CZ, P sin errores)
    adapt_convergence(...)     ADAPT al ir agregando operadores: error, varianza y S²
    shots_table(distances)     shots para σ = 0,5 mHa con cada forma de agrupar
    noisy_ansatz_bias(d)       sesgo EXACTO con el modelo de ruido (matriz densidad), por ansatz
    isolated_atoms_energy()    energía de 4 átomos de H separados (límite físico, sin FCI)
    plot_evolution(evo)        gráfico del resumen final
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .adapt import ADAPT_FERMIONIC, ADAPT_FERMIONIC_13
from .ansatz import build_ansatz, cz_count
from .config import BASIS, BASIS_GATES, SEED
from .estimator import measurement_groups
from .hamiltonian import molecule, number_operators, qubit_hamiltonian
from .optimize import optimize_params

P_CZ = 1e-3   # error por CZ del modelo de ruido (noise.py)

CLASS_INFO = {
    "I": "constante (no se mide)",
    "Z": "Z_p: energía de ocupar el espín-orbital p",
    "ZZ": "Z_p Z_q: repulsión entre dos ocupaciones (Coulomb e intercambio)",
    "XY2": "2 X/Y (+ cadenas Z): un electrón salta de un orbital a otro",
    "XY4": "4 X/Y (+ cadenas Z): dos electrones saltan a la vez",
}
_ORDER = list(CLASS_INFO)


# --------------------------------------------------------------------------- el hamiltoniano
def pauli_class(label: str) -> str:
    """Clase de una cadena de Pauli (ver CLASS_INFO): cuántas X/Y tiene, o cuántas Z si es diagonal."""
    n_xy = sum(ch in "XY" for ch in label)
    if n_xy:
        return f"XY{n_xy}"
    return {0: "I", 1: "Z", 2: "ZZ"}.get(label.count("Z"), "Z" * label.count("Z"))


def readable(label: str) -> str:
    """'IIIZIIIZ' (Qiskit es little-endian: el último carácter es el qubit 0) -> 'Z0 Z4'."""
    return " ".join(f"{ch}{q}" for q, ch in enumerate(reversed(label)) if ch != "I") or "I"


def hamiltonian_structure(Hs: dict) -> pd.DataFrame:
    """Por clase de término: cuántos hay, Σ|c| a cada distancia y el más grande a la última."""
    ds = list(Hs)
    rows: dict[str, dict] = {}
    for d in ds:
        for lab, c in Hs[d].to_list():
            k = pauli_class(lab)
            r = rows.setdefault(k, {"clase": k, "qué representa": CLASS_INFO.get(k, ""), "términos": 0})
            if d == ds[0]:
                r["términos"] += 1
            key = f"Σ|c| a {d:g} Å (Ha)"
            r[key] = r.get(key, 0.0) + abs(c)
            if d == ds[-1] and abs(c) > r.get("_max", -1):
                r["_max"] = abs(c)
                r[f"el más grande a {d:g} Å"] = f"{np.real(c):+.4f} · {readable(lab)}"
    out = sorted(rows.values(), key=lambda r: _ORDER.index(r["clase"]) if r["clase"] in _ORDER else 99)
    return pd.DataFrame(out).drop(columns="_max")


def top_terms(H, n: int = 12) -> pd.DataFrame:
    lst = sorted(H.to_list(), key=lambda t: -abs(t[1]))[:n]
    return pd.DataFrame([{"término": readable(lab), "clase": pauli_class(lab), "c (Ha)": float(np.real(c))}
                         for lab, c in lst])


def hamiltonian_latex(H, n: int = 10) -> str:
    """Ĥ = Σ c_k P_k con sus n términos más grandes, en LaTeX."""
    lst = sorted(H.to_list(), key=lambda t: -abs(t[1]))
    parts = []
    for i, (lab, c) in enumerate(lst[:n]):
        ops = "\\,".join(f"{ch}_{{{q}}}" for q, ch in enumerate(reversed(lab)) if ch != "I") or "I"
        c = float(np.real(c))
        sign = ("-" if c < 0 else "") if i == 0 else ("-" if c < 0 else "+")
        parts.append(f"{sign}{abs(c):.4f}\\,{ops}")
    return (f"\\hat H = {' '.join(parts)} + \\cdots"
            f"\\qquad ({len(H)}\\ \\text{{términos en total}})")


def isolated_atoms_energy(n_atoms: int = 4, basis: str = BASIS) -> float:
    """n átomos de H infinitamente separados: n × E(H). Con un electrón, UHF es exacto en la base."""
    from pyscf import gto, scf
    mol = gto.M(atom="H 0 0 0", basis=basis, spin=1, verbose=0)
    return n_atoms * float(scf.UHF(mol).kernel())


# --------------------------------------------------------------------------- lo mínimo y ADAPT
def _references(distances):
    ms = {d: molecule(d) for d in distances}
    Hs = {d: qubit_hamiltonian(m) for d, m in ms.items()}
    S2 = number_operators(ms[distances[0]]["norb"])["S2"].to_matrix(sparse=True)
    return ms, Hs, S2


def minimum_table(distances=(1.8, 2.0, 2.4, 3.0), reps=(1, 2), restarts=None) -> pd.DataFrame:
    """Lo mínimo (la primera versión, UCCSD con Trotter): por ansatz y distancia, error SIN ruido
    contra FCI (con el estado exacto, solo para elegir el ansatz), ⟨S²⟩, CZ y la probabilidad de
    que ninguna CZ falle con 0,1 % por CZ y qubits en línea.

    restarts: arranques al azar extra por cantidad de pasos de Trotter (además del warm start y de
    θ = 0). UCCSD ×1 los necesita para encontrar su mejor mínimo; UCCSD ×2 llega a FCI sin ellos,
    y con 52 parámetros son lo más caro de la celda."""
    distances = list(distances)
    restarts = {1: 2, 2: 0} if restarts is None else restarts
    ms, Hs, S2 = _references(distances)
    rows = []
    for r in reps:
        ans = build_ansatz(ms[distances[0]], reps=r)      # el circuito no depende de d (solo los θ)
        cz = (cz_count(ans.circuit), cz_count(ans.circuit, "line"))
        x0 = None
        for d in distances:
            o = optimize_params(ans, Hs[d], x0=x0, n_restarts=restarts.get(r, 0))
            x0 = o["x"]
            psi = ans.fast.state(o["x"])
            rows.append({"ansatz": f"UCCSD ×{r}", "d": d, "parámetros": ans.num_parameters,
                         "CZ todos con todos": cz[0], "CZ en línea": cz[1],
                         "P(ninguna CZ falla), línea": (1 - P_CZ) ** cz[1],
                         "error sin ruido (mHa)": 1e3 * (o["energy"] - ms[d]["e_fci"]),
                         "S²": float(np.real(np.vdot(psi, S2 @ psi))),
                         "|CCSD − FCI| (mHa)": 1e3 * abs(ms[d]["e_ccsd"] - ms[d]["e_fci"])})
    return pd.DataFrame(rows)


def exc_label(exc, norb: int = 4) -> str:
    """((1, 5), (3, 7)) -> '1α 1β → 3α 3β' (qubits 0..norb-1 = α, norb..2norb-1 = β)."""
    name = lambda q: f"{q % norb}{'α' if q < norb else 'β'}"
    src, dst = exc
    return " ".join(map(name, src)) + " → " + " ".join(map(name, dst))


def adapt_list(k: int) -> list:
    """Los primeros k operadores en el orden en que ADAPT los eligió (k ≤ 13), o la lista final
    de 15 con los pares de espín (ADAPT_FERMIONIC)."""
    return ADAPT_FERMIONIC_13[:k] if k <= len(ADAPT_FERMIONIC_13) else ADAPT_FERMIONIC[:k]


def adapt_convergence(distances, ks=(8, 10, 12, 13, 15)) -> pd.DataFrame:
    """ADAPT al ir agregando operadores. Por (k, d): error sin ruido, Var(H) y S², con el estado
    exacto. La varianza no usa FCI: es el criterio con el que cortamos la lista."""
    distances = list(distances)
    ms, Hs, S2 = _references(distances)
    Hm = {d: H.to_matrix(sparse=True) for d, H in Hs.items()}
    rows = []
    for k in ks:
        ans = build_ansatz(ms[distances[0]], excitations=adapt_list(k))
        cz = cz_count(ans.circuit)
        x0 = None
        for d in distances:
            o = optimize_params(ans, Hs[d], x0=x0)
            x0 = o["x"]
            psi = ans.fast.state(o["x"])
            hpsi = Hm[d] @ psi
            e = float(np.real(np.vdot(psi, hpsi)))
            rows.append({"operadores": k, "d": d, "CZ": cz, "error (mHa)": 1e3 * (e - ms[d]["e_fci"]),
                         "Var(H) (Ha²)": float(np.real(np.vdot(hpsi, hpsi))) - e ** 2,
                         "S²": float(np.real(np.vdot(psi, S2 @ psi)))})
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- shots
def _neyman_K(H, psi, method: str, grouping: str) -> tuple[int, float]:
    """(grupos, K) con K = (Σ_g √Var_g)²: con el reparto de Neyman, σ² = K / N."""
    groups = measurement_groups(H, method=method, grouping=grouping)
    tot = 0.0
    for g in groups:
        v = g.to_matrix(sparse=True) @ psi
        tot += np.sqrt(max(float(np.vdot(v, v).real - np.vdot(psi, v).real ** 2), 0.0))
    return len(groups), tot ** 2


GROUPINGS = [("qwc", "coloring", "QWC, coloreo (lo mínimo)"),
             ("qwc", "sorted", "QWC, sorted insertion (la curva)"),
             ("gc", "coloring", "GC, coloreo"),
             ("gc", "sorted", "GC, sorted insertion")]


def shots_table(distances, sigma: float = 0.5e-3, shots_ref: float = 8e6) -> pd.DataFrame:
    """Shots para σ objetivo con cada forma de agrupar, en el estado del ansatz de la curva.
    Solo para PLANIFICAR: la σ que se reporta sale siempre de las muestras."""
    rows, ans = [], None
    for d in distances:
        m = molecule(d)
        H = qubit_hamiltonian(m)
        ans = ans or build_ansatz(m, excitations=ADAPT_FERMIONIC)
        psi = ans.fast.state(optimize_params(ans, H)["x"])
        for grouping, method, name in GROUPINGS:
            n, K = _neyman_K(H, psi, method, grouping)
            rows.append({"d": d, "agrupamiento": name, "grupos": n,
                         "shots para σ = 0,5 mHa (M)": K / sigma ** 2 / 1e6,
                         "σ con 8 M de shots (mHa)": 1e3 * np.sqrt(K / shots_ref)})
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- ruido
def _line_isa(circ):
    from qiskit import transpile
    from qiskit.transpiler import CouplingMap
    n = circ.num_qubits
    return transpile(circ, basis_gates=BASIS_GATES, coupling_map=CouplingMap.from_line(n),
                     initial_layout=list(range(n)), optimization_level=3, seed_transpiler=SEED)


def noisy_ansatz_bias(d: float = 2.4, verbose: bool = True) -> tuple[pd.DataFrame, dict]:
    """Los tres candidatos a ansatz con el modelo de ruido de noise.py (qubits en línea), sin shots:
    probabilidades exactas con matriz densidad (mitigation.exact_probabilities). Mide el SESGO que
    deja el ruido, que con shots no se va. Devuelve la tabla y lo que reusa la mitigación
    (molécula, Ĥ₅ y θ* del compacto, del MISMO proceso)."""
    import time
    from qiskit.quantum_info import Statevector
    from . import mitigation as mt
    from .noise import to_line
    from .reduction import build_compact_ansatz, optimize_compact, reduced_number_operators, tapering
    from .routing import LINE_ORDER, hamiltonian_in_order, routed_ansatz

    m = molecule(d)
    H = qubit_hamiltonian(m)
    H5 = tapering(m, H).apply(H)
    S2_5 = reduced_number_operators(m)["S2"]
    comp = build_compact_ansatz(m)
    oc = optimize_compact(comp, H5, S2_5)
    ad = build_ansatz(m, excitations=ADAPT_FERMIONIC)
    oa = optimize_params(ad, H)
    qc_r, order, _ = routed_ansatz(m, ADAPT_FERMIONIC, LINE_ORDER)
    uc = build_ansatz(m, reps=2)
    ou = optimize_params(uc, H)
    cases = [("compacto 5 qubits (L = 8)", to_line(comp.circuit.assign_parameters(oc["x"])), H5, oc["energy"]),
             ("ADAPT 15 + ruteo fSWAP", _line_isa(qc_r.assign_parameters(oa["x"])),
              hamiltonian_in_order(m, order), oa["energy"]),
             ("UCCSD ×2 (lo mínimo)", _line_isa(uc.circuit.assign_parameters(ou["x"])), H, ou["energy"])]
    f = m["e_fci"]
    rows = []
    for name, isa, Hx, e_exact in cases:
        t0 = time.time()
        Hx = Hx.apply_layout(isa.layout) if isa.layout is not None else Hx
        e_ideal = float(Statevector(isa).expectation_value(Hx).real)   # sin ruido: chequea el ruteo
        groups = measurement_groups(Hx)
        U, W, _ = mt.mitigation_tables(groups)
        e_noisy = mt.exact_estimate(mt.exact_probabilities(isa, groups), Hx, U, W)["E"]
        ops = isa.count_ops()
        n1 = int(sum(ops.get(g, 0) for g in ("rz", "sx", "x")))
        n2 = int(ops.get("cz", 0))
        rows.append({"ansatz": name, "qubits": isa.num_qubits, "CZ (línea)": n2, "compuertas de 1 qubit": n1,
                     "P(ningún error)": (1 - P_CZ) ** n2 * (1 - 5e-5) ** n1,
                     "error sin ruido (mHa)": 1e3 * (e_ideal - f),
                     "error con ruido, sin mitigar (mHa)": 1e3 * (e_noisy - f)})
        if verbose:
            print(f"{name:28s} {n2:5d} CZ  sin ruido {1e3*(e_ideal-f):+8.3f} mHa  "
                  f"con ruido {1e3*(e_noisy-f):+8.1f} mHa  ({time.time()-t0:.0f} s)")
    rows.append({"ansatz": "HF (sin circuito)", "qubits": 8, "CZ (línea)": 0, "compuertas de 1 qubit": 4,
                 "P(ningún error)": 1.0, "error sin ruido (mHa)": 1e3 * (m["e_hf"] - f),
                 "error con ruido, sin mitigar (mHa)": 1e3 * (m["e_hf"] - f)})
    ctx = {"d": d, "m": m, "H": H, "H5": H5, "S2_5": S2_5, "compact": comp, "x_compact": oc["x"]}
    return pd.DataFrame(rows), ctx


# --------------------------------------------------------------------------- resumen final
def plot_evolution(evo: pd.DataFrame, chem_acc_mHa: float = 1.6):
    """Dos paneles por etapa: CZ con qubits en línea y error a 2,4 Å (sin ruido y con ruido)."""
    import matplotlib.pyplot as plt
    x = np.arange(len(evo))
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 4.4))
    bars = a1.bar(x, evo["CZ en línea"], color=["0.6"] + ["tab:blue"] * (len(evo) - 1))
    a1.bar_label(bars, labels=[f"{v:.0f}" for v in evo["CZ en línea"]], padding=2, fontsize=9)
    a1.set_yscale("log")
    a1.set_ylim(5, 5e3)
    a1.set_ylabel("CZ con qubits en línea (escala log)")
    a1.set_title("Costo del circuito")
    a2.axhspan(0, chem_acc_mHa, color="tab:green", alpha=0.2, label="precisión química (1,6 mHa)")
    a2.plot(x, evo["error sin ruido (mHa)"].abs().clip(lower=1e-3), "o-", label="sin ruido (θ*, exacto)")
    a2.plot(x, evo["error con ruido (mHa)"].abs(), "s-", color="tab:red", label="con ruido 0,1 %/CZ")
    a2.axhline(evo["|CCSD − FCI| (mHa)"].iloc[0], color="tab:orange", ls="--", label="CCSD")
    a2.set_yscale("log")
    a2.set_ylabel("|E − E_FCI| a 2,4 Å (mHa, escala log)")
    a2.set_title("Error a 2,4 Å")
    a2.legend(fontsize=8, loc="upper right")
    for a in (a1, a2):
        a.set_xticks(x)
        a.set_xticklabels(evo["etapa"], fontsize=8)
    fig.tight_layout()
    return fig
