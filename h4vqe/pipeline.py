"""Une todo: para cada d -> hamiltoniano, ansatz, θ*, energía medida con shots.

Dueño: Joaquín (integración). Los demás cambian sus módulos, no este flujo,
salvo que se acuerde en el grupo.
"""
from __future__ import annotations

import time
from functools import lru_cache

import numpy as np

from .adapt import ADAPT_FERMIONIC, ADAPT_FERMIONIC_13
from .ansatz import Ansatz, build_ansatz, cz_count
from .config import ANSATZ, ANSATZ_REPS, DISTANCES, SHOTS_PER_POINT, SEED
from .estimator import measure_energy
from .hamiltonian import molecule, qubit_hamiltonian
from .optimize import optimize_params


def curve_ansatz_kwargs() -> dict:
    """Argumentos de build_ansatz para el ansatz de la curva (config.ANSATZ)."""
    if ANSATZ == "adapt":
        return {"excitations": ADAPT_FERMIONIC, "reps": 1}
    if ANSATZ == "uccsd":
        return {"reps": ANSATZ_REPS}
    raise ValueError(f"config.ANSATZ desconocido: {ANSATZ!r}")


@lru_cache(maxsize=8)
def _routed_line_cz(excitations: tuple, norb: int, nelec: tuple) -> int:
    from .routing import LINE_ORDER, line_cz_count, routed_ansatz
    qc, _, _ = routed_ansatz({"norb": norb, "nelec": nelec}, list(excitations), LINE_ORDER)
    return line_cz_count(qc)


def line_cz_fswap(m: dict, ans: Ansatz) -> int | None:
    """CZ con qubits en línea usando el ruteo fSWAP de routing.py (None si no aplica).

    El circuito ruteado prepara el MISMO estado con los modos reordenados (se mide con
    routing.hamiltonian_in_order). La curva se mide con el circuito sin rutear y Ĥ original
    porque Ĥ reordenado necesita ~1,9× más grupos QWC; la energía esperada es la misma.
    Depende solo de la lista de excitaciones, así que se calcula una vez por lista.
    """
    if ans.meta.get("kind", "fermionic") != "fermionic" or not ans.excitations:
        return None
    exc = tuple((tuple(a), tuple(b)) for a, b in ans.excitations)
    return _routed_line_cz(exc, m["norb"], tuple(m["nelec"]))


def run_point(d: float, x0=None, shots: int = SHOTS_PER_POINT, seed: int = SEED,
              ansatz_kwargs=None) -> tuple[dict, dict]:
    t0 = time.time()
    m = molecule(d)
    H = qubit_hamiltonian(m)
    ans = build_ansatz(m, **(curve_ansatz_kwargs() | (ansatz_kwargs or {})))
    if x0 is not None and len(x0) != ans.num_parameters:
        x0 = None
    opt = optimize_params(ans, H, x0=x0, n_restarts=2)
    bound = ans.circuit.assign_parameters(opt["x"])
    meas = measure_energy(bound, H, shots, seed=seed)
    row = {
        "d": d, "E_vqe": meas["energy"], "sigma": meas["sigma"], "shots": meas["shots"],
        "n_groups": meas["n_groups"],
        "cz_all": cz_count(ans.circuit), "cz_line": cz_count(ans.circuit, "line"),
        "cz_line_fswap": line_cz_fswap(m, ans),
        "n_params": ans.num_parameters, "E_opt_exact": opt["energy"],
        "E_hf": m["e_hf"], "E_fci": m["e_fci"], "E_ccsd": m["e_ccsd"],
        "ccsd_converged": m["ccsd_converged"], "seconds": time.time() - t0,
    }
    ctx = {"m": m, "H": H, "ansatz": ans, "opt": opt}
    return row, ctx


# Paso entre semillas de puntos distintos. El estimador usa seed y seed+1 (piloto y Neyman)
# y run_reliability usa offsets < 1000: con un paso de 1000 ninguna semilla se repite.
SEED_STRIDE = 1000


def run_curve(distances=DISTANCES, shots: int = SHOTS_PER_POINT, verbose: bool = True,
              ansatz_kwargs=None, seed: int = SEED):
    """Barre las distancias con warm start (θ* del punto anterior como x0).

    Cada punto usa su PROPIA semilla. Con la misma semilla en todos, el ruido de shots
    queda correlacionado entre distancias vecinas (vimos ~3σ seguidos en 0,7–0,9 Å) y el
    Monte Carlo de R_e, que supone errores independientes, subestima la σ.
    """
    rows, ctxs, x0 = [], {}, None
    for i, d in enumerate(distances):
        row, ctx = run_point(d, x0=x0, shots=shots, seed=seed + SEED_STRIDE * i,
                             ansatz_kwargs=ansatz_kwargs)
        x0 = ctx["opt"]["x"]
        rows.append(row)
        ctxs[d] = ctx
        if verbose:
            print(f"d={d:4.2f} Å  E={row['E_vqe']:.6f} ± {1e3*row['sigma']:.2f} mHa  "
                  f"|E-FCI|={1e3*abs(row['E_vqe']-row['E_fci']):.2f} mHa  "
                  f"|CCSD-FCI|={1e3*abs(row['E_ccsd']-row['E_fci']):.2f} mHa  "
                  f"CZ={row['cz_all']}  ({row['seconds']:.0f}s)")
    return rows, ctxs


def depth_references(m: dict) -> list[tuple[str, str, Ansatz]]:
    """(rol, nombre, ansatz) para la prueba de profundidad, además del ansatz de la curva.

    Con ADAPT: el paso anterior al corte de ADAPT (12 operadores) y UCCSD ×ANSATZ_REPS,
    mucho más profundo. Con UCCSD: UCCSD ×1.
    """
    if ANSATZ == "adapt":
        return [("menos profundo", "ADAPT 12 ops", build_ansatz(m, excitations=ADAPT_FERMIONIC_13[:12])),
                ("más profundo", f"UCCSD ×{ANSATZ_REPS}", build_ansatz(m, reps=ANSATZ_REPS))]
    return [("menos profundo", "UCCSD ×1", build_ansatz(m, reps=1))]


def run_reliability(ctxs: dict, distances, shots: int = SHOTS_PER_POINT, var_shots: int = 2_000_000,
                    echo_shots: int = 1_000_000, sym_shots: int = 200_000, seed: int = SEED,
                    verbose: bool = True) -> list[dict]:
    """Diagnósticos de confiabilidad MEDIDOS con shots, sin usar FCI.

    Para cada d compara: HF (θ = 0), un ansatz menos profundo, el de la curva (ctx) y, con
    ADAPT, uno más profundo (ver depth_references). Para cada estado:
      - energía medida E ± σ (convergencia con la profundidad),
      - Var(H) medida por grupos: ⟨(H−c)²⟩ − (⟨H⟩−c)²,
      - Var(H) por eco de Loschmidt (mucho más sensible),
      - ⟨N_α⟩, ⟨N_β⟩, ⟨S²⟩ medidos.
    `var_exacta` (con el estado) queda solo como chequeo cruzado de las mediciones.
    """
    from .analysis import echo_variance, measured_symmetries, measured_variance, variance
    from .hamiltonian import number_operators

    out = []
    for j, d in enumerate(distances):
        c = ctxs[d]
        m, H = c["m"], c["H"]
        Hm = H.to_matrix(sparse=True)
        ops = number_operators(m["norb"])
        refs = depth_references(m)
        shallow = refs[0][2]
        states = [("HF", "HF", shallow, np.zeros(shallow.num_parameters))]
        for rol, name, ans in refs:
            x = optimize_params(ans, H, n_restarts=2)["x"]
            states.append((rol, f"{name} ({ans.num_parameters} θ)", ans, x))
        states.append(("curva", f"curva ({c['ansatz'].num_parameters} θ)", c["ansatz"], c["opt"]["x"]))
        order = {"HF": 0, "menos profundo": 1, "curva": 2, "más profundo": 3}
        states.sort(key=lambda st: order[st[0]])
        for k, (rol, name, ans, x) in enumerate(states):
            t0 = time.time()
            # semilla propia por (distancia, estado), separada de las de run_curve
            sd = seed + 100 * SEED_STRIDE + SEED_STRIDE * j + 100 * k
            circ = ans.circuit.assign_parameters(x)
            e = measure_energy(circ, H, shots, seed=sd)
            row = {"d": d, "rol": rol, "estado": name, "E": e["energy"], "sigma": e["sigma"],
                   "cz": 0 if rol == "HF" else cz_count(ans.circuit)}
            row |= measured_variance(circ, H, var_shots, e["energy"], e["sigma"], seed=sd + 7)
            row |= echo_variance(circ, H, echo_shots, seed=sd + 11)
            row |= measured_symmetries(circ, ops, sym_shots, seed=sd + 13)
            row["var_exacta"] = variance(ans.fast.state(x), Hm)
            row["seconds"] = time.time() - t0
            out.append(row)
            if verbose:
                print(f"d={d:4.2f} {name:20s} E={row['E']:.6f}±{1e3*row['sigma']:.2f} mHa  "
                      f"Var(grupos)={row['var']:.1e}±{row['var_sigma']:.0e}  "
                      f"Var(eco)={row['var_echo']:.1e} (≤{row['var_echo_up95']:.0e})  "
                      f"N_α={row['N_alpha']:.3f} N_β={row['N_beta']:.3f} S²={row['S2']:.3f}  "
                      f"({row['seconds']:.0f}s)")
    return out
