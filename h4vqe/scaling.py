"""Extra «¿y la nitrogenasa?»: cómo crecen los recursos de VQE al alargar la cadena H_N.

Dueño: Santino.

Cadenas lineales equiespaciadas en STO-3G (como H4): N átomos -> N orbitales, 2N qubits,
N electrones (uno por orbital, igual que el modelo (54e, 54o) de FeMoco de Reiher et al.).

Qué se cuenta para cada N (y qué es exacto):
    términos de Ĥ, grupos QWC, Σ|c_k|      -> se construye Ĥ y se cuenta (N ≤ 14)
    parámetros de UCCSD                     -> fórmula exacta (cualquier N)
    CZ de UCCSD con el circuito compacto   -> fórmula exacta (cualquier N), ver uccsd_cost
    shots por evaluación de energía         -> cota (Σ|c|)²/ε² para todo N y valor exacto
                                               con el estado FCI y los grupos reales (N ≤ 10/12)
    E_HF, E_CCSD, E_FCI                     -> PySCF (N ≤ 14)
La extrapolación a FeMoco (N = 54, 76) usa las fórmulas exactas donde existen y un ajuste
log-log donde no (términos, grupos, shots).
"""
from __future__ import annotations

import time
from math import comb

import numpy as np
from pyscf import ao2mo, cc, fci, gto, scf
from pyscf.fci import cistring
from qiskit.quantum_info import SparsePauliOp

from .config import BASIS, CHEM_ACC
from .estimator import split_constant
from .hamiltonian import qubit_hamiltonian

D_EQ = 0.894        # R_e de H4 de la curva principal (Å)
D_STRETCH = 2.0     # zona estirada (Å)
P_CZ = 1e-3         # ruido.py: error por CZ
FEMOCO = {"Reiher 2017": (54, 27, 27), "Li y Chan 2019": (76, 57, 56)}  # (norb, n_alfa, n_beta)


# ---------------------------------------------------------------- molécula y referencias

def chain_geometry(n: int, d: float) -> str:
    return "; ".join(f"H 0 0 {i * d:.6f}" for i in range(n))


def chain_molecule(n: int, d: float, basis: str = BASIS, run_fci: bool = True) -> dict:
    """Integrales en base de orbitales RHF + E_HF, E_FCI (con su vector) y E_CCSD."""
    mol = gto.M(atom=chain_geometry(n, d), basis=basis, unit="Angstrom", verbose=0)
    hf = scf.RHF(mol)
    hf.conv_tol = 1e-10
    hf.max_cycle = 300
    hf.kernel()
    if not hf.converged:                       # segundo intento con Newton (zona estirada)
        hf = hf.newton()
        hf.kernel()
    C = hf.mo_coeff
    norb = C.shape[1]
    h1 = C.T @ hf.get_hcore() @ C
    g = ao2mo.restore(1, ao2mo.kernel(mol, C), norb)
    out = {"N": n, "d": d, "norb": norb, "nelec": (n // 2, n // 2), "h": h1, "g": g,
           "constante": mol.energy_nuc(), "e_hf": hf.e_tot, "hf_converged": bool(hf.converged)}
    if run_fci:
        t0 = time.time()
        solver = fci.direct_spin1.FCI()
        solver.max_space = 8
        solver.conv_tol = 1e-10
        e, ci = solver.kernel(h1, g, norb, out["nelec"], ecore=mol.energy_nuc())
        out.update(e_fci=float(e), ci=ci, t_fci=time.time() - t0)
    t0 = time.time()
    ccsd = cc.CCSD(hf)
    ccsd.conv_tol = 1e-9
    ccsd.max_cycle = 500
    ccsd.kernel()
    out.update(e_ccsd=float(ccsd.e_tot), ccsd_converged=bool(ccsd.converged),
               t_ccsd=time.time() - t0)
    return out


def fci_support(m: dict) -> tuple[np.ndarray, np.ndarray]:
    """Vector FCI de PySCF en la base computacional de los qubits (α en 0..N-1, β en N..2N-1).

    Devuelve (índices, amplitudes) solo sobre el sector N_α = N_β = N/2. Las convenciones de
    orden de PySCF y de Jordan-Wigner difieren a lo sumo en una fase global (el número de
    electrones de cada espín es fijo); se verifica con <ψ|Ĥ|ψ> = E_FCI en scripts/scaling.py.
    """
    norb, (na, nb) = m["norb"], m["nelec"]
    sa = np.asarray(cistring.make_strings(range(norb), na), dtype=np.int64)
    sb = np.asarray(cistring.make_strings(range(norb), nb), dtype=np.int64)
    idx = (sa[:, None] | (sb[None, :] << norb)).ravel()
    return idx, np.asarray(m["ci"], dtype=complex).ravel()


# ---------------------------------------------------------------- mediciones

def qwc_sorted_insertion(rest: SparsePauliOp) -> list[np.ndarray]:
    """Mismo algoritmo que estimator.measurement_groups(H, "sorted", "qwc"), vectorizado.

    El de qiskit (`group_commuting`) arma una matriz T×T×n densa y para H14 pide 20 GiB.
    """
    code = rest.paulis.z.astype(np.int8) + 2 * rest.paulis.x.astype(np.int8)  # 0 I, 1 Z, 2 X, 3 Y
    T, n = code.shape
    bases = np.zeros((T, n), dtype=np.int8)
    members: list[list[int]] = []
    for k in np.argsort(-np.abs(rest.coeffs), kind="stable"):
        p = code[k]
        G = len(members)
        if G:
            B = bases[:G]
            ok = ((B == p) | (B == 0) | (p == 0)).all(axis=1)
            j = int(np.argmax(ok))
            if ok[j]:
                nz = p != 0
                bases[j, nz] = p[nz]
                members[j].append(int(k))
                continue
        bases[G] = p
        members.append([int(k)])
    return [np.array(mm) for mm in members]


def _pauli_masks(rest: SparsePauliOp):
    assert not np.any(rest.paulis.phase), "se esperan Paulis sin fase (la fase va en c_k)"
    nq = rest.num_qubits
    pow2 = np.left_shift(np.int64(1), np.arange(nq, dtype=np.int64))
    x = rest.paulis.x.astype(np.int64) @ pow2
    z = rest.paulis.z.astype(np.int64) @ pow2
    ny = (rest.paulis.x & rest.paulis.z).sum(axis=1)
    return x, z, (1j) ** ny


def group_moments(rest: SparsePauliOp, groups, idx: np.ndarray, amps: np.ndarray):
    """<P_g> y Var(P_g) del estado (idx, amps) para cada grupo P_g = Σ_{k∈g} c_k P_k.

    P|b> = i^{nY} (-1)^{popcount(b & z)} |b XOR x>. Se acumula φ_g = P_g ψ en un buffer del
    tamaño de 2^n y <P_g²> = ||φ_g||², <P_g> = <ψ|φ_g>.
    """
    x, z, iy = _pauli_masks(rest)
    c = rest.coeffs
    buf = np.zeros(1 << rest.num_qubits, dtype=complex)
    mean, var = np.zeros(len(groups)), np.zeros(len(groups))
    for gi, g in enumerate(groups):
        touched = []
        for k in g:
            sign = 1.0 - 2.0 * (np.bitwise_count(idx & z[k]) & 1)   # uint8 -> float
            tgt = idx ^ x[k]
            buf[tgt] += (c[k] * iy[k]) * sign * amps
            touched.append(tgt)
        m1 = float(np.real(np.vdot(amps, buf[idx])))
        m2 = 0.0
        for tgt in touched:          # cada tgt no repite índices; entre listas sí -> poner en 0
            v = buf[tgt]
            m2 += float(np.vdot(v, v).real)
            buf[tgt] = 0.0
        mean[gi], var[gi] = m1, max(m2 - m1 * m1, 0.0)
    return mean, var


def shots_bound(rest: SparsePauliOp, eps: float = CHEM_ACC) -> float:
    """Cota de peor caso: M ≤ (Σ_{k≠I} |c_k|)² / ε² (√Var_g ≤ Σ_{k∈g} |c_k|)."""
    return float(np.sum(np.abs(rest.coeffs)) ** 2 / eps ** 2)


def shots_exact(var_groups: np.ndarray, eps: float = CHEM_ACC) -> float:
    """Shots con reparto de Neyman (n_g ∝ σ_g) para σ = ε: M = (Σ_g σ_g)² / ε²."""
    return float(np.sum(np.sqrt(var_groups)) ** 2 / eps ** 2)


# ---------------------------------------------------------------- ansatz: fórmulas exactas

def _sum_gaps(lo: np.ndarray, hi: np.ndarray) -> int:
    """Σ_{l ∈ lo, h ∈ hi, l < h} (h − l − 1)."""
    d = hi[None, :] - lo[:, None] - 1
    return int(d[d >= 0].sum())


def _sum_gaps_pairs(v: np.ndarray) -> int:
    """Σ_{v_i < v_j} (v_j − v_i − 1) sobre pares de un mismo conjunto."""
    d = v[None, :] - v[:, None] - 1
    return int(np.triu(d, 1).sum())


def uccsd_cost(norb: int, na: int, nb: int) -> dict:
    """Parámetros y CZ de UCCSD×1 con `ansatz.excitation_circuit` (todos con todos).

    El circuito compacto usa 4 CNOT por simple y 14 por doble, más 2 por cada qubit
    espectador cuya ocupación cambia el signo de Jordan-Wigner. Con los 2 o 4 índices de la
    excitación ordenados p1 < p2 (< p3 < p4), esos qubits son los que están estrictamente
    entre p1 y p2 (y entre p3 y p4): (p2 − p1 − 1) [+ (p4 − p3 − 1)]. Las simples y dobles
    αα/ββ caen en un mismo bloque de espín; en las αβ, todo α (0..N-1) va antes que β.
    Verificado contra excitation_circuit para N = 4 y 6 (scripts/scaling.py).
    """
    occ = [np.arange(na), np.arange(nb)]
    vac = [np.arange(na, norb), np.arange(nb, norb)]
    n_s = [len(occ[s]) * len(vac[s]) for s in range(2)]
    gap_s = [_sum_gaps(occ[s], vac[s]) for s in range(2)]          # simples: a − i − 1
    singles = sum(n_s)
    cz = 4 * singles + 2 * sum(gap_s)
    # dobles αβ: i < a (α) < j < b (β)  ->  (a − i − 1) + (b − j − 1)
    d_ab = n_s[0] * n_s[1]
    cz += 14 * d_ab + 2 * (gap_s[0] * n_s[1] + gap_s[1] * n_s[0])
    # dobles del mismo espín: i < j < a < b  ->  (j − i − 1) + (b − a − 1)
    d_same = 0
    for s in range(2):
        no, nv = len(occ[s]), len(vac[s])
        cnt = comb(no, 2) * comb(nv, 2)
        d_same += cnt
        cz += 14 * cnt + 2 * (_sum_gaps_pairs(occ[s]) * comb(nv, 2)
                              + _sum_gaps_pairs(vac[s]) * comb(no, 2))
    return {"params": singles + d_ab + d_same, "singles": singles,
            "doubles": d_ab + d_same, "cz_cota": int(cz)}


def uccsd_cz_transpiled(norb: int, coupling: str = "all") -> int:
    """CZ de UCCSD×1 como en el desafío central: PauliEvolutionGate por excitación y
    `ansatz.cz_count` (transpile nivel 3). Cancela compuertas entre cadenas de Pauli vecinas,
    así que da menos que uccsd_cost (×1,25 en H4, ×1,9 en H10). Medido hasta N = 14."""
    from qiskit import QuantumCircuit
    from qiskit.circuit.library import PauliEvolutionGate

    from .ansatz import cz_count, excitation_generator, hf_occupation, uccsd_excitations
    from .hamiltonian import annihilation_ops

    m = {"norb": norb, "nelec": (norb // 2, norb // 2)}
    nq = 2 * norb
    a = annihilation_ops(nq)
    qc = QuantumCircuit(nq)
    qc.x(hf_occupation(m))
    for e in uccsd_excitations(m):
        qc.append(PauliEvolutionGate(excitation_generator(e, nq, a), time=0.1), range(nq))
    return cz_count(qc, coupling)


def log10_p_no_error(n_cz: int, p: float = P_CZ) -> float:
    """log10 de P(ninguna CZ falla) = (1 − p)^n_CZ. Ignora el error de 1 qubit (20× menor),
    así que es una cota superior de la probabilidad de circuito limpio."""
    return float(n_cz * np.log10(1.0 - p))


# ---------------------------------------------------------------- un punto completo

def chain_resources(n: int, d: float, exact_shots: bool = True, verbose: bool = True) -> dict:
    """Todo lo de la tabla para una cadena H_n a distancia d."""
    t0 = time.time()
    m = chain_molecule(n, d)
    H = qubit_hamiltonian(m)
    const, rest = split_constant(H)
    t_h = time.time() - t0 - m["t_fci"] - m["t_ccsd"]
    t1 = time.time()
    groups = qwc_sorted_insertion(rest)
    t_g = time.time() - t1
    norb, (na, nb) = m["norb"], m["nelec"]
    cost = uccsd_cost(norb, na, nb)
    row = {
        "N": n, "d": d, "qubits": 2 * norb,
        "dim_fci": comb(norb, na) * comb(norb, nb),
        "statevector_bytes": 16 * 2 ** (2 * norb),
        "terminos": len(rest), "grupos_qwc": len(groups),
        "suma_abs_c": float(np.sum(np.abs(rest.coeffs))),
        **cost,
        "log10_p_sin_error": log10_p_no_error(cost["cz_cota"]),
        "shots_cota": shots_bound(rest),
        "e_hf": m["e_hf"], "e_fci": m["e_fci"], "e_ccsd": m["e_ccsd"],
        "ccsd_menos_fci_mHa": 1e3 * (m["e_ccsd"] - m["e_fci"]),
        "ccsd_converged": m["ccsd_converged"], "hf_converged": m["hf_converged"],
        "t_fci_s": m["t_fci"], "t_ccsd_s": m["t_ccsd"], "t_ham_s": t_h, "t_grupos_s": t_g,
    }
    if exact_shots:
        t1 = time.time()
        idx, amps = fci_support(m)
        mean, var = group_moments(rest, groups, idx, amps)
        row["e_check_mHa"] = 1e3 * (const + mean.sum() - m["e_fci"])   # tiene que dar ~0
        row["shots_fci"] = shots_exact(var)
        row["t_shots_s"] = time.time() - t1
    if verbose:
        print(f"H{n:<2d} d={d:.3f}: {row['terminos']} términos, {row['grupos_qwc']} grupos, "
              f"{cost['params']} parámetros, {cost['cz_cota']} CZ (cota), "
              f"CCSD−FCI = {row['ccsd_menos_fci_mHa']:+.2f} mHa  "
              f"(FCI {m['t_fci']:.1f}s, grupos {t_g:.1f}s)", flush=True)
    return row


# ---------------------------------------------------------------- ajuste y extrapolación

def powerlaw_fit(N, y):
    """log y = a + b log N. Devuelve (a, b, cov)."""
    (b, a), cov = np.polyfit(np.log(N), np.log(y), 1, cov=True)
    return a, b, cov


def powerlaw_predict(a, b, cov, N):
    """Predicción y banda de 1σ (solo la estadística del ajuste, no el error de modelo)."""
    x = np.log(N)
    J = np.array([x, 1.0])                       # d(log y)/d(b, a), en el orden de polyfit
    s = float(np.sqrt(J @ cov @ J))
    mu = a + b * x
    return float(np.exp(mu)), float(np.exp(mu - s)), float(np.exp(mu + s))


FIT_MIN_N = 6        # los ajustes usan N ≥ 6 (H4 todavía no está en el régimen asintótico)
FIT_COLUMNS = ["terminos", "grupos_qwc", "suma_abs_c", "cz_transpilado", "cz_cota",
               "shots_cota", "shots_fci"]


def extrapolate(df, fit_min_n: int = FIT_MIN_N):
    """Ajuste log-log de cada magnitud contra N (N ≥ fit_min_n) y extrapolación a FeMoco.

    Devuelve (tabla por (d, modelo de FeMoco), líneas de texto con los exponentes). Parámetros y
    CZ de la cota salen de la fórmula exacta; términos, grupos, CZ transpiladas y shots, del ajuste
    (banda de 1σ: solo la estadística del ajuste, no el error de modelo).
    """
    import pandas as pd

    lines, rows = [], []
    for d in sorted(df.d.unique()):
        s = df[np.isclose(df.d, d) & (df.N >= fit_min_n)]
        lines.append(f"d = {d} Å  (ajuste log-log con N = {sorted(s.N.astype(int).tolist())})")
        fits = {}
        for col in FIT_COLUMNS:
            ss = s.dropna(subset=[col]) if col in s else s.iloc[0:0]
            if len(ss) < 3:
                continue
            a, b, cov = powerlaw_fit(ss.N.values, ss[col].values.astype(float))
            fits[col] = (a, b, cov)
            lines.append(f"  {col:15s} ~ N^{b:.2f} ± {np.sqrt(cov[0, 0]):.2f}   (N = {sorted(ss.N.astype(int).tolist())})")
        for name, (norb, na, nb) in FEMOCO.items():
            cost = uccsd_cost(norb, na, nb)
            r = {"d": d, "modelo": name, "N_orb": norb, "qubits": 2 * norb,
                 "dim_fci": comb(norb, na) * comb(norb, nb), **cost}
            for col, (a, b, cov) in fits.items():
                if col == "cz_cota":          # exacto por fórmula, no hace falta ajustar
                    continue
                mu, lo, hi = powerlaw_predict(a, b, cov, norb)
                r[col], r[col + "_lo"], r[col + "_hi"] = mu, lo, hi
            if "cz_transpilado" in r:
                r["log10_p_sin_error"] = log10_p_no_error(r["cz_transpilado"])
            r["log10_p_sin_error_cota"] = log10_p_no_error(cost["cz_cota"])
            rows.append(r)
        if "shots_fci" in fits:
            sb = s.dropna(subset=["shots_fci"])
            ratio = (sb.shots_fci / sb.shots_cota).values
            lines.append(f"  shots exactos / cota (N = {sb.N.astype(int).tolist()}): {np.round(ratio, 3).tolist()}")
    return pd.DataFrame(rows), lines


def plot_scaling(df, path=None):
    """(a) recursos de VQE contra N, con los ajustes extrapolados a FeMoco; (b) error de CCSD."""
    import matplotlib.pyplot as plt

    eq = df[np.isclose(df.d, D_EQ)].sort_values("N")
    st = df[np.isclose(df.d, D_STRETCH)].sort_values("N")
    colors = ["tab:blue", "tab:orange", "tab:green", "tab:red"]
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(13, 4.8), gridspec_kw={"width_ratios": [1.5, 1]})
    Nx = np.arange(4, 81, 2)
    cz = [uccsd_cost(n, n // 2, n // 2)["cz_cota"] for n in Nx]
    series = [("términos de Ĥ", eq.terminos, "o"),
              ("grupos QWC", eq.grupos_qwc, "s"),
              ("CZ de UCCSD ×1 (transpilado)", eq.cz_transpilado, "D"),
              ("shots por energía (exactos, σ = 1,6 mHa)", eq.shots_fci, "^")]
    for (lab, y, mk), col in zip(series, colors):
        ok = y.notna()
        n, y = eq.N[ok], y[ok].astype(float)
        m = n >= FIT_MIN_N
        a, b, _ = powerlaw_fit(n[m].values, y[m].values)
        ax.plot(n, y, mk, color=col, ms=7, label=f"{lab}  ~N^{b:.1f}", zorder=3)
        ax.plot(Nx, np.exp(a) * Nx.astype(float) ** b, "--", color=col, lw=1.3)
        if "CZ" in lab:   # la fórmula exacta del circuito compacto: cota superior para todo N
            ax.plot(Nx, cz, ":", color=col, lw=1.3, label="CZ: cota exacta (fórmula)")
    for n, txt, ha, dx in [(54, "FeMoco (54e, 54o)\n108 qubits", "right", -3),
                           (76, "FeMoco (113e, 76o)\n152 qubits", "left", 3)]:
        ax.axvline(n, color="0.4", lw=1, ls=":")
        ax.annotate(txt, (n, 3e14), xytext=(dx, 0), textcoords="offset points", ha=ha, va="top", fontsize=8.5)
    ax.axvline(4, color="0.4", lw=1, ls=":")
    ax.annotate("H₄ (este trabajo)\n8 qubits", (4, 3e14), xytext=(3, 0), textcoords="offset points",
                ha="left", va="top", fontsize=8.5)
    ax.set(xscale="log", yscale="log", xlabel="N (átomos de H = orbitales; qubits = 2N)",
           ylabel="cantidad", ylim=(10, 1e15), xlim=(3.6, 110))
    ax.set_xticks([4, 6, 8, 10, 14, 20, 30, 54, 76])
    ax.set_xticklabels(["4", "6", "8", "10", "14", "20", "30", "54", "76"])
    ax.minorticks_off()
    ax.grid(True, color="0.9")
    ax.legend(loc="upper left", bbox_to_anchor=(0.02, 0.86), fontsize=8, frameon=False)
    ax.set_title(f"(a) Recursos de VQE (puntos: medidos; rayas: ajuste con N ≥ {FIT_MIN_N})", fontsize=10, loc="left")
    for (lab, s), col, mk in zip([(f"d = {D_EQ} Å (equilibrio): CCSD por encima de FCI", eq),
                                  (f"d = {D_STRETCH} Å (estirada): CCSD por debajo de FCI", st)],
                                 colors, "os"):
        err = np.abs(s.ccsd_menos_fci_mHa)
        bx.plot(s.N, err, "-" + mk, color=col, lw=2, ms=7, label=lab)
        bad = ~s.ccsd_converged.astype(bool)
        if bad.any():
            bx.plot(s.N[bad], err[bad], mk, ms=9, mfc="white", mec=col, mew=2, zorder=3)
            for n_, e_ in zip(s.N[bad], err[bad]):
                bx.annotate("no converge", (n_, e_), xytext=(0, 9), textcoords="offset points",
                            ha="center", fontsize=8.5)
    bx.axhline(1e3 * CHEM_ACC, color="0.4", lw=1, ls="--")
    bx.text(4, 1e3 * CHEM_ACC * 1.25, "precisión química (1,6 mHa)", ha="left", va="bottom", fontsize=8.5)
    bx.set(yscale="log", xlabel="N", ylabel="|E_CCSD − E_FCI| (mHa)", ylim=(3e-3, 1e3))
    bx.set_xticks(sorted(df.N.unique()))
    bx.grid(True, color="0.9")
    bx.legend(loc="lower right", fontsize=8, frameon=False)
    bx.set_title("(b) Error de CCSD contra FCI", fontsize=10, loc="left")
    fig.tight_layout()
    if path:
        fig.savefig(path, dpi=150)
    return fig
