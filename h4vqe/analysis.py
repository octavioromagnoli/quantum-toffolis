"""Tabla, R_e, D_e, comparación contra CCSD, gráficos y diagnósticos de confiabilidad.

Dueño: Joaquín.

Todo lo que se reporta sale de medir circuitos con shots (SamplerV2 de Aer).
Lo que usa el estado exacto (`variance`) queda solo como chequeo cruzado.
"""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
from qiskit import QuantumCircuit, transpile
from qiskit.circuit.library import PauliEvolutionGate
from qiskit.quantum_info import SparsePauliOp
from qiskit.synthesis import SuzukiTrotter
from qiskit_aer.primitives import SamplerV2 as AerSampler
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize_scalar
from scipy.stats import beta

from .config import BASIS_GATES, CHEM_ACC, STRETCHED_FROM, SEED
from .estimator import measure_energy, split_constant


# ---------------------------------------------------------------- tabla y gráfico

def results_table(rows: list[dict]) -> pd.DataFrame:
    df = pd.DataFrame(rows).sort_values("d").reset_index(drop=True)
    df["err_vqe_mHa"] = 1e3 * (df["E_vqe"] - df["E_fci"]).abs()
    df["err_ccsd_mHa"] = 1e3 * (df["E_ccsd"] - df["E_fci"]).abs()
    df["sigma_mHa"] = 1e3 * df["sigma"]
    stretched = df["d"] >= STRETCHED_FROM - 1e-9
    df["beats_ccsd"] = np.where(stretched, df["err_vqe_mHa"] + df["sigma_mHa"] < df["err_ccsd_mHa"], np.nan)
    df["chem_acc"] = df["err_vqe_mHa"] < 1e3 * CHEM_ACC
    return df


def plot_curve(df: pd.DataFrame, path: str | None = None):
    import matplotlib.pyplot as plt
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))
    ax1.plot(df.d, df.E_fci, "k-", label="FCI (exacta)")
    ax1.plot(df.d, df.E_ccsd, "s--", ms=4, label="CCSD")
    ax1.plot(df.d, df.E_hf, ":", label="HF")
    ax1.errorbar(df.d, df.E_vqe, yerr=df.sigma, fmt="o", ms=4, capsize=3, label="VQE (shots)")
    ax1.set_xlabel("d (Å)"); ax1.set_ylabel("E (Ha)"); ax1.legend(); ax1.set_title("H₄ lineal, STO-3G")
    ax2.errorbar(df.d, df.err_vqe_mHa, yerr=df.sigma_mHa, fmt="o", capsize=3, label="|VQE − FCI|")
    ax2.plot(df.d, df.err_ccsd_mHa, "s--", label="|CCSD − FCI|")
    ax2.axhline(1e3 * CHEM_ACC, color="g", ls=":", label="precisión química")
    ax2.axvspan(STRETCHED_FROM, df.d.max(), color="0.9")
    ax2.set_yscale("log"); ax2.set_xlabel("d (Å)"); ax2.set_ylabel("error (mHa)"); ax2.legend()
    fig.tight_layout()
    if path:
        fig.savefig(path, dpi=150)
    return fig


# ---------------------------------------------------------------- R_e y D_e
#
# Modelos para el mínimo: (tipo, grado, cantidad de puntos alrededor del punto más bajo).
# Elegimos el default aplicando cada modelo a las energías FCI de la MISMA grilla y
# comparando con el mínimo FCI verdadero (minimizando PySCF): el cuadrático de 5 puntos
# (el que usábamos) tiene un sesgo de +8 mÅ en R_e, mayor que su σ estadística.
# El cuártico sobre 7 puntos (0,7–1,1 Å) tiene sesgo < 1 mÅ y < 20 µHa en E_min.
FIT_MODELS = {
    "cuadrático 3": ("poly", 2, 3),
    "cuadrático 5": ("poly", 2, 5),
    "cúbico 6": ("poly", 3, 6),
    "cúbico 7": ("poly", 3, 7),
    "cuártico 7": ("poly", 4, 7),
    "cuártico 8": ("poly", 4, 8),
    "spline 0,6–1,6 Å": ("spline", 3, None),
}
DEFAULT_FIT = "cuártico 7"


def _window(d, E, n):
    """Índices de los n puntos de la grilla más cercanos al punto de E más baja."""
    k = int(np.argmin(E))
    order = np.argsort(np.abs(d - d[k]) + 1e-9 * np.arange(len(d)), kind="stable")
    return np.sort(order[:n])


def _fit_once(d, E, sigma, model):
    kind, deg, n = FIT_MODELS[model]
    if kind == "spline":
        m = (d >= 0.6 - 1e-9) & (d <= 1.6 + 1e-9)
        if m.sum() < 4:
            return np.nan, np.nan
        cs = CubicSpline(d[m], E[m])
        k = int(np.argmin(E[m]))
        lo, hi = d[m][max(k - 1, 0)], d[m][min(k + 1, m.sum() - 1)]
        r = minimize_scalar(cs, bounds=(lo, hi), method="bounded")
        return float(r.x), float(cs(r.x))
    n = min(n, len(d))
    deg = min(deg, n - 1)
    if deg < 2:
        return np.nan, np.nan
    i = _window(d, E, n)
    c = np.polyfit(d[i], E[i], deg, w=1 / np.maximum(sigma[i], 1e-12))
    roots = np.roots(np.polyder(c))
    roots = roots[np.abs(roots.imag) < 1e-9].real
    roots = roots[(roots > d[i].min()) & (roots < d[i].max()) & (np.polyval(np.polyder(c, 2), roots) > 0)]
    if len(roots) == 0:
        return np.nan, np.nan
    r = roots[np.argmin(np.polyval(c, roots))]
    return float(r), float(np.polyval(c, r))


def fit_minimum(d, E, sigma, model: str = DEFAULT_FIT, n_mc: int = 2000, seed: int = SEED) -> dict:
    """R_e y E(R_e) con el modelo `model` (ver FIT_MODELS).
    Incertidumbre por Monte Carlo: se re-muestrean las energías con sus σ y se re-ajusta
    (incluida la elección de la ventana de puntos)."""
    d, E, sigma = (np.asarray(v, float) for v in (d, E, sigma))
    re, emin = _fit_once(d, E, sigma, model)
    rng = np.random.default_rng(seed)
    mc = np.array([_fit_once(d, E + rng.normal(0, sigma), sigma, model) for _ in range(n_mc)])
    with np.errstate(all="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        s_re, s_e = float(np.nanstd(mc[:, 0])), float(np.nanstd(mc[:, 1]))
    kind, _, n = FIT_MODELS[model]
    used = d[_window(d, E, min(n, len(d)))] if kind == "poly" else d[(d >= 0.6 - 1e-9) & (d <= 1.6 + 1e-9)]
    return {"model": model, "R_e": re, "R_e_sigma": s_re, "E_min": emin,
            "E_min_sigma": s_e, "fit_points": used.tolist(),
            "mc_failures": int(np.isnan(mc[:, 0]).sum())}


def dissociation_energy(E_far: float, sigma_far: float, fit: dict) -> dict:
    """D_e = E(3,0 Å) − E(R_e). Los dos son independientes (circuitos distintos)."""
    De = E_far - fit["E_min"]
    return {"D_e": De, "D_e_sigma": float(np.hypot(sigma_far, fit["E_min_sigma"]))}


def fci_minimum(d0: float = 0.89) -> dict:
    """Mínimo FCI verdadero (minimizando PySCF en d). SOLO como referencia para comparar."""
    from .hamiltonian import molecule
    r = minimize_scalar(lambda x: molecule(x)["e_fci"], bracket=(d0 - 0.05, d0, d0 + 0.05), tol=1e-6)
    return {"R_e": float(r.x), "E_min": float(r.fun)}


def compare_fits(df: pd.DataFrame, fci_ref: dict | None = None, n_mc: int = 1000) -> pd.DataFrame:
    """Aplica todos los modelos a VQE (con σ) y a FCI en la misma grilla.

    sesgo = (modelo aplicado a FCI) − (mínimo FCI verdadero): mide el error del MÉTODO de
    ajuste, independiente del VQE. Sirve para elegir el modelo sin mirar los datos del VQE.
    """
    d, far = df.d.to_numpy(), df[df.d == df.d.max()].iloc[0]
    tiny = np.full(len(df), 1e-6)
    out = []
    for model in FIT_MODELS:
        v = fit_minimum(d, df.E_vqe, df.sigma, model, n_mc=n_mc)
        if np.isnan(v["R_e"]):
            continue
        de = dissociation_energy(far.E_vqe, far.sigma, v)
        row = {"modelo": model, "n_puntos": len(v["fit_points"]),
               "R_e (Å)": v["R_e"], "σ R_e (mÅ)": 1e3 * v["R_e_sigma"],
               "D_e (mHa)": 1e3 * de["D_e"], "σ D_e (mHa)": 1e3 * de["D_e_sigma"]}
        if fci_ref is not None:
            f_re, f_e = _fit_once(d, df.E_fci.to_numpy(), tiny, model)
            row["sesgo R_e (mÅ)"] = 1e3 * (f_re - fci_ref["R_e"])
            row["sesgo E_min (µHa)"] = 1e6 * (f_e - fci_ref["E_min"])
            row["R_e − R_e(FCI) (mÅ)"] = 1e3 * (v["R_e"] - fci_ref["R_e"])
            row["D_e − D_e(FCI) (mHa)"] = 1e3 * (de["D_e"] - (far.E_fci - fci_ref["E_min"]))
        out.append(row)
    return pd.DataFrame(out)


# ---------------------------------------------------------------- confiabilidad (medida)

def variance(psi: np.ndarray, Hm) -> float:
    """⟨H²⟩ − ⟨H⟩² con el estado exacto. SOLO chequeo cruzado de las versiones medidas."""
    hpsi = Hm @ psi
    e = np.real(np.vdot(psi, hpsi))
    return float(np.real(np.vdot(hpsi, hpsi)) - e ** 2)


def measured_variance(circuit: QuantumCircuit, H: SparsePauliOp, shots: int, E: float, sigma_E: float,
                      seed: int = SEED) -> dict:
    """Var(H) = ⟨(H−c)²⟩ − (⟨H⟩−c)², con ⟨(H−c)²⟩ MEDIDO con shots (grupos QWC del estimador).

    c = término constante de H: minimiza (casi) la σ de la estimación (probado a 2,0 Å:
    σ es la mitad que con c = E). E y σ_E vienen de la medición de energía (independiente).
    """
    c, _ = split_constant(H)
    K = H - c * SparsePauliOp("I" * H.num_qubits)
    K2 = K.dot(K).simplify(atol=1e-12)
    r = measure_energy(circuit, K2, shots, seed=seed)
    var = r["energy"] - (E - c) ** 2
    sig = float(np.hypot(r["sigma"], 2 * abs(E - c) * sigma_E))
    return {"var": float(var), "var_sigma": sig, "shots": r["shots"], "n_groups": r["n_groups"]}


def echo_circuit(circuit: QuantumCircuit, H: SparsePauliOp, t: float, trotter_reps: int = 4) -> QuantumCircuit:
    """U(θ) → e^{−iHt} (Suzuki-Trotter de 2º orden) → U(θ)† → medir todo."""
    _, H0 = split_constant(H)  # la constante solo da una fase global
    ev = PauliEvolutionGate(H0, time=t, synthesis=SuzukiTrotter(order=2, reps=trotter_reps))
    qc = circuit.copy()
    qc.append(ev, range(H.num_qubits))
    qc.compose(circuit.inverse(), inplace=True)
    qc.measure_all()
    return transpile(qc, basis_gates=BASIS_GATES + ["measure"], optimization_level=1, seed_transpiler=SEED)


def echo_variance(circuit: QuantumCircuit, H: SparsePauliOp, shots: int, t: float = 1.0,
                  trotter_reps: int = 4, seed: int = SEED) -> dict:
    """Eco de Loschmidt: 1 − |⟨ψ|e^{−iHt}|ψ⟩|² = Var(H)·t² + O(t⁴).

    Se prepara ψ, se evoluciona y se des-prepara: la probabilidad de NO volver a |0…0⟩
    estima Var·t². Es mucho más sensible que medir ⟨H²⟩ por grupos, porque cuenta eventos
    raros (≈ 0 para un autoestado) en vez de restar dos números grandes.
    Por Mandelstam-Tamm, 1−P₀ ≤ Var·t²: subestima un poco (≈ 8 % a t = 1 para UCCSD×1).
    Cota superior del 95 % por Clopper-Pearson (con 0 eventos, ≈ 3/shots).
    Error de Trotter (2º orden, t = 1, a 2,0 Å, estado convergido): 1−P₀ = 1,3e-6 con 2 pasos
    y 9e-8 con 4 (∝ 1/n⁴). Con 4 pasos queda debajo de la resolución de 1M shots.
    """
    qc = echo_circuit(circuit, H, t, trotter_reps)
    counts = AerSampler(seed=seed).run([(qc, None, int(shots))]).result()[0].data.meas.get_counts()
    k = shots - counts.get("0" * H.num_qubits, 0)
    p = k / shots
    up = beta.ppf(0.95, k + 1, shots - k) if k < shots else 1.0
    return {"var_echo": p / t ** 2, "var_echo_sigma": np.sqrt(max(p * (1 - p), 1 / shots) / shots) / t ** 2,
            "var_echo_up95": up / t ** 2, "echo_events": int(k), "echo_shots": int(shots),
            "echo_cz": int(qc.count_ops().get("cz", 0))}


def measured_symmetries(circuit: QuantumCircuit, ops: dict, shots: int, seed: int = SEED) -> dict:
    """⟨N_α⟩, ⟨N_β⟩, ⟨S²⟩ MEDIDOS con shots (mismo estimador que la energía)."""
    out = {}
    for name, op in ops.items():
        r = measure_energy(circuit, op, shots, seed=seed)
        out[name], out[name + "_sigma"] = r["energy"], r["sigma"]
    return out
