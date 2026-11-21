"""Tabla, R_e, D_e, comparación contra CCSD, gráficos y diagnósticos.

Dueño: Joaquín.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .config import CHEM_ACC, STRETCHED_FROM, SEED


def results_table(rows: list[dict]) -> pd.DataFrame:
    df = pd.DataFrame(rows).sort_values("d").reset_index(drop=True)
    df["err_vqe_mHa"] = 1e3 * (df["E_vqe"] - df["E_fci"]).abs()
    df["err_ccsd_mHa"] = 1e3 * (df["E_ccsd"] - df["E_fci"]).abs()
    df["sigma_mHa"] = 1e3 * df["sigma"]
    stretched = df["d"] >= STRETCHED_FROM - 1e-9
    df["beats_ccsd"] = np.where(stretched, df["err_vqe_mHa"] + df["sigma_mHa"] < df["err_ccsd_mHa"], np.nan)
    df["chem_acc"] = df["err_vqe_mHa"] < 1e3 * CHEM_ACC
    return df


def fit_minimum(d, E, sigma, n_points: int = 5, n_mc: int = 2000, seed: int = SEED):
    """R_e y E(R_e) con un ajuste cuadrático sobre los n_points más bajos.
    Incertidumbre por Monte Carlo sobre las barras de error."""
    d, E, sigma = map(np.asarray, (d, E, sigma))
    idx = np.argsort(E)[:n_points]
    idx = np.sort(idx)
    x, y, s = d[idx], E[idx], sigma[idx]

    def one(yy):
        c = np.polyfit(x, yy, 2, w=1 / np.maximum(s, 1e-12))
        r = -c[1] / (2 * c[0])
        return r, np.polyval(c, r)

    re, emin = one(y)
    rng = np.random.default_rng(seed)
    mc = np.array([one(y + rng.normal(0, s)) for _ in range(n_mc)])
    return {"R_e": re, "R_e_sigma": mc[:, 0].std(), "E_min": emin, "E_min_sigma": mc[:, 1].std(),
            "fit_points": x.tolist()}


def dissociation_energy(E_far: float, sigma_far: float, fit: dict) -> dict:
    De = E_far - fit["E_min"]
    return {"D_e": De, "D_e_sigma": float(np.hypot(sigma_far, fit["E_min_sigma"]))}


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


def variance(psi: np.ndarray, Hm) -> float:
    """⟨H²⟩ − ⟨H⟩² (diagnóstico: → 0 cerca de un autoestado)."""
    hpsi = Hm @ psi
    e = np.real(np.vdot(psi, hpsi))
    return float(np.real(np.vdot(hpsi, hpsi)) - e ** 2)
