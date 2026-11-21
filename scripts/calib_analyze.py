"""Analiza results/calib_sigma/runs.jsonl: tabla de calibración e histogramas de z.

    python scripts/calib_analyze.py [--quick]
Escribe results/calib_sigma/summary.md y hist_z.png.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "calib_sigma")


def std_ci(z, level=0.95):
    """IC de la desviación estándar suponiendo normalidad (χ² con n−1 g.l.)."""
    n = len(z)
    s = np.std(z, ddof=1)
    a = (1 - level) / 2
    return s, s * np.sqrt((n - 1) / stats.chi2.ppf(1 - a, n - 1)), s * np.sqrt((n - 1) / stats.chi2.ppf(a, n - 1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args()
    path = os.path.join(OUT, "runs_quick.jsonl" if args.quick else "runs.jsonl")
    with open(path) as f:
        df = pd.DataFrame([json.loads(l) for l in f if l.strip().endswith("}")])

    rows = []
    for (cfg, d), g in df.groupby(["cfg", "d"]):
        z = g["z"].to_numpy()
        s, lo, hi = std_ci(z)
        sb, lob, hib = std_ci(g["z_boot"].to_numpy())
        rows.append({
            "config": cfg, "d (Å)": d, "reps": len(z), "shots": int(g["shots"].median()),
            "σ media (mHa)": 1e3 * g["sigma"].mean(),
            "media z": f"{z.mean():+.3f} ± {z.std(ddof=1)/np.sqrt(len(z)):.3f}",
            "std z [IC 95%]": f"{s:.3f} [{lo:.3f}, {hi:.3f}]",
            "std z_boot [IC 95%]": f"{sb:.3f} [{lob:.3f}, {hib:.3f}]",
            "σ_boot/σ": f"{(g['sigma_boot']/g['sigma']).mean():.4f}",
            "std(E)/σ media": f"{g['E'].std(ddof=1)/g['sigma'].mean():.3f}",
            "frac(abs z>2)": f"{np.mean(np.abs(z) > 2):.3f}", "n(abs z>3)": int(np.sum(np.abs(z) > 3)),
            "KS p": f"{stats.kstest(z, 'norm').pvalue:.3f}",
            "s/rep": f"{g['seconds'].mean():.1f}",
        })
    tab = pd.DataFrame(rows)
    # Referencia: con 60 repeticiones, el IC 95% de std(z) si la verdadera es 1.
    lo60, hi60 = [np.sqrt(stats.chi2.ppf(q, 59) / 59) for q in (0.025, 0.975)]
    md = ["# Calibración de σ", "",
          "z = (E − E_exacta)/σ con E_exacta = <ψ|H|ψ> del circuito medido. Si σ está bien, z ~ N(0,1):",
          "media 0, std 1 y abs(z)>2 en ~4,6 % de los casos.", "",
          tab.to_markdown(index=False), "",
          f"Referencia: con 60 repeticiones y std verdadera 1, la std muestral cae en "
          f"[{lo60:.2f}, {hi60:.2f}] el 95 % de las veces."]
    with open(os.path.join(OUT, "summary_quick.md" if args.quick else "summary.md"), "w",
              encoding="utf-8") as f:
        f.write("\n".join(md) + "\n")
    print("\n".join(md))

    groups = list(df.groupby(["cfg", "d"]))
    fig, axs = plt.subplots(1, len(groups), figsize=(3.6 * len(groups), 3.2), squeeze=False)
    x = np.linspace(-4, 4, 200)
    for ax, ((cfg, d), g) in zip(axs[0], groups):
        ax.hist(g["z"], bins=np.linspace(-4, 4, 25), density=True, alpha=0.6, label="z")
        ax.plot(x, stats.norm.pdf(x), "k-", lw=1, label="N(0,1)")
        ax.set_title(f"{cfg}, d={d} Å\nn={len(g)}, std={g['z'].std(ddof=1):.2f}", fontsize=9)
        ax.set_xlabel("z = (E − E_exacta)/σ")
    axs[0][0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "hist_z_quick.png" if args.quick else "hist_z.png"), dpi=130)


if __name__ == "__main__":
    main()
