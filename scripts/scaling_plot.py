"""Gráfico del extra de escalado: results/scaling/escalado.png (lee recursos.csv).

El gráfico lo arma h4vqe.scaling.plot_scaling, el mismo que usa el notebook.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import matplotlib
matplotlib.use("Agg")
import pandas as pd

from h4vqe.scaling import plot_scaling

OUT = ROOT / "results" / "scaling"


def main():
    plot_scaling(pd.read_csv(OUT / "recursos.csv"), OUT / "escalado.png")
    print("guardado", OUT / "escalado.png")


if __name__ == "__main__":
    main()
