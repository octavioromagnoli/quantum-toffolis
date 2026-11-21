"""Extra de escalado H4 -> H_N -> FeMoco. Corre los conteos y guarda results/scaling/.

    python scripts/scaling.py                 # N = 4..14, shots exactos hasta N = 12
    python scripts/scaling.py --nmax 12 --exact-max 10

Guarda cada (N, d) apenas termina (results/scaling/recursos.csv) y saltea los ya hechos.
Después ajusta y extrapola a FeMoco (results/scaling/extrapolacion.csv y resumen.txt).
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd

from h4vqe.config import CHEM_ACC
from h4vqe.scaling import (D_EQ, D_STRETCH, chain_resources, extrapolate, log10_p_no_error,
                           uccsd_cz_transpiled)

OUT = ROOT / "results" / "scaling"
CSV = OUT / "recursos.csv"


def run_counts(nmax: int, exact_max: int) -> pd.DataFrame:
    OUT.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(CSV) if CSV.exists() else pd.DataFrame()
    df = df.rename(columns={"cz_uccsd": "cz_cota"})
    for n in range(4, nmax + 1, 2):
        for d in (D_EQ, D_STRETCH):
            have = df[(df.N == n) & np.isclose(df.d, d)] if len(df) else df
            # saltear lo ya hecho, salvo que falten los shots exactos que ahora se piden
            if len(have) and (n > exact_max or have["shots_fci"].notna().all()):
                continue
            row = chain_resources(n, d, exact_shots=n <= exact_max)
            df = pd.concat([df.drop(have.index), pd.DataFrame([row])], ignore_index=True)
            df.to_csv(CSV, index=False)
    cz_csv = OUT / "cz_transpilado.csv"
    cz = pd.read_csv(cz_csv) if cz_csv.exists() else pd.DataFrame(columns=["N", "cz_transpilado"])
    for n in range(4, nmax + 1, 2):
        if n not in set(cz.N):
            t0 = time.time()
            c = uccsd_cz_transpiled(n)
            print(f"H{n}: {c} CZ de UCCSD transpilado ({time.time() - t0:.0f} s)", flush=True)
            cz = pd.concat([cz, pd.DataFrame([{"N": n, "cz_transpilado": c}])], ignore_index=True)
            cz.to_csv(cz_csv, index=False)
    df = df.drop(columns=["cz_transpilado"], errors="ignore").merge(
        cz.astype({"N": int}), on="N", how="left")
    df["log10_p_sin_error"] = [log10_p_no_error(c) for c in df.cz_transpilado]
    df = df.sort_values(["d", "N"]).reset_index(drop=True)
    df.to_csv(CSV, index=False)
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmax", type=int, default=14)
    ap.add_argument("--exact-max", type=int, default=12)
    args = ap.parse_args()
    df = run_counts(args.nmax, args.exact_max)
    ext, lines = extrapolate(df)
    ext.to_csv(OUT / "extrapolacion.csv", index=False)
    cols = ["N", "d", "qubits", "terminos", "grupos_qwc", "params", "cz_transpilado", "cz_cota",
            "log10_p_sin_error", "shots_cota", "shots_fci", "ccsd_menos_fci_mHa",
            "ccsd_converged", "t_fci_s"]
    with pd.option_context("display.width", 200, "display.max_columns", 30,
                           "display.float_format", "{:.4g}".format):
        txt = [str(df[[c for c in cols if c in df]]), *lines, "", str(ext)]
    (OUT / "resumen.txt").write_text("\n".join(txt), encoding="utf-8")
    print("\n".join(txt))
    print(f"\nε = {1e3 * CHEM_ACC} mHa para los shots (σ = precisión química)")


if __name__ == "__main__":
    main()
