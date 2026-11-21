"""Confirmación CON SHOTS de la mitigación (issue #3, bloque 2). Las energías salen de SamplerV2.

Ansatz compacto de Octavio (5 qubits, 16 CZ), modelo de noise.py + error de lectura de 2 % por
qubit. La matriz de lectura se CALIBRA con shots (noise.readout_calibration), no se usa la exacta.

    python scripts/mitigation_shots.py [--quick]     # ~10 min con 3 procesos (--quick: ~1 min)
"""
from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")

import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # _compact.py (también en los procesos hijos)

import numpy as np

QUICK = "--quick" in sys.argv
DISTS = [1.8, 2.4, 3.0]
P_RO = 0.02
S = 0.1 if QUICK else 1.0
SHOTS = {"crudo": int(10e6 * S), "RO": int(14e6 * S), "PS+RO": int(14e6 * S),
         "ZNE lineal + PS+RO": int(54e6 * S), "eco + ZNE lineal + RO": int(58e6 * S)}
ECHO_SHOTS = int(2e6 * S)
CAL_SHOTS = int(20e6 * S)


def run(d: float) -> dict:
    from _compact import compact_point
    from h4vqe import mitigation as mt
    from h4vqe.estimator import measurement_groups
    from h4vqe.noise import readout_calibration
    t0 = time.time()
    pt = compact_point(d)  # θ* guardado: el mismo circuito que mitigation_study.py
    m, H, ops, isa, E0 = pt["m"], pt["H"], pt["ops"], pt["isa"], pt["E0"]
    sym = [(ops["N_alpha"], 2), (ops["N_beta"], 2)]
    mk = lambda sd: mt.sampler(sd, p_ro=P_RO)
    A = readout_calibration(H.num_qubits, CAL_SHOTS, seed=int(d * 1000) + 5, make_sampler=mk)
    kw = {"p_ro": P_RO}
    seed = int(d * 1000)
    out = {"d": d, "E_fci": m["e_fci"], "E_ccsd": m["e_ccsd"], "E_circ": E0,
           "p_ro_cal": [float(a[1, 0]) for a in A], "rows": {}}
    r = {
        "crudo": mt.measure_energy_mitigated(isa, H, SHOTS["crudo"], seed=seed + 1, **kw),
        "RO": mt.measure_energy_mitigated(isa, H, SHOTS["RO"], seed=seed + 2, readout=A, **kw),
        "PS+RO": mt.measure_energy_mitigated(isa, H, SHOTS["PS+RO"], seed=seed + 3, readout=A,
                                             symmetries=sym, **kw),
        "ZNE lineal + PS+RO": mt.zne_mitigated(isa, H, SHOTS["ZNE lineal + PS+RO"], lams=(1, 3),
                                               seed=seed + 4, readout=A, symmetries=sym, **kw),
        "eco + ZNE lineal + RO": mt.zne_echo(isa, H, SHOTS["eco + ZNE lineal + RO"], lams=(1, 3),
                                             echo_shots=ECHO_SHOTS, seed=seed + 5, readout=A, **kw),
    }
    for k, v in r.items():
        out["rows"][k] = {"E": v["energy"], "sigma": v["sigma"], "shots": v["shots"],
                          "err_vs_fci_mHa": 1e3 * (v["energy"] - m["e_fci"]), "sigma_mHa": 1e3 * v["sigma"]}
    # Error sistemático de la calibración de lectura: misma energía exacta (sin shots) con A
    # calibrada contra A verdadera, para PS+RO (λ = 1).
    groups = measurement_groups(H)
    P = mt.exact_probabilities(isa, groups, p_ro=P_RO)
    e_cal = mt.exact_estimate(P, H, *mt.mitigation_tables(groups, sym, A)[:2])["E"]
    e_true = mt.exact_estimate(P, H, *mt.mitigation_tables(groups, sym, mt.readout_matrices(5, P_RO))[:2])["E"]
    out["readout_cal_systematic_mHa"] = 1e3 * (e_cal - e_true)
    out["seconds"] = time.time() - t0
    return out


def main():
    with ProcessPoolExecutor(len(DISTS)) as ex:
        res = list(ex.map(run, DISTS))
    lines = [f"Shots con ruido (5 qubits, 16 CZ, lectura {P_RO:.0%}); error = E − FCI en mHa ± σ"]
    for r in res:
        lines.append(f"\nd = {r['d']} Å   (CCSD − FCI = {1e3*(r['E_ccsd']-r['E_fci']):+.1f} mHa; "
                     f"p_ro calibrada = {np.mean(r['p_ro_cal']):.4f}; sistemático de la calibración de "
                     f"lectura = {r['readout_cal_systematic_mHa']:+.3f} mHa; {r['seconds']:.0f} s)")
        for k, v in r["rows"].items():
            lines.append(f"  {k:24s} {v['err_vs_fci_mHa']:+9.2f} ± {v['sigma_mHa']:.2f} mHa   "
                         f"({v['shots']/1e6:.0f} M shots)")
    txt = "\n".join(lines)
    print(txt)
    name = "mitigation_shots_quick" if QUICK else "mitigation_shots"
    with open(os.path.join(ROOT, "results", name + ".txt"), "w", encoding="utf-8") as f:
        f.write(txt + "\n")
    with open(os.path.join(ROOT, "results", name + ".json"), "w", encoding="utf-8") as f:
        json.dump(res, f, indent=1)


if __name__ == "__main__":
    main()
