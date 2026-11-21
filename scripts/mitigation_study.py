"""Cuánto mejora cada técnica de mitigación (issue #3, bloque 2), SIN shots: sesgo exacto y costo.

Circuito: ansatz compacto de Octavio (5 qubits, 16 CZ, en línea), modelo de ruido de noise.py
(despolarizante 0,1 % por CZ y 0,005 % por compuerta de 1 qubit), con y sin error de lectura.
Para cada técnica:
  sesgo = E(infinitos shots, con ruido y mitigación) − E(circuito sin ruido), en mHa
  shots para σ = 0,5 mHa = (Σ_i |γ_i| √K_i)² / σ²  (K_i = σ²·N con Neyman en el nivel λ_i;
  γ_i = coeficientes de Richardson; reparto óptimo de shots entre niveles)
Las probabilidades salen de matriz densidad (exactas). Las energías reportadas salen de shots.

    python scripts/mitigation_study.py            # ~2 min
"""
from __future__ import annotations

import json
import os
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import numpy as np
from qiskit.quantum_info import SparsePauliOp, Statevector

from _compact import compact_point
from h4vqe import mitigation as mt
from h4vqe.estimator import measurement_groups, split_constant

DISTS = [1.8, 2.4, 3.0]
P_RO = [0.0, 0.02]
SIG = 0.5e-3


def gammas(lams):
    lams = np.asarray(lams, float)
    return np.array([np.prod([lj / (lj - li) for j, lj in enumerate(lams) if j != i])
                     for i, li in enumerate(lams)])


def grouping_costs(H: SparsePauliOp, psi: np.ndarray) -> dict:
    """K = (Σ_g √V_g)² con el estado exacto, para cada agrupamiento (sin ruido)."""
    out = {}
    for name, method, grouping in [("QWC coloreo", "coloring", "qwc"), ("QWC sorted", "sorted", "qwc"),
                                   ("GC coloreo", "coloring", "gc"), ("GC sorted", "sorted", "gc")]:
        gs = measurement_groups(H, method=method, grouping=grouping)
        sv = Statevector(psi)
        V = [max(sv.expectation_value(g.dot(g)).real - sv.expectation_value(g).real ** 2, 0) for g in gs]
        out[name] = {"groups": len(gs), "K": float(np.sum(np.sqrt(V)) ** 2)}
    return out


def main():
    results = {}
    for d in DISTS:
        pt = compact_point(d)  # θ* guardado: el mismo circuito que mitigation_shots.py
        m, H, ops, bound, isa, E0 = pt["m"], pt["H"], pt["ops"], pt["bound"], pt["isa"], pt["E0"]
        n = H.num_qubits
        sym = [(ops["N_alpha"], 2), (ops["N_beta"], 2)]
        groups = measurement_groups(H)
        F = mt.echo_fidelity_exact(isa)
        emix = split_constant(H)[0]
        res = {"E_circ_minus_fci_mHa": 1e3 * (E0 - m["e_fci"]), "cz": int(isa.count_ops()["cz"]),
               "F_echo": F, "grouping": grouping_costs(H, Statevector(bound).data), "rows": {}}
        print(f"\n=== d = {d} Å  (5 qubits, {res['cz']} CZ; E_circuito − FCI = {res['E_circ_minus_fci_mHa']:.1e} mHa;"
              f" F eco = {F:.4f})")
        for p_ro in P_RO:
            P = {lam: mt.exact_probabilities(isa, groups, p_ro=p_ro, lam=lam) for lam in (1, 3, 5)}
            # como noise.zne: plegando solo el ansatz (las rotaciones de base quedan sin escalar)
            P_old = {lam: mt.exact_probabilities(isa, groups, p_ro=p_ro, lam=lam, fold_rotations=False)
                     for lam in (1, 3, 5)}
            ro = mt.readout_matrices(n, p_ro) if p_ro > 0 else None
            tabs = {"crudo": mt.mitigation_tables(groups),
                    "RO": mt.mitigation_tables(groups, readout=ro),
                    "PS": mt.mitigation_tables(groups, symmetries=sym),
                    "PS+RO": mt.mitigation_tables(groups, symmetries=sym, readout=ro)}
            est = {k: {lam: mt.exact_estimate(P[lam], H, U, W) for lam in (1, 3, 5)}
                   for k, (U, W, _) in tabs.items()}
            k0 = "RO" if p_ro > 0 else "crudo"
            est["viejo"] = {lam: mt.exact_estimate(P_old[lam], H, *tabs[k0][:2]) for lam in (1, 3, 5)}
            rows = {}

            def add(name, base, lams, rescale=False):
                g = gammas(lams)
                E = np.array([est[base][l]["E"] for l in lams])
                K = np.array([est[base][l]["K"] for l in lams])
                if rescale:  # eco de Octavio: E' = (E − (1 − F^λ) E_mix) / F^λ, σ' = σ / F^λ
                    fl = F ** np.asarray(lams, float)
                    E, K = (E - (1 - fl) * emix) / fl, K / fl ** 2
                e = float(g @ E)
                cost = float(np.sum(np.abs(g) * np.sqrt(K)) ** 2)
                rows[name] = {"bias_mHa": 1e3 * (e - E0), "shots_M_sigma0.5": cost / SIG ** 2 / 1e6}

            for k in tabs:
                add(k, k, (1,))
            for k in ("RO", "PS+RO") if p_ro > 0 else ("crudo", "PS"):
                add(f"ZNE lineal (1,3) + {k}", k, (1, 3))
                add(f"ZNE Richardson (1,3,5) + {k}", k, (1, 3, 5))
            add(f"eco + ZNE lineal (1,3) + {k0}  (plegando R·C)", k0, (1, 3), rescale=True)
            add(f"ZNE lineal (1,3) plegando solo C  [noise.zne]", "viejo", (1, 3))
            add(f"Richardson (1,3,5) plegando solo C  [noise.zne]", "viejo", (1, 3, 5))
            add(f"eco + lineal plegando solo C  [noise.zne]", "viejo", (1, 3), rescale=True)
            res["rows"][f"p_ro={p_ro}"] = rows
            print(f"  -- lectura p_ro = {p_ro:.0%} --")
            print(f"  {'técnica':48s} {'sesgo (mHa)':>12s} {'shots σ=0,5':>12s}")
            for name, r in rows.items():
                print(f"  {name:48s} {r['bias_mHa']:+12.3f} {r['shots_M_sigma0.5']:10.1f} M")
        print("  -- agrupamiento sin ruido (5 qubits), shots para σ = 0,5 mHa --")
        for name, r in res["grouping"].items():
            print(f"  {name:12s} {r['groups']:3d} grupos   K = {r['K']:.2f}   {r['K']/SIG**2/1e6:5.1f} M")
        results[d] = res
    os.makedirs(os.path.join(ROOT, "results"), exist_ok=True)
    with open(os.path.join(ROOT, "results", "mitigation_study.json"), "w", encoding="utf-8") as f:
        json.dump(results, f, indent=1, ensure_ascii=False)


if __name__ == "__main__":
    main()
