"""Reproduce los números del diferencial (ver docs/ideas/diferencial.md).

    python scripts/diferencial_demo.py            # completo (~10 min en una laptop)
    python scripts/diferencial_demo.py --quick    # 5 distancias y menos shots (~2 min)

1. CZ y precisión: UCCSD del equipo contra el ansatz compacto de 5 qubits.
2. Curva completa con el compacto: energía MEDIDA con SamplerV2 (estimator.py del equipo), σ, shots, CZ,
   y <S^2>, <N> también medidos con shots.
3. Nuestro modelo de ruido (noise.py) en la zona estirada: crudo, ZNE (plegado global, λ = 1,3,5) y eco + ZNE (λ = 1,3).
4. Argumento de confiabilidad sin FCI: límite de 4 átomos de H aislados.
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pyscf import gto, scf  # noqa: E402
from qiskit.primitives import StatevectorEstimator  # noqa: E402
from qiskit_aer.primitives import EstimatorV2 as AerEstimator  # noqa: E402

from h4vqe import config  # noqa: E402
from h4vqe.ansatz import build_ansatz, cz_count  # noqa: E402
from h4vqe.estimator import measure_energy  # noqa: E402
from h4vqe.hamiltonian import molecule, qubit_hamiltonian  # noqa: E402
from h4vqe.noise import measure_energy_noisy, noise_model, to_line, zne  # noqa: E402
from h4vqe.optimize import optimize_params  # noqa: E402
from h4vqe.reduction import (build_compact_ansatz, optimize_compact, reduced_number_operators,  # noqa: E402
                             tapering)

QUICK = "--quick" in sys.argv
DISTS = [0.9, 1.8, 2.2, 2.6, 3.0] if QUICK else config.DISTANCES
SHOTS = 1_000_000 if QUICK else config.SHOTS_PER_POINT
SHOTS_SYM = 200_000
NOISY_D = [1.8, 3.0] if QUICK else [1.8, 2.4, 3.0]
SHOTS_NOISY = 2_000_000 if QUICK else 40_000_000
LAYERS = 8
mh = lambda x: 1e3 * x
out = {"quick": QUICK}

# ----------------------------------------------------------------------------- 1. CZ
print("== 1. Costo en CZ (transpilado a cz/sx/x/rz, optimization_level=3)")
m = molecule(2.0)
H8 = qubit_hamiltonian(m)
cz_rows = []
for reps in (1, 2):
    a = build_ansatz(m, reps=reps)
    o = optimize_params(a, H8)
    cz_rows.append({"ansatz": f"UCCSD x{reps} (8q)", "params": a.num_parameters, "cz_all": cz_count(a.circuit),
                    "cz_line": cz_count(a.circuit, "line"), "err_2.0_mHa": mh(o["energy"] - m["e_fci"]),
                    "_circ": a.circuit.assign_parameters(o["x"])})
for L in (6, 7, 8):
    a = build_compact_ansatz(m, layers=L)
    H5 = tapering(m, H8).apply(H8)
    o = optimize_compact(a, H5, reduced_number_operators(m)["S2"])
    cz_rows.append({"ansatz": f"compacto L={L} (5q)", "params": a.num_parameters, "cz_all": cz_count(a.circuit),
                    "cz_line": cz_count(a.circuit, "line"), "err_2.0_mHa": mh(o["energy"] - m["e_fci"]),
                    "_circ": a.circuit.assign_parameters(o["x"]), "_H": H5})
# sesgo exacto por ruido a 2,0 Å (valor de expectación del simulador con ruido, solo para comparar)
est = AerEstimator(options={"backend_options": {"noise_model": noise_model(), "method": "density_matrix"}})
for r in cz_rows:
    c = r.pop("_circ")
    H = r.pop("_H", H8)
    isa = to_line(c, optimization_level=3)
    e = float(est.run([(isa, H.apply_layout(isa.layout))], precision=0).result()[0].data.evs)
    r["err_ruido_2.0_mHa"] = mh(e - m["e_fci"])
    print(f"  {r['ansatz']:22s} params={r['params']:3d}  CZ todos={r['cz_all']:4d}  CZ línea={r['cz_line']:5d}  "
          f"err sin ruido={r['err_2.0_mHa']:8.4f} mHa  err con ruido={r['err_ruido_2.0_mHa']:7.1f} mHa")
print(f"  (referencias a 2,0 Å: HF {mh(m['e_hf'] - m['e_fci']):.1f} mHa, CCSD {mh(m['e_ccsd'] - m['e_fci']):.1f} mHa)")
out["cz"] = cz_rows

# ----------------------------------------------------------------------------- 2. curva
print(f"\n== 2. Curva con el compacto L={LAYERS} (5 qubits, 16 CZ), {SHOTS:,} shots por punto (x3 si d >= 1,8 Å)")
rows, x0, ctx = [], None, {}
t0 = time.time()
for d in DISTS:
    m = molecule(d)
    H8 = qubit_hamiltonian(m)
    tp = tapering(m, H8)
    H5 = tp.apply(H8)
    ops = reduced_number_operators(m)
    ans = build_compact_ansatz(m, layers=LAYERS)
    opt = optimize_compact(ans, H5, ops["S2"], x0=x0, n_random=3)
    x0 = opt["x"]
    bound = ans.circuit.assign_parameters(opt["x"])
    # presupuesto según necesidad: 3x en la zona estirada (σ ≈ 0,5 mHa: precisión química con margen)
    shots_d = SHOTS * (3 if (d >= config.STRETCHED_FROM and not QUICK) else 1)
    r = measure_energy(bound, H5, shots_d, seed=config.SEED + int(100 * d))
    s2 = measure_energy(bound, ops["S2"], SHOTS_SYM, seed=1 + int(100 * d))
    na = measure_energy(bound, ops["N_alpha"], SHOTS_SYM, seed=2 + int(100 * d))
    nb = measure_energy(bound, ops["N_beta"], SHOTS_SYM, seed=3 + int(100 * d))
    row = {"d": d, "E": r["energy"], "sigma": r["sigma"], "shots": r["shots"], "n_groups": r["n_groups"],
           "cz_line": cz_count(ans.circuit, "line"), "E_fci": m["e_fci"], "E_ccsd": m["e_ccsd"],
           "E_hf": m["e_hf"], "E_circuit_exact": opt["energy"],
           "S2": s2["energy"], "S2_sigma": s2["sigma"], "N": na["energy"] + nb["energy"],
           "N_sigma": float(np.hypot(na["sigma"], nb["sigma"]))}
    err, eccsd = abs(row["E"] - row["E_fci"]), abs(row["E_ccsd"] - row["E_fci"])
    row["beats_ccsd"] = bool(err + row["sigma"] < eccsd) if d >= config.STRETCHED_FROM else None
    row["chem_acc"] = bool(err < config.CHEM_ACC)
    rows.append(row)
    ctx[d] = (bound, H5)
    print(f"  d={d:4.2f}  E={row['E']:.6f} ± {mh(row['sigma']):.2f} mHa  |E-FCI|={mh(err):5.2f}  "
          f"|CCSD-FCI|={mh(eccsd):5.2f}  S2={row['S2']:+.3f}±{row['S2_sigma']:.3f}  N={row['N']:.3f}  "
          f"CZ={row['cz_line']}  le gana a CCSD={row['beats_ccsd']}")
st = [r for r in rows if r["d"] >= config.STRETCHED_FROM]
print(f"  zona estirada: le gana a CCSD (con σ) en {sum(r['beats_ccsd'] for r in st)}/{len(st)}; "
      f"precisión química en {sum(r['chem_acc'] for r in st)}/{len(st)}; peor |E-FCI|+σ = "
      f"{max(mh(abs(r['E'] - r['E_fci']) + r['sigma']) for r in st):.2f} mHa  ({time.time() - t0:.0f}s)")
out["curva"] = rows

# ----------------------------------------------------------------------------- 3. ruido
print("\n== 3a. Sesgo por ruido de cada técnica (valor de expectación exacto con ruido: sin shots, solo"
      " para elegir la técnica; NO son energías reportadas)")


def fold_cz_only(c, lam):
    from qiskit import QuantumCircuit
    o = QuantumCircuit(*c.qregs, *c.cregs)
    for ins in c.data:
        for _ in range(lam if ins.operation.name == "cz" else 1):
            o.append(ins.operation, ins.qubits, ins.clbits)
    return o


from h4vqe.noise import echo_fidelity, fold_global, identity_coefficient, richardson  # noqa: E402
bias = []
for d in NOISY_D:
    if d not in ctx:
        continue
    bound, H5 = ctx[d]
    f = molecule(d)["e_fci"]
    isa = to_line(bound)
    ev = lambda c: float(est.run([(c, H5)], precision=0).result()[0].data.evs)
    g = [ev(fold_global(isa, k)) for k in (0, 1, 2)]
    czf = [ev(fold_cz_only(isa, l)) for l in (1, 3, 5)]
    F = echo_fidelity(isa, 2_000_000)["F"]
    emix = identity_coefficient(H5)
    resc = [(e - (1 - F ** l) * emix) / F ** l for e, l in zip(g, (1, 3, 5))]
    b = {"d": d, "crudo": g[0], "ZNE solo CZ (1,3,5)": richardson((1, 3, 5), czf, [0] * 3)[0],
         "ZNE global lineal (1,3)": richardson((1, 3), g[:2], [0] * 2)[0],
         "ZNE global Richardson (1,3,5)": richardson((1, 3, 5), g, [0] * 3)[0],
         "eco solo": resc[0], "eco + ZNE lineal (1,3)": richardson((1, 3), resc[:2], [0] * 2)[0]}
    bias.append(b)
    print(f"  d={d}: " + " | ".join(f"{k} {mh(v - f):+.2f}" for k, v in b.items() if k != "d") + "  (mHa)")
out["sesgo_mitigacion"] = bias

print(f"\n== 3b. Ruido con shots: eco + ZNE lineal (λ = 1,3), {SHOTS_NOISY:,} shots repartidos 3:1")
noisy = []
for d in NOISY_D:
    if d not in ctx:
        continue
    bound, H5 = ctx[d]
    f = molecule(d)["e_fci"]
    isa = to_line(bound)
    raw = measure_energy_noisy(isa, H5, SHOTS_NOISY // 4, seed=11)
    z2 = zne(isa, H5, SHOTS_NOISY, lams=(1, 3), rescale=True)
    rr = {"d": d, "raw": raw, "eco_zne13": z2, "E_fci": f, "E_ccsd": molecule(d)["e_ccsd"]}
    noisy.append(rr)
    print(f"  d={d}: crudo {mh(raw['energy'] - f):6.2f} ± {mh(raw['sigma']):.2f} mHa ({raw['shots']:,} shots) | "
          f"eco+ZNE(1,3) {mh(z2['E'] - f):6.2f} ± {mh(z2['sigma']):.2f} mHa (F={z2['F']:.4f}, "
          f"{z2['shots']:,} shots) | CCSD {mh(rr['E_ccsd'] - f):+.2f} mHa")
out["ruido"] = noisy

# ----------------------------------------------------------------------------- 4. confiabilidad
eH = scf.UHF(gto.M(atom="H 0 0 0", basis=config.BASIS, spin=1, verbose=0)).run().e_tot
far = rows[-1]
print(f"\n== 4. Límite de átomos separados (sin FCI): 4 E_H = {4 * eH:.6f} Ha (1 electrón: exacto)")
print(f"  a {far['d']} Å: VQE - 4E_H = {mh(far['E'] - 4 * eH):+.2f} ± {mh(far['sigma']):.2f} mHa | "
      f"CCSD - 4E_H = {mh(far['E_ccsd'] - 4 * eH):+.2f} mHa  <- CCSD 'liga' 4 átomos casi separados: no físico")
out["4E_H"] = 4 * eH

path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "docs", "ideas",
                    "diferencial_resultados" + ("_quick" if QUICK else "") + ".json")
json.dump(out, open(path, "w"), indent=1, default=float)
print("\nguardado", path)
