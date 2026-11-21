"""Extra "comparar ansätze CON ruido" (issue #2): cuánto pesa la cantidad de CZ.

Mismo punto estirado, tres ansätze que SIN ruido llegan a ≤ 0,4 mHa de FCI:
  compacto 5q (L = 8, 16 CZ) · ADAPT 15 + ruteo fSWAP (8q) · UCCSD ×2 (transpiler común).
Ruido: modelo de noise.py (qubits en línea, 0,1 % por CZ, 0,005 % por compuerta de 1 qubit).
Cada circuito se transpila a la línea con layout inicial trivial; si el transpiler agrega SWAP,
Ĥ se mapea con H.apply_layout(isa.layout) (se verifica sin ruido que da la energía exacta).
ZNE: plegado global con λ = 1, 3 (noise.zne), solo para los circuitos de < 500 CZ.

    python scripts/noise_ansatze.py [d] [shots]     # default 2,4 Å y 500 000 (~10 min)
Dueño: Gonzalo.
"""
from __future__ import annotations

import json
import os
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from qiskit import transpile
from qiskit.quantum_info import Statevector
from qiskit.transpiler import CouplingMap

from h4vqe.adapt import ADAPT_FERMIONIC
from h4vqe.ansatz import build_ansatz
from h4vqe.config import BASIS_GATES, SEED
from h4vqe.hamiltonian import molecule, qubit_hamiltonian
from h4vqe.noise import measure_energy_noisy, to_line, zne
from h4vqe.optimize import optimize_params
from h4vqe.reduction import build_compact_ansatz, optimize_compact, reduced_number_operators, tapering
from h4vqe.routing import LINE_ORDER, hamiltonian_in_order, routed_ansatz


def line_isa(circ):
    n = circ.num_qubits
    return transpile(circ, basis_gates=BASIS_GATES, coupling_map=CouplingMap.from_line(n),
                     initial_layout=list(range(n)), optimization_level=3, seed_transpiler=SEED)


def main(d: float = 2.4, shots: int = 500_000):
    m = molecule(d)
    H = qubit_hamiltonian(m)
    fci = m["e_fci"]
    cases = []
    tp = tapering(m, H)
    H5 = tp.apply(H)
    c = build_compact_ansatz(m, layers=8)
    o = optimize_compact(c, H5, reduced_number_operators(m)["S2"], n_random=4)
    cases.append(("compacto 5q L=8", c.circuit.assign_parameters(o["x"]), H5, o["energy"]))
    a = build_ansatz(m, excitations=ADAPT_FERMIONIC)
    oa = optimize_params(a, H, n_restarts=2)
    qc, order, _ = routed_ansatz(m, ADAPT_FERMIONIC, LINE_ORDER)
    cases.append(("ADAPT 15 + fSWAP", qc.assign_parameters(oa["x"]), hamiltonian_in_order(m, order),
                  oa["energy"]))
    u = build_ansatz(m, reps=2)
    ou = optimize_params(u, H, n_restarts=1)
    cases.append(("UCCSD ×2", u.circuit.assign_parameters(ou["x"]), H, ou["energy"]))

    rows = []
    for name, circ, Hc, e_exact in cases:
        isa = line_isa(circ)
        Hphys = Hc.apply_layout(isa.layout)
        check = Statevector(isa).expectation_value(Hphys).real
        assert abs(check - e_exact) < 1e-8, (name, check, e_exact)
        r = measure_energy_noisy(isa, Hphys, shots, seed=SEED)
        ops = isa.count_ops()
        row = {"ansatz": name, "d": d, "cz": int(ops.get("cz", 0)),
               "compuertas_1q": int(sum(ops.get(g, 0) for g in ("sx", "x", "rz"))),
               "err_sin_ruido_mHa": 1e3 * (e_exact - fci),
               "err_con_ruido_mHa": 1e3 * (r["energy"] - fci), "sigma_con_ruido_mHa": 1e3 * r["sigma"],
               "shots": r["shots"]}
        if row["cz"] < 500:
            # ZNE sobre el circuito de to_line (nivel 1), porque noise.fold_global necesita que
            # el inverso re-transpilado tenga las mismas compuertas.
            isa1 = to_line(circ)
            try:
                z = zne(isa1, Hc.apply_layout(isa1.layout), 4 * shots, lams=(1, 3), seed=SEED + 5)
                row |= {"err_zne_mHa": 1e3 * (z["E"] - fci), "sigma_zne_mHa": 1e3 * z["sigma"],
                        "cz_zne": int(isa1.count_ops().get("cz", 0)), "shots_zne": z["shots"]}
            except AssertionError:
                row["zne"] = "falló: el inverso re-transpilado no tiene las mismas compuertas"
        rows.append(row)
        print({k: (round(v, 2) if isinstance(v, float) else v) for k, v in row.items()}, flush=True)
    print(f"CCSD − FCI = {1e3 * (m['e_ccsd'] - fci):.1f} mHa")
    out = os.path.join(ROOT, "scripts", f"noise_ansatze_{d}.json")
    with open(out, "w", encoding="utf-8") as f:
        json.dump(rows, f, ensure_ascii=False, indent=1)
    print("escrito", out)


if __name__ == "__main__":
    main(float(sys.argv[1]) if len(sys.argv) > 1 else 2.4,
         int(sys.argv[2]) if len(sys.argv) > 2 else 500_000)
