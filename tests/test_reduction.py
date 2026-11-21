"""Chequeos de la reducción a 5 qubits, el ansatz compacto y el ruido.
Correr con:  python tests/test_reduction.py   (o python -m pytest -q tests)"""
import os
import sys

import numpy as np
from qiskit.quantum_info import Operator, Statevector

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from h4vqe.ansatz import cz_count  # noqa: E402
from h4vqe.estimator import measure_energy  # noqa: E402
from h4vqe.hamiltonian import molecule, qubit_hamiltonian  # noqa: E402
from h4vqe.noise import fold_global, measure_energy_noisy, to_line  # noqa: E402
from h4vqe.reduction import (build_compact_ansatz, optimize_compact, reduced_number_operators,  # noqa: E402
                             tapering)


def test_reduction_and_compact():
    for d in (0.9, 2.0, 3.0):
        m = molecule(d)
        H = qubit_hamiltonian(m)
        tp = tapering(m, H)
        H5 = tp.apply(H)
        assert H5.num_qubits == 5
        # 1. mismo E0 que FCI y mismo <HF|H|HF> que PySCF
        M = H5.to_matrix()
        assert abs(np.linalg.eigvalsh(M)[0] - m["e_fci"]) < 1e-9
        assert abs(M[tp.hf_index, tp.hf_index].real - m["e_hf"]) < 1e-9
        ans = build_compact_ansatz(m, layers=8)
        # 2. simulador numpy == circuito Qiskit; gradiente adjunto == diferencias finitas
        x = np.random.default_rng(1).uniform(-1, 1, ans.num_parameters)
        fid = abs(np.vdot(Statevector(ans.circuit.assign_parameters(x)).data, ans.fast.state(x))) ** 2
        assert fid > 1 - 1e-9, fid
        Hm = H5.to_matrix(sparse=True)
        e0, g = ans.fast.energy_and_grad(x, Hm)
        k, h = 7, 1e-6
        xp = x.copy(); xp[k] += h
        assert abs((ans.fast.energy_and_grad(xp, Hm)[0] - e0) / h - g[k]) < 1e-4
        # 3. VQE: variacional, exacto y en el sector correcto (singlete, N = 4)
        ops = reduced_number_operators(m)
        opt = optimize_compact(ans, H5, ops["S2"])
        assert opt["energy"] >= m["e_fci"] - 1e-9
        assert opt["energy"] - m["e_fci"] < 1e-5, 1e3 * (opt["energy"] - m["e_fci"])
        sv = Statevector(ans.fast.state(opt["x"]))
        assert abs(sv.expectation_value(ops["S2"]).real) < 1e-3
        assert abs(sv.expectation_value(ops["N_alpha"]).real - 2) < 1e-6
    # 4. 16 CZ en línea, sin SWAP
    assert cz_count(ans.circuit, "line") == 16
    # 5. medición con shots (estimator.py del equipo) consistente con la exacta
    bound = ans.circuit.assign_parameters(opt["x"])
    r = measure_energy(bound, H5, 200_000)
    assert abs(r["energy"] - opt["energy"]) < 4 * r["sigma"], (r, opt["energy"])
    # 6. el plegado global no cambia la unitaria y el ruido sube la energía
    isa = to_line(bound)
    assert Operator(fold_global(isa, 1).remove_final_measurements(inplace=False)).equiv(Operator(isa))
    rn = measure_energy_noisy(isa, H5, 200_000)
    assert rn["energy"] - opt["energy"] > 5e-3, rn          # sesgo de ~13 mHa con nuestro modelo de ruido
    print("OK reducción + compacto + ruido", 1e3 * (opt["energy"] - m["e_fci"]), r, rn)


if __name__ == "__main__":
    test_reduction_and_compact()
