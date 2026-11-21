"""Chequeos de cordura rápidos. Correr con:  python -m pytest -q tests  (o python tests/test_smoke.py)"""
import os
import sys

import numpy as np
from qiskit.quantum_info import Statevector

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from h4vqe.ansatz import build_ansatz  # noqa: E402
from h4vqe.estimator import measure_energy  # noqa: E402
from h4vqe.hamiltonian import molecule, qubit_hamiltonian, number_operators  # noqa: E402
from h4vqe.optimize import optimize_params  # noqa: E402

D = 2.0


def test_all():
    m = molecule(D)
    H = qubit_hamiltonian(m)
    ans = build_ansatz(m)

    # 1. HF: el estado |00110011> da la energía HF de PySCF
    e_hf = Statevector(ans.fast.psi0).expectation_value(H).real
    assert abs(e_hf - m["e_hf"]) < 1e-8, (e_hf, m["e_hf"])

    # 2. El simulador numpy coincide con el circuito de Qiskit
    x = np.random.default_rng(0).normal(0, 0.2, ans.num_parameters)
    sv = Statevector(ans.circuit.assign_parameters(x))
    fid = abs(np.vdot(sv.data, ans.fast.state(x))) ** 2
    assert fid > 1 - 1e-9, fid

    # 3. Gradiente analítico == diferencias finitas
    Hm = H.to_matrix(sparse=True)
    e0, g = ans.fast.energy_and_grad(x, Hm)
    k, h = 3, 1e-6
    xp = x.copy(); xp[k] += h
    assert abs((ans.fast.energy_and_grad(xp, Hm)[0] - e0) / h - g[k]) < 1e-4

    # 4. Variacional: VQE >= FCI, y conserva N_alfa = N_beta = 2
    opt = optimize_params(ans, H)
    assert opt["energy"] >= m["e_fci"] - 1e-9
    psi = ans.fast.state(opt["x"])
    ops = number_operators(m["norb"])
    sv = Statevector(psi)
    assert abs(sv.expectation_value(ops["N_alpha"]).real - 2) < 1e-8
    assert abs(sv.expectation_value(ops["N_beta"]).real - 2) < 1e-8

    # 5. Energía medida con shots consistente con la exacta (dentro de 4σ)
    r = measure_energy(ans.circuit.assign_parameters(opt["x"]), H, 100_000)
    assert abs(r["energy"] - opt["energy"]) < 4 * r["sigma"] + 1e-6, (r, opt["energy"])
    print("OK", D, opt["energy"] - m["e_fci"], r)


if __name__ == "__main__":
    test_all()
