"""El ansatz ruteado para qubits en línea da EXACTAMENTE la misma energía (dueño: Gonzalo).

Correr con:  python tests/test_routing.py
"""
import os
import sys

import numpy as np
from qiskit.quantum_info import Statevector

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from h4vqe.adapt import ADAPT_FERMIONIC  # noqa: E402
from h4vqe.ansatz import build_ansatz  # noqa: E402
from h4vqe.hamiltonian import molecule, qubit_hamiltonian  # noqa: E402
from h4vqe.routing import LINE_ORDER, hamiltonian_in_order, routed_ansatz  # noqa: E402


def test_routing():
    m = molecule(2.0)
    H = qubit_hamiltonian(m)
    # Con el orden identidad, hamiltonian_in_order da el mismo Ĥ
    diff = (hamiltonian_in_order(m, list(range(8))) - H).simplify(atol=1e-9)
    assert np.allclose(diff.coeffs, 0, atol=1e-9)
    # Mismos θ -> misma energía, midiendo con Ĥ en el orden final de modos
    ans = build_ansatz(m, excitations=ADAPT_FERMIONIC)
    x = np.random.default_rng(0).normal(0, 0.3, ans.num_parameters)
    e_ref = Statevector(ans.fast.state(x)).expectation_value(H).real
    qc, final_order, _ = routed_ansatz(m, ADAPT_FERMIONIC, LINE_ORDER)
    e = Statevector(qc.assign_parameters(x)).expectation_value(
        hamiltonian_in_order(m, final_order)).real
    assert abs(e - e_ref) < 1e-10, (e, e_ref)
    print("OK routing", e - e_ref, final_order)


if __name__ == "__main__":
    test_routing()
