"""Chequeos del extra de FakeFez y de la corrección de lectura. Correr con:  python -m pytest -q tests"""
import os
import sys
from functools import reduce

import numpy as np
from qiskit_aer.primitives import EstimatorV2 as AerEstimator

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from h4vqe.hamiltonian import molecule, qubit_hamiltonian  # noqa: E402
from h4vqe.hardware import best_line, fez, line_noise_model, to_fez  # noqa: E402
from h4vqe.noise import readout_corrected_values, to_line  # noqa: E402
from h4vqe.reduction import build_compact_ansatz, tapering  # noqa: E402


def test_readout_correction():
    """Σ_x p_medida(x) v'(x) = Σ_x p_real(x) v(x) con p_medida = (⊗ A_q) p_real."""
    rng = np.random.default_rng(0)
    n = 5
    A = [np.array([[1 - a, b], [a, 1 - b]]) for a, b in rng.uniform(0, 0.1, (n, 2))]
    full = reduce(np.kron, A[::-1])                 # índice x = Σ b_q 2^q
    p = rng.dirichlet(np.ones(2 ** n))
    v = rng.normal(size=2 ** n)
    assert abs(readout_corrected_values(v, A) @ (full @ p) - v @ p) < 1e-12


def test_line_noise_model_equals_full():
    """El modelo reducido a la línea da el mismo valor esperado que el chip completo de 156 qubits."""
    m = molecule(2.4)
    H8 = qubit_hamiltonian(m)
    H5 = tapering(m, H8).apply(H8)
    ans = build_compact_ansatz(m, layers=8)
    x = np.random.default_rng(1).normal(0, 0.3, ans.num_parameters)
    isa = to_line(ans.circuit.assign_parameters(x))
    line, info = best_line(isa)
    assert info["expected_errors"] <= info["median_expected_errors"]

    def ev(c, H, nm):
        est = AerEstimator(options={"backend_options": {"noise_model": nm, "method": "density_matrix"}})
        return float(est.run([(c, H)], precision=0).result()[0].data.evs)

    tc = to_fez(line)([isa])[0]
    e_full = ev(tc, H5.apply_layout(tc.layout), fez()[1])
    e_line = ev(isa, H5, line_noise_model(tuple(line)))
    assert abs(e_full - e_line) < 1e-10, (e_full, e_line)


if __name__ == "__main__":
    test_readout_correction()
    test_line_noise_model_equals_full()
    print("OK")
