"""Chequeos rápidos de analysis.py (ajuste del mínimo y diagnósticos medidos).

    python tests/test_analysis.py      (o python -m pytest -q tests)
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from h4vqe.analysis import (echo_variance, fit_minimum, measured_symmetries,  # noqa: E402
                            measured_variance, variance)
from h4vqe.ansatz import build_ansatz  # noqa: E402
from h4vqe.config import DISTANCES  # noqa: E402
from h4vqe.estimator import measure_energy  # noqa: E402
from h4vqe.hamiltonian import molecule, number_operators, qubit_hamiltonian  # noqa: E402


def test_fit_recovers_known_minimum():
    # Curva anarmónica tipo Morse con mínimo conocido en la grilla real.
    d = np.array(DISTANCES)
    re, De, a = 0.888, 0.3, 1.9
    E = De * (1 - np.exp(-a * (d - re))) ** 2 - De
    f = fit_minimum(d, E, np.full(len(d), 1e-6), n_mc=50)
    assert abs(f["R_e"] - re) < 2e-3, f
    assert abs(f["E_min"] + De) < 1e-4, f


def test_measured_diagnostics_hf():
    # Estado HF a 2,0 Å: Var grande y conocida, N_α = N_β = 2 exactos.
    m = molecule(2.0)
    H = qubit_hamiltonian(m)
    ans = build_ansatz(m)
    x = np.zeros(ans.num_parameters)
    circ = ans.circuit.assign_parameters(x)
    exact = variance(ans.fast.state(x), H.to_matrix(sparse=True))
    e = measure_energy(circ, H, 200_000)
    v = measured_variance(circ, H, 300_000, e["energy"], e["sigma"])
    assert abs(v["var"] - exact) < 4 * v["var_sigma"], (v, exact)
    ev = echo_variance(circ, H, 50_000, t=0.5)
    assert abs(ev["var_echo"] - exact) < 0.15 * exact, (ev, exact)  # O(t⁴) + shots
    s = measured_symmetries(circ, number_operators(m["norb"]), 20_000)
    assert abs(s["N_alpha"] - 2) < 1e-9 and abs(s["N_beta"] - 2) < 1e-9, s
    print("OK", exact, v["var"], ev["var_echo"], s)


if __name__ == "__main__":
    test_fit_recovers_known_minimum()
    test_measured_diagnostics_hf()
