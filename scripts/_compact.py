"""Punto del ansatz compacto (5 qubits, L=8) con θ* guardado en disco.

optimize_compact tiene semilla, pero varios arranques llegan a la misma energía (10⁻¹⁵) y cuál
gana depende del redondeo (cantidad de hilos de BLAS). La energía sin ruido es la misma, pero el
circuito cambia y con él el sesgo por ruido (~0,4 mHa). Guardando θ*, mitigation_study.py y
mitigation_shots.py usan EXACTAMENTE el mismo circuito.
"""
from __future__ import annotations

import os

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE = os.path.join(ROOT, "results", "compact_params")


def compact_point(d: float, layers: int = 8) -> dict:
    from qiskit.quantum_info import Statevector
    from h4vqe.hamiltonian import molecule, qubit_hamiltonian
    from h4vqe.noise import to_line
    from h4vqe.reduction import build_compact_ansatz, optimize_compact, reduced_number_operators, tapering
    m = molecule(d)
    H8 = qubit_hamiltonian(m)
    H = tapering(m, H8).apply(H8)
    ops = reduced_number_operators(m)
    ans = build_compact_ansatz(m, layers=layers)
    path = os.path.join(CACHE, f"d{d:.2f}_L{layers}.npy")
    if os.path.exists(path):
        x = np.load(path)
    else:
        x = optimize_compact(ans, H, ops["S2"])["x"]
        os.makedirs(CACHE, exist_ok=True)
        np.save(path, x)
    bound = ans.circuit.assign_parameters(x)
    return {"m": m, "H": H, "ops": ops, "bound": bound, "isa": to_line(bound),
            "E0": float(Statevector(bound).expectation_value(H).real)}
