"""Punto del ansatz compacto (5 qubits, L=8) con θ* Y EL HAMILTONIANO guardados en disco.

Por qué guardar también H: el signo de los orbitales moleculares de PySCF depende de la cantidad
de hilos de BLAS (a 1,8 y 3,0 Å, con 1–2 hilos sale un gauge y con ≥4 otro). Con otro gauge cambia
el signo de algunos términos de H, y un θ* optimizado en un proceso puede dar cientos de mHa de error
en otro. Guardando (θ*, H₅, N_α, N_β, S²) juntos, mitigation_study.py y mitigation_shots.py usan
exactamente el mismo circuito y el mismo H, en cualquier máquina y con cualquier cantidad de hilos.
(Además, optimize_compact tiene empates entre arranques que se rompen según el redondeo: con otro θ*
de la misma energía, el sesgo por ruido cambia ~0,5 mHa.)
"""
from __future__ import annotations

import os
import pickle

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE = os.path.join(ROOT, "results", "compact_params")


def compact_point(d: float, layers: int = 8) -> dict:
    from qiskit.quantum_info import Statevector
    from h4vqe.hamiltonian import molecule, qubit_hamiltonian
    from h4vqe.noise import to_line
    from h4vqe.reduction import build_compact_ansatz, optimize_compact, reduced_number_operators, tapering
    m = molecule(d)  # solo energías de referencia (invariantes de gauge) y estructura del ansatz
    ans = build_compact_ansatz(m, layers=layers)
    path = os.path.join(CACHE, f"d{d:.2f}_L{layers}.pkl")
    if os.path.exists(path):
        with open(path, "rb") as f:
            c = pickle.load(f)
        x, H, ops = c["x"], c["H"], c["ops"]
    else:
        H8 = qubit_hamiltonian(m)
        H = tapering(m, H8).apply(H8)
        ops = reduced_number_operators(m)
        x = optimize_compact(ans, H, ops["S2"])["x"]
        os.makedirs(CACHE, exist_ok=True)
        with open(path, "wb") as f:
            pickle.dump({"x": x, "H": H, "ops": ops}, f)
    bound = ans.circuit.assign_parameters(x)
    E0 = float(Statevector(bound).expectation_value(H).real)
    assert abs(E0 - m["e_fci"]) < 1e-8, f"θ* guardado no corresponde a H ({1e3*(E0-m['e_fci']):.3f} mHa)"
    return {"m": m, "H": H, "ops": ops, "bound": bound, "isa": to_line(bound), "E0": E0}
