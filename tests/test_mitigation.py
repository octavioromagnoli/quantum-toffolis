"""Chequeos de mitigation.py (~30 s):  python tests/test_mitigation.py"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from qiskit.quantum_info import Statevector  # noqa: E402

from h4vqe import mitigation as mt  # noqa: E402
from h4vqe.estimator import measurement_groups, split_constant, value_table  # noqa: E402
from h4vqe.hamiltonian import molecule, qubit_hamiltonian  # noqa: E402
from h4vqe.noise import to_line  # noqa: E402
from h4vqe.reduction import build_compact_ansatz, optimize_compact, reduced_number_operators, tapering  # noqa: E402


def test_all():
    m = molecule(2.4)
    H8 = qubit_hamiltonian(m)
    H = tapering(m, H8).apply(H8)
    ops = reduced_number_operators(m)
    ans = build_compact_ansatz(m, layers=8)
    bound = ans.circuit.assign_parameters(optimize_compact(ans, H, ops["S2"], n_random=0)["x"])
    isa = to_line(bound)
    E0 = Statevector(bound).expectation_value(H).real
    G = measurement_groups(H)
    sym = [(ops["N_alpha"], 2), (ops["N_beta"], 2)]
    A = mt.readout_matrices(H.num_qubits, 0.03)

    # 1. Sin ruido: crudo y post-selección dan la energía exacta (y PS no descarta nada)
    P = mt.exact_probabilities(isa, G, p2=0, p1=0)
    for s in (None, sym):
        U, W, keep = mt.mitigation_tables(G, s)
        assert abs(mt.exact_estimate(P, H, U, W)["E"] - E0) < 1e-9
    U, W, keep = mt.mitigation_tables(G, sym)
    assert abs((P * keep).sum(1) - 1).max() < 1e-9

    # 2. Solo error de lectura: la corrección lo saca exacto (con y sin post-selección)
    P = mt.exact_probabilities(isa, G, p2=0, p1=0, p_ro=0.03)
    assert mt.exact_estimate(P, H, *mt.mitigation_tables(G)[:2])["E"] - E0 > 0.05
    for s in (None, sym):
        assert abs(mt.exact_estimate(P, H, *mt.mitigation_tables(G, s, A)[:2])["E"] - E0) < 1e-9

    # 3. Sin mitigación, el estimador de cociente es el de estimator.py (misma E y misma σ)
    rng = np.random.default_rng(0)
    hists = rng.integers(0, 1000, size=(len(G), 2 ** H.num_qubits)).astype(float)
    U, W, _ = mt.mitigation_tables(G)
    E, s, _ = mt.combine_ratio(0.0, hists, U, W)
    V = np.array([value_table(g) for g in G])
    n = hists.sum(1)
    mean = (hists * V).sum(1) / n
    var = (hists * (V - mean[:, None]) ** 2).sum(1) / (n - 1)
    assert abs(E - mean.sum()) < 1e-12 and abs(s - np.sqrt(np.sum(var / n))) < 1e-12

    # 4. Con ruido: PS baja el sesgo y la ZNE (plegando R·C) lo deja en < 0,05 mHa
    Pl = {l: mt.exact_probabilities(isa, G, lam=l) for l in (1, 3, 5)}
    raw = mt.exact_estimate(Pl[1], H, *mt.mitigation_tables(G)[:2])["E"] - E0
    ps = mt.exact_estimate(Pl[1], H, *mt.mitigation_tables(G, sym)[:2])["E"] - E0
    assert 0 < ps < raw
    from h4vqe.noise import richardson
    Es = [mt.exact_estimate(Pl[l], H, *mt.mitigation_tables(G)[:2])["E"] for l in (1, 3, 5)]
    assert abs(richardson((1, 3, 5), Es, [0, 0, 0])[0] - E0) < 5e-5
    print("OK", f"crudo {1e3*raw:+.2f} mHa, PS {1e3*ps:+.2f} mHa, Richardson "
          f"{1e3*(richardson((1, 3, 5), Es, [0, 0, 0])[0] - E0):+.4f} mHa")


if __name__ == "__main__":
    test_all()
