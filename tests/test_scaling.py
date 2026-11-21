"""Chequeos del extra de escalado: para N = 4 todo tiene que coincidir con el desafío central."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
from qiskit.quantum_info import Statevector

from h4vqe.ansatz import excitation_circuit, uccsd_excitations
from h4vqe.estimator import measurement_groups, split_constant
from h4vqe.hamiltonian import molecule, qubit_hamiltonian
from h4vqe.scaling import (chain_molecule, fci_support, group_moments, qwc_sorted_insertion,
                           uccsd_cost)


def test_h4_matches_central():
    d = 1.4
    m = chain_molecule(4, d)
    ref = molecule(d)
    assert abs(m["e_fci"] - ref["e_fci"]) < 1e-8
    assert abs(m["e_ccsd"] - ref["e_ccsd"]) < 1e-6
    H = qubit_hamiltonian(m)
    H_ref = qubit_hamiltonian(ref)
    assert len(H) == len(H_ref) == 185
    _, rest = split_constant(H)
    mine = qwc_sorted_insertion(rest)
    theirs = measurement_groups(H, "sorted", "qwc")
    assert len(mine) == len(theirs)          # 37 con sorted insertion (46 con coloreo)
    assert all(sorted(rest[g].paulis.to_labels()) == sorted(t.paulis.to_labels())
               for g, t in zip(mine, theirs))


def test_group_moments_vs_statevector():
    m = chain_molecule(4, 2.0)
    H = qubit_hamiltonian(m)
    const, rest = split_constant(H)
    groups = qwc_sorted_insertion(rest)
    idx, amps = fci_support(m)
    mean, var = group_moments(rest, groups, idx, amps)
    assert abs(const + mean.sum() - m["e_fci"]) < 1e-8
    psi = np.zeros(2 ** 8, dtype=complex)
    psi[idx] = amps
    sv = Statevector(psi)
    for g, v in zip(groups[:10], var[:10]):
        P = rest[g]
        e1 = sv.expectation_value(P).real
        e2 = sv.expectation_value(P.dot(P).simplify()).real
        assert abs((e2 - e1 ** 2) - v) < 1e-10


def test_uccsd_cost_formula():
    for norb in (4, 6):
        m = {"norb": norb, "nelec": (norb // 2, norb // 2)}
        exc = uccsd_excitations(m)
        cost = uccsd_cost(norb, norb // 2, norb // 2)
        assert cost["params"] == len(exc)
        nq = 2 * norb
        cx = sum(excitation_circuit(e, nq, 0.3).count_ops().get("cx", 0) for e in exc)
        assert cost["cz_cota"] == cx, (norb, cost["cz_cota"], cx)
    assert uccsd_cost(4, 2, 2)["params"] == 26
    assert uccsd_cost(54, 27, 27)["params"] == 779_301


if __name__ == "__main__":
    test_h4_matches_central()
    test_group_moments_vs_statevector()
    test_uccsd_cost_formula()
    print("OK")
