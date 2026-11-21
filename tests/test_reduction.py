"""Chequeos de la reducción a 5 qubits, el ansatz compacto y el ruido.
Correr con:  python tests/test_reduction.py   (o python -m pytest -q tests)"""
import os
import sys

import numpy as np
from qiskit.quantum_info import Operator, Statevector

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from scipy.linalg import expm  # noqa: E402

from h4vqe.ansatz import cz_count, excitation_generator, uccsd_excitations  # noqa: E402
from h4vqe.estimator import measure_energy  # noqa: E402
from h4vqe.hamiltonian import molecule, qubit_hamiltonian  # noqa: E402
from h4vqe.noise import fold_global, measure_energy_noisy, to_line  # noqa: E402
from h4vqe.reduction import (build_compact_ansatz, optimize_compact, reduced_number_operators,  # noqa: E402
                             tapering)


def test_six_qubits():
    """6 qubits (solo paridades): misma física, todo el pool de excitaciones se reduce y anda con
    orbitales rotados. En 5 qubits solo se reducen las excitaciones simétricas ante la inversión."""
    m = molecule(2.0)
    H = qubit_hamiltonian(m)
    tp = tapering(m, H, n_qubits=6)
    H6 = tp.apply(H)
    assert H6.num_qubits == 6
    M = H6.to_matrix()
    assert abs(np.linalg.eigvalsh(M)[0] - m["e_fci"]) < 1e-9
    assert abs(M[tp.hf_index, tp.hf_index].real - m["e_hf"]) < 1e-9
    # 1. el pool de UCCSD/ADAPT: 26/26 en 6 qubits, 14/26 en 5
    pool = uccsd_excitations(m)
    tp5 = tapering(m)
    G8 = [excitation_generator(e, 8) for e in pool]
    G6 = [tp.apply(G) for G in G8]
    n5 = 0
    for G in G8:
        try:
            tp5.apply(G)
            n5 += 1
        except ValueError:
            pass
    assert (len(pool), n5) == (26, 14), (len(pool), n5)
    # 2. un ansatz de excitaciones reducido da la misma energía que en 8 qubits
    theta = np.random.default_rng(3).uniform(-0.5, 0.5, len(pool))
    psi8 = np.zeros(256, dtype=complex)
    psi8[sum(1 << q for q in (0, 1, 4, 5))] = 1
    psi6 = np.zeros(64, dtype=complex)
    psi6[tp.hf_index] = 1
    for t, g8, g6 in zip(theta, G8, G6):
        psi8 = expm(-1j * t * g8.to_matrix()) @ psi8
        psi6 = expm(-1j * t * g6.to_matrix()) @ psi6
    e8 = np.vdot(psi8, H.to_matrix() @ psi8).real
    e6 = np.vdot(psi6, M @ psi6).real
    assert abs(e8 - e6) < 1e-10, (e8, e6)
    # 3. orbitales rotados (enlaces localizados): 6 qubits sigue exacto; 5 no se puede
    R = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
    U = np.block([[R, np.zeros((2, 2))], [np.zeros((2, 2)), R]])
    g = np.einsum("pqrs,pi,qj,rk,sl->ijkl", m["g"], U, U, U, U)
    Hr = qubit_hamiltonian(dict(m, h=U.T @ m["h"] @ U, g=g))
    assert abs(np.linalg.eigvalsh(tapering(m, Hr, n_qubits=6).apply(Hr).to_matrix())[0] - m["e_fci"]) < 1e-9
    try:
        tp5.apply(Hr)
        raise AssertionError("con orbitales localizados no debería poder reducirse a 5 qubits")
    except ValueError:
        pass
    # 4. el compacto de 6 qubits: 25 CZ en línea y el simulador numpy == circuito
    ans = build_compact_ansatz(m, n_qubits=6)
    assert ans.circuit.num_qubits == 6 and cz_count(ans.circuit, "line") == 25
    x = np.random.default_rng(1).uniform(-1, 1, ans.num_parameters)
    fid = abs(np.vdot(Statevector(ans.circuit.assign_parameters(x)).data, ans.fast.state(x))) ** 2
    assert fid > 1 - 1e-9, fid
    print("OK 6 qubits", len(H6), "términos; pool 26/26 (5q: 14/26); E8 == E6:", e8 - e6)


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
    test_six_qubits()
    test_reduction_and_compact()
