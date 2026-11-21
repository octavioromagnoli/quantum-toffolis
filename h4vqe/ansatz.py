"""Ansatz que conserva el número de partículas (UCCSD con Trotter, 1 paso por defecto).

Dueño: Gonzalo.

Cada excitación k es una compuerta exp(-i θ_k G_k) con G_k = i(T_k - T_k†).
Los términos de Pauli de un mismo G_k conmutan entre sí, así que el circuito
implementa EXACTAMENTE Π_k exp(-i θ_k G_k). Eso permite optimizar con un
simulador numpy rápido (`FastAnsatz`) y después medir el MISMO circuito con shots.

Interfaz que usan los otros módulos (no romperla):
    build_ansatz(m, ...) -> Ansatz   con .circuit (parametrizado), .fast, .excitations
    cz_count(circuit, coupling="all"|"line") -> int
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.circuit import ParameterVector
from qiskit.circuit.library import PauliEvolutionGate
from qiskit.quantum_info import SparsePauliOp
from qiskit.transpiler import CouplingMap

from .config import BASIS_GATES, SEED
from .hamiltonian import annihilation_ops


def hf_occupation(m: dict) -> list[int]:
    norb = m["norb"]
    na, nb = m["nelec"]
    return list(range(na)) + [norb + i for i in range(nb)]


def uccsd_excitations(m: dict) -> list[tuple[tuple[int, ...], tuple[int, ...]]]:
    """Simples y dobles ocupado -> vacío que conservan S_z."""
    norb = m["norb"]
    na, nb = m["nelec"]
    occ = [list(range(na)), [norb + i for i in range(nb)]]
    vac = [list(range(na, norb)), [norb + i for i in range(nb, norb)]]
    ex = [((i,), (a,)) for s in range(2) for i in occ[s] for a in vac[s]]
    ex += [((i, j), (a, b)) for i in occ[0] for j in occ[1] for a in vac[0] for b in vac[1]]
    for s in range(2):
        for x, i in enumerate(occ[s]):
            for j in occ[s][x + 1:]:
                for y, a in enumerate(vac[s]):
                    for b in vac[s][y + 1:]:
                        ex.append(((i, j), (a, b)))
    return ex


def qubit_annihilation_ops(nq: int) -> list[SparsePauliOp]:
    """Q_p = (X_p + iY_p)/2, SIN cadena de Z (excitaciones de qubit, QEB)."""
    return [SparsePauliOp.from_sparse_list([("X", [p], 0.5), ("Y", [p], 0.5j)], nq)
            for p in range(nq)]


def excitation_generator(exc, nq: int, a=None, kind: str = "fermionic") -> SparsePauliOp:
    """G = i(T - T†), con T = a†_a a†_b a_j a_i.

    kind="fermionic": a_p de Jordan-Wigner (con cadenas Z).
    kind="qeb": excitación de qubit, Q_p sin cadenas Z (Yordanov et al., 2020/2021).
    Conserva N_alfa y N_beta igual, pero toca solo 2 o 4 qubits, así que
    cuesta menos CZ, sobre todo con qubits en línea.
    """
    if a is None:
        a = annihilation_ops(nq) if kind == "fermionic" else qubit_annihilation_ops(nq)
    src, dst = exc
    T = SparsePauliOp("I" * nq)
    for p in dst:
        T = T.dot(a[p].adjoint())
    for p in reversed(src):
        T = T.dot(a[p])
    return (1j * (T - T.adjoint())).simplify(atol=1e-12)


def _multiplexor_weights(n_ctrl: int) -> tuple[list[int], np.ndarray]:
    """RY multicontrolada (todos los controles en 1) con 2^n RY y 2^n CNOT (código Gray).

    Devuelve el orden de controles de las CNOT y los pesos w_k: RY(w_k φ) antes de la
    CNOT k. Con X^m antes de una RY(α), el ángulo efectivo es (−1)^m α, así que el
    ángulo total para los controles c es Σ_k (−1)^{m_k(c)} w_k φ. Se pide φ si c = 1…1
    y 0 si no (sistema lineal de ±1, invertible).
    """
    n = 2 ** n_ctrl
    order = [(k & -k).bit_length() - 1 if k < n else n_ctrl - 1 for k in range(1, n + 1)]
    M = np.zeros((n, n))
    for c in range(n):
        m = 0
        for k in range(n):
            M[c, k] = (-1) ** m
            m ^= (c >> order[k]) & 1
    f = np.zeros(n)
    f[n - 1] = 1.0
    return order, np.linalg.solve(M, f)


def excitation_circuit(exc, nq: int, theta, kind: str = "fermionic") -> QuantumCircuit:
    """Circuito compacto que implementa EXACTAMENTE exp(−iθG) para una excitación.

    Simple ((i,), (a,)):     CNOT(a→i), RY(2σθ) en a controlada por i, CNOT(a→i).
    Doble ((i, j), (a, b)):  escalera CNOT(i→j), CNOT(a→b), CNOT(a→i) que lleva el par
    |s1⟩ = |ij ocupados⟩, |s2⟩ = |ab ocupados⟩ a dos estados que difieren solo en a;
    después, RY(2σθ) en a controlada por i = 1, j = 0, b = 0 (8 CNOT); después, la
    escalera inversa. Total: 14 CNOT para la doble y 4 para la simple, actuando SOLO
    sobre los qubits de la excitación, lo que se rutea mucho mejor en línea que
    PauliEvolutionGate.
    Fermiónica: el signo σ de T|s1⟩ depende de la paridad de algunos qubits espectadores
    (cadenas Z de Jordan-Wigner). Se implementa con CNOT(q→a) antes y después de la
    rotación (X en el blanco invierte el ángulo de la RY).
    """
    src, dst = exc
    G = excitation_generator(exc, nq, kind=kind).to_matrix(sparse=True)
    occ = list(src)
    s1 = sum(1 << q for q in src)
    s2 = sum(1 << q for q in dst)
    spect = [q for q in range(nq) if q not in src + dst]
    # -iθG = θ(T − T†): s1 -> cos θ s1 + σ sin θ s2, con σ = <s2|(-iG)|s1> = <s2|T|s1>
    sigma0 = np.real(-1j * G[s2, s1])
    parity = [q for q in spect
              if np.real(-1j * G[s2 | (1 << q), s1 | (1 << q)]) * sigma0 < 0]
    if len(src) == 1:
        (i,), (a,) = src, dst
        ladder = [(a, i)]
        ctrls, vals, target = [i], [1], a
    else:
        (i, j), (a, b) = src, dst
        ladder = [(i, j), (a, b), (a, i)]
        ctrls, vals, target = [i, j, b], [1, 0, 0], a
    order, w = _multiplexor_weights(len(ctrls))
    qc = QuantumCircuit(nq)
    for c, t in ladder:
        qc.cx(c, t)
    neg = [c for c, v in zip(ctrls, vals) if v == 0]
    for q in parity:
        qc.cx(q, target)
    if neg:
        qc.x(neg)
    for k, wk in enumerate(w):
        qc.ry(2.0 * sigma0 * wk * theta, target)
        qc.cx(ctrls[order[k]], target)
    if neg:
        qc.x(neg)
    for q in parity:
        qc.cx(q, target)
    for c, t in reversed(ladder):
        qc.cx(c, t)
    return qc


class FastAnsatz:
    """Simulador numpy exacto de Π_k exp(-i θ_k G_k)|HF> (solo para OPTIMIZAR)."""

    def __init__(self, generators: list[SparsePauliOp], hf_index: int, nq: int):
        self.G = [g.to_matrix(sparse=False) for g in generators]
        self.G2 = [G @ G for G in self.G]
        for G, G2 in zip(self.G, self.G2):  # exp cerrada requiere G^3 = G
            assert np.allclose(G2 @ G, G, atol=1e-10), "G^3 != G"
        self.psi0 = np.zeros(2 ** nq, dtype=complex)
        self.psi0[hf_index] = 1.0

    @classmethod
    def from_matrices(cls, G: list, G2: list, psi0: np.ndarray) -> "FastAnsatz":
        """Igual que el constructor, pero con matrices ya calculadas (lo usa ADAPT)."""
        obj = cls.__new__(cls)
        obj.G, obj.G2, obj.psi0 = list(G), list(G2), psi0
        return obj

    def state(self, theta) -> np.ndarray:
        psi = self.psi0
        for t, G, G2 in zip(theta, self.G, self.G2):
            # exp(-i t G) = I + (cos t - 1) G^2 - i sin t G
            psi = psi + (np.cos(t) - 1.0) * (G2 @ psi) - 1j * np.sin(t) * (G @ psi)
        return psi

    def energy_and_grad(self, theta, Hm: np.ndarray):
        """Energía y gradiente analítico (barrido hacia atrás)."""
        n = len(theta)
        states = [self.psi0]
        for t, G, G2 in zip(theta, self.G, self.G2):
            p = states[-1]
            states.append(p + (np.cos(t) - 1.0) * (G2 @ p) - 1j * np.sin(t) * (G @ p))
        psi = states[-1]
        lam = Hm @ psi
        e = float(np.real(np.vdot(psi, lam)))
        grad = np.zeros(n)
        for k in reversed(range(n)):
            # dE/dθ_k = 2 Re <λ| (-i G_k) |ψ_k>, con ψ_k el estado después de la compuerta k
            grad[k] = 2.0 * np.real(np.vdot(lam, -1j * (self.G[k] @ states[k + 1])))
            t, G, G2 = theta[k], self.G[k], self.G2[k]
            # λ <- U_k† λ
            lam = lam + (np.cos(t) - 1.0) * (G2 @ lam) + 1j * np.sin(t) * (G @ lam)
        return e, grad


@dataclass
class Ansatz:
    circuit: QuantumCircuit          # parametrizado, en compuertas nativas (sin acoplamiento)
    excitations: list
    fast: FastAnsatz
    params: ParameterVector
    meta: dict = field(default_factory=dict)

    @property
    def num_parameters(self) -> int:
        return len(self.params)


def build_ansatz(m: dict, excitations=None, reps: int = 1, kind: str = "fermionic") -> Ansatz:
    """UCCSD (o la lista de excitaciones que se pase, p. ej. la de ADAPT).

    kind="fermionic" (default) o "qeb" (excitaciones de qubit, ver excitation_generator).
    """
    norb = m["norb"]
    nq = 2 * norb
    excitations = list(excitations or uccsd_excitations(m)) * reps
    a = annihilation_ops(nq) if kind == "fermionic" else qubit_annihilation_ops(nq)
    gens = [excitation_generator(e, nq, a, kind) for e in excitations]
    theta = ParameterVector("θ", len(gens))
    qc = QuantumCircuit(nq)
    occ = hf_occupation(m)
    for q in occ:
        qc.x(q)
    for t, G in zip(theta, gens):
        qc.append(PauliEvolutionGate(G, time=t), range(nq))
    native = transpile(qc, basis_gates=BASIS_GATES, optimization_level=1, seed_transpiler=SEED)
    hf_index = sum(1 << q for q in occ)
    return Ansatz(native, excitations, FastAnsatz(gens, hf_index, nq), theta,
                  {"reps": reps, "n_qubits": nq, "kind": kind})


def cz_count(circuit: QuantumCircuit, coupling: str = "all", seed: int = SEED) -> int:
    """CZ después de transpilar. coupling='line' = qubits en línea (modelo de ruido)."""
    kw = {}
    if coupling == "line":
        kw["coupling_map"] = CouplingMap.from_line(circuit.num_qubits)
    t = transpile(circuit, basis_gates=BASIS_GATES, optimization_level=3,
                  seed_transpiler=seed, **kw)
    return int(t.count_ops().get("cz", 0))
