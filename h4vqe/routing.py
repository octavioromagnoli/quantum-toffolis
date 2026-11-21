"""Ruteo manual para qubits EN LÍNEA con fSWAP y orden de modos seguido clásicamente.

Dueño: Gonzalo.

Idea (Co-ADAPT, Ramôa et al., arXiv:2601.20681, 2026; redes fSWAP, Kivlichan et al. 2018):
    - El transpiler con SWAP comunes multiplica las CZ ~2,8× en línea (185 → ~520).
    - Un fSWAP entre qubits vecinos intercambia dos MODOS fermiónicos (incluido el signo
      de Jordan-Wigner) y cuesta 2 CZ (es iSWAP salvo fases locales). Un SWAP cuesta 3.
    - Antes de cada excitación juntamos sus modos en qubits contiguos con fSWAP y NO
      deshacemos la permutación: llevamos clásicamente el orden `order[q] = modo en el
      qubit q`. Con los modos contiguos, la excitación fermiónica no necesita CNOT de
      paridad (no hay espectadores en el medio).
    - Al final, el estado está escrito en Jordan-Wigner con el orden final de modos, así
      que la energía se mide con `hamiltonian_in_order(m, order_final)`.

Los θ son los MISMOS que los del ansatz sin rutear (las fSWAP solo reetiquetan modos),
así que se optimiza con FastAnsatz como siempre y se mide este circuito.
"""
from __future__ import annotations

import numpy as np
from qiskit import QuantumCircuit
from qiskit.circuit import ParameterVector
from qiskit.quantum_info import SparsePauliOp

from .ansatz import _multiplexor_weights, excitation_circuit, excitation_generator, hf_occupation
from .hamiltonian import annihilation_ops


def fswap(qc: QuantumCircuit, q: int) -> None:
    """fSWAP entre q y q+1: SWAP seguido de CZ (2 CZ después de transpilar)."""
    qc.swap(q, q + 1)
    qc.cz(q, q + 1)


def hamiltonian_in_order(m: dict, order: list[int]) -> SparsePauliOp:
    """Ĥ en Jordan-Wigner con el espín-orbital order[q] en el qubit q.

    Misma fórmula que hamiltonian.qubit_hamiltonian (espín-orbital = s*norb + p), pero
    con los modos reordenados. Con order = range(2*norb) da el mismo Ĥ.
    """
    norb = m["norb"]
    nq = 2 * norb
    pos = np.argsort(order)                     # pos[modo] = qubit
    a_q = annihilation_ops(nq)
    a = [a_q[pos[so]] for so in range(nq)]
    ad = [op.adjoint() for op in a]
    h, g = m["h"], m["g"]
    hop = [[[ad[s * norb + p].dot(a[s * norb + q]).simplify() for q in range(norb)]
            for p in range(norb)] for s in range(2)]
    terms = [SparsePauliOp("I" * nq, m["constante"])]
    for s in range(2):
        for p in range(norb):
            for q in range(norb):
                terms.append(h[p, q] * hop[s][p][q])
    for s in range(2):
        for t in range(2):
            for p in range(norb):
                for q in range(norb):
                    second = SparsePauliOp.sum([g[p, q, r, u] * hop[t][r][u]
                                                for r in range(norb) for u in range(norb)])
                    terms.append(0.5 * hop[s][p][q].dot(second))
                    if s == t:
                        for u in range(norb):
                            terms.append(-0.5 * g[p, q, q, u] * hop[s][p][u])
    return SparsePauliOp.sum(terms).simplify(atol=1e-10).chop(1e-10)


def line_double_circuit(exc, nq: int, theta) -> QuantumCircuit:
    """Doble fermiónica sobre 4 qubits CONTIGUOS q0<q1<q2<q3 con CNOT casi solo entre vecinos.

    Escalera (q1→q0), (q2→q3), (q2→q1): todas entre vecinos. Lleva |s1⟩, |s2⟩ a dos
    estados que difieren solo en q2. Después, RY sobre q2 controlada por q1 (4 CNOT,
    vecino), q3 (2 CNOT, vecino) y q0 (2 CNOT a distancia 2). Después, la escalera
    inversa. Con modos contiguos no hay espectadores en el medio: el signo de JW es fijo.
    """
    src, dst = exc
    q0, q1, q2, q3 = sorted(src + dst)
    G = excitation_generator(exc, nq, kind="fermionic").to_matrix(sparse=True)
    s1 = sum(1 << q for q in src)
    s2 = sum(1 << q for q in dst)
    sigma = np.real(-1j * G[s2, s1])           # s1 -> cos θ s1 + σ sin θ s2
    ladder = [(q1, q0), (q2, q3), (q2, q1)]
    bits = {q: (s1 >> q) & 1 for q in (q0, q1, q2, q3)}
    for c, t in ladder:
        bits[t] ^= bits[c]
    ctrls = [q1, q3, q0]                       # bit 0 (4 CNOT), bit 1 (2), bit 2 (2)
    neg = [c for c in ctrls if bits[c] == 0]
    sign = -1.0 if bits[q2] == 1 else 1.0      # s1 en |1> del blanco -> girar al revés
    order, w = _multiplexor_weights(3)
    qc = QuantumCircuit(nq)
    for c, t in ladder:
        qc.cx(c, t)
    if neg:
        qc.x(neg)
    for k, wk in enumerate(w):
        qc.ry(2.0 * sign * sigma * wk * theta, q2)
        qc.cx(ctrls[order[k]], q2)
    if neg:
        qc.x(neg)
    for c, t in reversed(ladder):
        qc.cx(c, t)
    return qc


def _gather(qc: QuantumCircuit, order: list[int], modes: list[int]) -> int:
    """Junta `modes` en qubits contiguos con fSWAP vecinos (sin cruzarlos entre sí).

    Elige la ventana que minimiza los movimientos y devuelve la cantidad de fSWAP.
    """
    pos = sorted(order.index(md) for md in modes)
    k = len(pos)
    best = min(range(len(order) - k + 1),
               key=lambda w: sum(abs(p - (w + i)) for i, p in enumerate(pos)))
    n_sw = 0
    # Mover primero los que van a la izquierda (de izquierda a derecha) y después los que
    # van a la derecha (de derecha a izquierda), para no cruzar los modos.
    targets = [best + i for i in range(k)]
    for i in range(k):
        while pos[i] > targets[i]:
            q = pos[i] - 1
            fswap(qc, q)
            order[q], order[q + 1] = order[q + 1], order[q]
            pos[i] -= 1
            n_sw += 1
    for i in reversed(range(k)):
        while pos[i] < targets[i]:
            q = pos[i]
            fswap(qc, q)
            order[q], order[q + 1] = order[q + 1], order[q]
            pos[i] += 1
            n_sw += 1
    return n_sw


def routed_ansatz(m: dict, excitations: list, init_order: list[int] | None = None):
    """Circuito para qubits en línea + orden final de modos.

    Devuelve (circuit, final_order, n_fswap). Medir con hamiltonian_in_order(m, final_order).
    """
    nq = 2 * m["norb"]
    order = list(init_order or range(nq))
    theta = ParameterVector("θ", len(excitations))
    qc = QuantumCircuit(nq)
    for md in hf_occupation(m):
        qc.x(order.index(md))
    n_sw = 0
    for t, (src, dst) in zip(theta, excitations):
        n_sw += _gather(qc, order, list(src) + list(dst))
        pos = {md: order.index(md) for md in list(src) + list(dst)}
        exc = (tuple(pos[x] for x in src), tuple(pos[x] for x in dst))
        block = (line_double_circuit(exc, nq, t) if len(src) == 2
                 else excitation_circuit(exc, nq, t, "fermionic"))
        qc.compose(block, inplace=True)
    return qc, order, n_sw


# Mejor orden inicial encontrado para ADAPT_FERMIONIC (búsqueda sobre 72 órdenes:
# 24 permutaciones espaciales × {bloques, intercalado, bloques espejados}).
# Orden intercalado por orbital espacial: [α3 β3 α1 β1 α2 β2 α0 β0].
LINE_ORDER = [3, 7, 1, 5, 2, 6, 0, 4]


def line_cz_count(circuit: QuantumCircuit, seeds: int = 3) -> int:
    """CZ del circuito ya ruteado, con layout trivial en línea (el mínimo de varias semillas)."""
    from qiskit import transpile
    from qiskit.transpiler import CouplingMap

    from .config import BASIS_GATES
    cm = CouplingMap.from_line(circuit.num_qubits)
    return min(transpile(circuit, basis_gates=BASIS_GATES, coupling_map=cm,
                         initial_layout=list(range(circuit.num_qubits)),
                         optimization_level=3, seed_transpiler=s).count_ops().get("cz", 0)
               for s in range(seeds))
