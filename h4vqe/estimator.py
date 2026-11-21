"""Energía MEDIDA con SamplerV2 de qiskit_aer, con barra de error ±σ.

Dueño: Santino.

Método:
  1. Se agrupan los términos de Pauli que conmutan qubit a qubit (QWC).
  2. Cada grupo se mide en su base (H para X, S†H para Y) con N_g shots.
  3. Cada shot da UNA muestra del observable del grupo O_g = Σ c_k P_k,
     así que E_g = media muestral y σ_g² = var muestral / N_g
     (las covarianzas dentro del grupo quedan incluidas).
  4. E = constante + Σ_g E_g ;  σ² = Σ_g σ_g²  (grupos independientes).

Interfaz: measure_energy(circuit_bound, H, total_shots, ...) -> dict

Asignación de shots en dos etapas (piloto + Neyman), ver measure_energy.

Implementación: para cada grupo se precalcula una tabla con el valor de O_g para
cada bitstring posible (2^n entradas). Los conteos se hacen con np.bincount
sobre el arreglo de bits del SamplerV2 (mucho más rápido que get_counts con
millones de shots). Los circuitos de medición se transpilan una sola vez
(`prepare_measurement`) y se pueden reusar entre llamadas.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.quantum_info import SparsePauliOp
from qiskit_aer.primitives import SamplerV2 as AerSampler

from .config import BASIS_GATES, SEED


def split_constant(H: SparsePauliOp) -> tuple[float, SparsePauliOp]:
    ident = "I" * H.num_qubits
    const = sum(c for p, c in H.to_list() if p == ident)
    rest = SparsePauliOp.from_list([(p, c) for p, c in H.to_list() if p != ident])
    return float(np.real(const)), rest


def measurement_groups(H: SparsePauliOp) -> list[SparsePauliOp]:
    _, rest = split_constant(H)
    return rest.group_commuting(qubit_wise=True)


def group_basis(group: SparsePauliOp) -> str:
    """Base común del grupo: por qubit X, Y, Z o I (string little-endian como Pauli)."""
    nq = group.num_qubits
    basis = ["I"] * nq
    for p in group.paulis:
        lab = p.to_label()  # lab[-1-q] es el qubit q
        for q in range(nq):
            ch = lab[-1 - q]
            if ch != "I":
                basis[q] = ch
    return "".join(reversed(basis))


def rotate_and_measure(circ: QuantumCircuit, basis: str) -> QuantumCircuit:
    qc = circ.copy()
    nq = qc.num_qubits
    for q in range(nq):
        ch = basis[-1 - q]
        if ch == "X":
            qc.h(q)
        elif ch == "Y":
            qc.sdg(q)
            qc.h(q)
    qc.measure_all()
    return qc


def value_table(group: SparsePauliOp) -> np.ndarray:
    """table[x] = valor de O_g = Σ c_k P_k cuando se mide el bitstring x (bit q = qubit q),
    ya rotado a la base del grupo (cada P_k actúa como Z en su soporte)."""
    nq = group.num_qubits
    xs = np.arange(2 ** nq, dtype=np.uint64)
    table = np.zeros(2 ** nq)
    for p, c in zip(group.paulis, np.real(group.coeffs)):
        mask = np.uint64(int(np.dot(p.z | p.x, 1 << np.arange(nq, dtype=np.int64))))
        parity = np.bitwise_count(xs & mask) & 1
        table += c * (1 - 2 * parity.astype(float))
    return table


@dataclass
class MeasurementPlan:
    """Lo que no depende de los shots: grupos, circuitos transpilados y tablas de valores."""
    const: float
    groups: list
    circuits: list
    tables: list
    num_qubits: int

    @property
    def n_groups(self) -> int:
        return len(self.groups)


def prepare_measurement(circuit: QuantumCircuit, H: SparsePauliOp) -> MeasurementPlan:
    const, _ = split_constant(H)
    groups = measurement_groups(H)
    circs = [rotate_and_measure(circuit, group_basis(g)) for g in groups]
    circs = transpile(circs, basis_gates=BASIS_GATES + ["measure"], optimization_level=1,
                      seed_transpiler=SEED)
    return MeasurementPlan(const, groups, circs, [value_table(g) for g in groups],
                           circuit.num_qubits)


def allocate_shots(groups, total_shots: int, weights=None) -> list[int]:
    """Por defecto: N_g ∝ Σ|c_k| (proxy de la desviación del grupo).
    Con weights = σ_g estimadas -> asignación de Neyman (N_g ∝ σ_g)."""
    if weights is None:
        weights = np.array([np.sum(np.abs(g.coeffs)) for g in groups])
    w = np.asarray(weights, float)
    w = w / w.sum()
    shots = np.maximum(100, np.floor(w * total_shots)).astype(int)
    return shots.tolist()


def _bitarray_to_int(bitarray) -> np.ndarray:
    """BitArray del SamplerV2 -> entero por shot (bit q = qubit q)."""
    arr = bitarray.array  # (shots, nbytes) uint8, big-endian
    x = np.zeros(arr.shape[0], dtype=np.int64)
    for i in range(arr.shape[1]):
        x = (x << 8) | arr[:, i]
    return x


def _sample_counts(plan: MeasurementPlan, shots, seed) -> np.ndarray:
    """Corre un circuito por grupo; devuelve conteos por bitstring, shape (n_groups, 2^n)."""
    # Una semilla INDEPENDIENTE por grupo. Con una sola corrida y seed fija, Aer usa
    # semillas consecutivas por circuito y vimos errores correlacionados entre grupos
    # (σ subestimada, desvíos de 4σ). SeedSequence da flujos independientes.
    ss = seed if isinstance(seed, np.random.SeedSequence) else np.random.SeedSequence(seed)
    seeds = ss.generate_state(plan.n_groups)
    counts = np.zeros((plan.n_groups, 2 ** plan.num_qubits))
    for i, (c, n, sd) in enumerate(zip(plan.circuits, shots, seeds)):
        r = AerSampler(seed=int(sd)).run([(c, None, int(n))]).result()[0]
        counts[i] = np.bincount(_bitarray_to_int(r.data.meas), minlength=counts.shape[1])
    return counts


def _stats(plan: MeasurementPlan, counts: np.ndarray) -> np.ndarray:
    """(n, Σx, Σx²) por grupo a partir de los conteos."""
    T = np.array(plan.tables)
    return np.stack([counts.sum(1), (counts * T).sum(1), (counts * T ** 2).sum(1)], axis=1)


def _combine(const, stats):
    n, s1, s2 = stats[:, 0], stats[:, 1], stats[:, 2]
    mean = s1 / n
    svar = np.maximum(s2 - n * mean ** 2, 0) / (n - 1)
    return const + mean.sum(), np.sqrt(np.sum(svar / n)), np.sqrt(svar)


def bootstrap_sigma(plan: MeasurementPlan, stage_counts: list[np.ndarray], n_boot: int = 200,
                    seed: int = 0) -> float:
    """σ por bootstrap: se re-muestrean los shots de cada grupo y de cada etapa
    (multinomial con las frecuencias observadas, mismo N) y se recalcula E.
    Chequeo cruzado de la σ analítica, sin suponer normalidad."""
    rng = np.random.default_rng(seed)
    T = np.array(plan.tables)
    n_tot = sum(c.sum(1) for c in stage_counts)                    # (n_groups,)
    E_b = np.zeros(n_boot)
    for g in range(plan.n_groups):
        s = np.zeros(n_boot)
        for counts in stage_counts:
            n = int(counts[g].sum())
            if n == 0:
                continue
            p = counts[g] / n
            s += rng.multinomial(n, p / p.sum(), size=n_boot) @ T[g]
        E_b += s / n_tot[g]
    return float(np.std(E_b, ddof=1))


def measure_energy(circuit: QuantumCircuit, H: SparsePauliOp, total_shots: int,
                   seed: int = SEED, pilot_fraction: float = 0.1,
                   return_details: bool = False, plan: MeasurementPlan | None = None,
                   weights=None, n_boot: int = 0) -> dict:
    """circuit: ansatz con parámetros YA asignados, sin mediciones.

    Dos etapas:
      1. piloto (pilot_fraction del presupuesto) con N_g ∝ Σ|c_k|  -> estima σ_g;
      2. resto con asignación de Neyman N_g ∝ σ_g (minimiza σ total para N fijo).
    Las muestras de ambas etapas se juntan (misma distribución).
    pilot_fraction=1 -> una sola etapa con N_g ∝ Σ|c_k|.
    weights=array -> una sola etapa con N_g ∝ weights (p. ej. σ_g teóricas).

    plan: resultado de prepare_measurement(circuit, H) para no re-transpilar
          (si se pasa, circuit y H no se usan).
    n_boot > 0: agrega "sigma_boot" (σ por bootstrap con n_boot réplicas).
    """
    plan = plan or prepare_measurement(circuit, H)
    seed_pilot, seed_rest, seed_boot = np.random.SeedSequence(seed).spawn(3)
    if weights is not None:
        pilot_fraction, w_pilot = 1.0, weights
    else:
        w_pilot = None
    stage_counts = [_sample_counts(plan, allocate_shots(plan.groups, int(total_shots * pilot_fraction),
                                                        weights=w_pilot), seed_pilot)]
    stats = _stats(plan, stage_counts[0])
    rest = total_shots - int(stats[:, 0].sum())
    if rest > 0 and pilot_fraction < 1:
        _, _, sg = _combine(plan.const, stats)
        stage_counts.append(_sample_counts(plan, allocate_shots(plan.groups, rest, weights=sg + 1e-12),
                                           seed_rest))
        stats = stats + _stats(plan, stage_counts[1])
    E, sigma, sg = _combine(plan.const, stats)
    out = {"energy": float(E), "sigma": float(sigma), "shots": int(stats[:, 0].sum()),
           "n_groups": plan.n_groups}
    if n_boot:
        out["sigma_boot"] = bootstrap_sigma(plan, stage_counts, n_boot,
                                            int(seed_boot.generate_state(1)[0]))
    if return_details:
        out["groups"] = [{"n_terms": len(g), "shots": int(st[0]), "mean": st[1] / st[0], "std_sample": s}
                         for g, st, s in zip(plan.groups, stats, sg)]
    return out
