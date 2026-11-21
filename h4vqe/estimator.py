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
"""
from __future__ import annotations

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


def allocate_shots(groups, total_shots: int, weights=None) -> list[int]:
    """Por defecto: N_g ∝ Σ|c_k| (proxy de la desviación del grupo).
    Mejorable con varianzas estimadas (asignación de Neyman): N_g ∝ σ_g."""
    if weights is None:
        weights = np.array([np.sum(np.abs(g.coeffs)) for g in groups])
    w = np.asarray(weights, float)
    w = w / w.sum()
    shots = np.maximum(100, np.floor(w * total_shots)).astype(int)
    return shots.tolist()


def _group_samples(group: SparsePauliOp, counts: dict) -> tuple[np.ndarray, np.ndarray]:
    """Valores de O_g para cada bitstring observado y su multiplicidad."""
    nq = group.num_qubits
    z_masks = []
    for p in group.paulis:
        lab = p.to_label()
        z_masks.append(sum(1 << q for q in range(nq) if lab[-1 - q] != "I"))
    coeffs = np.real(group.coeffs)
    vals, mult = [], []
    for bits, n in counts.items():
        x = int(bits, 2)
        v = sum(c * (1 - 2 * (bin(x & mk).count("1") & 1)) for c, mk in zip(coeffs, z_masks))
        vals.append(v)
        mult.append(n)
    return np.array(vals), np.array(mult)


def _sample_groups(circuit, groups, shots, seed):
    """Corre un circuito por grupo y devuelve, por grupo, (n, Σx, Σx²) de las muestras."""
    circs = [rotate_and_measure(circuit, group_basis(g)) for g in groups]
    circs = transpile(circs, basis_gates=BASIS_GATES + ["measure"], optimization_level=1,
                      seed_transpiler=SEED)
    # Una semilla INDEPENDIENTE por grupo. Con una sola corrida y seed fija, Aer usa
    # semillas consecutivas por circuito y vimos errores correlacionados entre grupos
    # (σ subestimada, desvíos de 4σ). SeedSequence da flujos independientes.
    seeds = np.random.SeedSequence(seed).generate_state(len(circs))
    stats = []
    for g, c, n, sd in zip(groups, circs, shots, seeds):
        r = AerSampler(seed=int(sd)).run([(c, None, int(n))]).result()[0]
        vals, mult = _group_samples(g, r.data.meas.get_counts())
        stats.append(np.array([mult.sum(), np.sum(mult * vals), np.sum(mult * vals ** 2)], float))
    return np.array(stats)


def _combine(const, stats):
    n, s1, s2 = stats[:, 0], stats[:, 1], stats[:, 2]
    mean = s1 / n
    svar = np.maximum(s2 - n * mean ** 2, 0) / (n - 1)
    return const + mean.sum(), np.sqrt(np.sum(svar / n)), np.sqrt(svar)


def measure_energy(circuit: QuantumCircuit, H: SparsePauliOp, total_shots: int,
                   seed: int = SEED, pilot_fraction: float = 0.1,
                   return_details: bool = False) -> dict:
    """circuit: ansatz con parámetros YA asignados, sin mediciones.

    Dos etapas:
      1. piloto (pilot_fraction del presupuesto) con N_g ∝ Σ|c_k|  -> estima σ_g;
      2. resto con asignación de Neyman N_g ∝ σ_g (minimiza σ total para N fijo).
    Las muestras de ambas etapas se juntan (misma distribución).
    pilot_fraction=1 -> una sola etapa con N_g ∝ Σ|c_k|.
    """
    const, _ = split_constant(H)
    groups = measurement_groups(H)
    n_pilot = int(total_shots * pilot_fraction)
    stats = _sample_groups(circuit, groups, allocate_shots(groups, n_pilot), seed)
    rest = total_shots - int(stats[:, 0].sum())
    if rest > 0:
        _, _, sg = _combine(const, stats)
        stats = stats + _sample_groups(circuit, groups,
                                       allocate_shots(groups, rest, weights=sg + 1e-12), seed + 1)
    E, sigma, sg = _combine(const, stats)
    out = {"energy": float(E), "sigma": float(sigma), "shots": int(stats[:, 0].sum()),
           "n_groups": len(groups)}
    if return_details:
        out["groups"] = [{"n_terms": len(g), "shots": int(st[0]), "mean": st[1] / st[0], "std_sample": s}
                         for g, st, s in zip(groups, stats, sg)]
    return out
