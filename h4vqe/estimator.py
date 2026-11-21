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


def group_outcome_values(group: SparsePauliOp) -> np.ndarray:
    """Valor de O_g = Σ c_k P_k para cada uno de los 2^n resultados posibles (entero x = bitstring)."""
    nq = group.num_qubits
    x = np.arange(2 ** nq)
    bits = (x[:, None] >> np.arange(nq)) & 1                       # (2^n, n), bit q = qubit q
    masks = np.array([[lab[-1 - q] != "I" for q in range(nq)]
                      for lab in (p.to_label() for p in group.paulis)], dtype=int)  # (T, n)
    parity = (bits @ masks.T) & 1                                  # (2^n, T)
    return (1 - 2 * parity) @ np.real(group.coeffs)


def _outcome_hist(bit_array, nq: int) -> np.ndarray:
    """Histograma de resultados leyendo el BitArray con numpy (get_counts arma 1 string por shot: lento)."""
    a = bit_array.array.astype(np.int64)                            # (shots, n_bytes), big-endian
    ints = np.zeros(a.shape[0], dtype=np.int64)
    for k in range(a.shape[1]):
        ints = (ints << 8) | a[:, k]
    return np.bincount(ints, minlength=2 ** nq)


def _prepare(circuit, groups, seed=SEED):
    circs = [rotate_and_measure(circuit, group_basis(g)) for g in groups]
    return transpile(circs, basis_gates=BASIS_GATES + ["measure"], optimization_level=1,
                     seed_transpiler=seed)


def _sample_groups(circs, nq, shots, seed):
    """Corre un circuito por grupo y devuelve el histograma de resultados de cada uno."""
    # Una semilla INDEPENDIENTE por grupo. Con una sola corrida y seed fija, Aer usa
    # semillas consecutivas por circuito y vimos errores correlacionados entre grupos
    # (σ subestimada, desvíos de 4σ). SeedSequence da flujos independientes.
    seeds = np.random.SeedSequence(seed).generate_state(len(circs))
    hists = []
    for c, n, sd in zip(circs, shots, seeds):
        r = AerSampler(seed=int(sd)).run([(c, None, int(n))]).result()[0]
        hists.append(_outcome_hist(r.data.meas, nq))
    return np.array(hists)


def _combine(const, hists, values):
    """E = const + Σ_g media_g ; σ² = Σ_g var_g / N_g, a partir de los histogramas."""
    n = hists.sum(axis=1)
    mean = np.einsum("gx,gx->g", hists, values) / n
    svar = np.einsum("gx,gx->g", hists, (values - mean[:, None]) ** 2) / (n - 1)
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
    nq = H.num_qubits
    values = np.array([group_outcome_values(g) for g in groups])   # (G, 2^n)
    circs = _prepare(circuit, groups)                               # se transpila una sola vez
    n_pilot = int(total_shots * pilot_fraction)
    hists = _sample_groups(circs, nq, allocate_shots(groups, n_pilot), seed)
    rest = total_shots - int(hists.sum())
    if rest > 0:
        _, _, sg = _combine(const, hists, values)
        hists = hists + _sample_groups(circs, nq, allocate_shots(groups, rest, weights=sg + 1e-12),
                                       seed + 1)
    E, sigma, sg = _combine(const, hists, values)
    out = {"energy": float(E), "sigma": float(sigma), "shots": int(hists.sum()),
           "n_groups": len(groups)}
    if return_details:
        n = hists.sum(axis=1)
        out["groups"] = [{"n_terms": len(g), "shots": int(ng), "mean": float(h @ v / ng), "std_sample": s}
                         for g, h, v, ng, s in zip(groups, hists, values, n, sg)]
    return out
