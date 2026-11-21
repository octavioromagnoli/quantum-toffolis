"""Modelo de ruido propio, medición con ruido y herramientas de escalado de ruido (ZNE).

Dueño: Octavio.

La consigna ya no incluye ruido.py (los organizadores aclararon que quedó de una versión anterior y
que no entra). Usamos un modelo PROPIO, simple y controlado, con los parámetros que describía esa
versión: base cz, sx, x, rz; qubits en línea; error despolarizante de 0,1 % por cz y 0,005 % por
compuerta de 1 qubit. Sirve para medir cuánto daño hace cada compuerta y para probar la mitigación.

Lo que aprendimos (detalles en docs/ideas/diferencial.md):
  * El sesgo por ruido escala con TODAS las compuertas: el circuito compacto (16 CZ) tiene ~150
    compuertas de 1 qubit, y a 0,005 % c/u pesan casi tanto como las 16 CZ al 0,1 %.
    => ZNE plegando SOLO las CZ deja 4-5 mHa de sesgo; el plegado GLOBAL C (C† C)^k lo deja en
       0,01-0,05 mHa (Richardson con λ = 1, 3, 5).
  * El reescalado depolarizante con la fidelidad del eco C C† (sin datos clásicos) sobrecorrige
    (queda 1-3 mHa POR DEBAJO de FCI) porque el ruido es local, no global. Combinado con ZNE lineal
    (λ = 1, 3) queda en 0,02-0,07 mHa con menos amplificación de σ que Richardson.
"""
from __future__ import annotations

import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.quantum_info import SparsePauliOp
from qiskit.transpiler import CouplingMap
from qiskit_aer import AerSimulator
from qiskit_aer.noise import NoiseModel, depolarizing_error
from qiskit_aer.primitives import SamplerV2 as AerSampler

from .config import BASIS_GATES, SEED
from .estimator import allocate_shots, group_basis, measurement_groups, split_constant

P_CZ = 1e-3
P_1Q = 5e-5


# Copias locales de funciones internas de estimator.py: así noise.py solo depende de su interfaz
# pública (split_constant, measurement_groups, group_basis, allocate_shots), que no cambia.
def group_outcome_values(group: SparsePauliOp) -> np.ndarray:
    """Valor de O_g = Σ c_k P_k para cada uno de los 2^n resultados posibles (entero x = bitstring)."""
    nq = group.num_qubits
    bits = (np.arange(2 ** nq)[:, None] >> np.arange(nq)) & 1
    masks = np.array([[lab[-1 - q] != "I" for q in range(nq)]
                      for lab in (p.to_label() for p in group.paulis)], dtype=int)
    return (1 - 2 * ((bits @ masks.T) & 1)) @ np.real(group.coeffs)


def _outcome_hist(bit_array, n_bits: int) -> np.ndarray:
    """Histograma de resultados leyendo el BitArray con numpy."""
    a = bit_array.array.astype(np.int64)                            # (shots, n_bytes), big-endian
    ints = np.zeros(a.shape[0], dtype=np.int64)
    for k in range(a.shape[1]):
        ints = (ints << 8) | a[:, k]
    return np.bincount(ints, minlength=2 ** n_bits)


def _combine(const, hists, values):
    """E = const + Σ_g media_g ; σ² = Σ_g var_g / N_g, a partir de los histogramas."""
    n = hists.sum(axis=1)
    mean = np.einsum("gx,gx->g", hists, values) / n
    svar = np.einsum("gx,gx->g", hists, (values - mean[:, None]) ** 2) / (n - 1)
    return const + mean.sum(), np.sqrt(np.sum(svar / n)), np.sqrt(svar)


def noise_model(p2: float = P_CZ, p1: float = P_1Q, noisy_rz: bool = True) -> NoiseModel:
    nm = NoiseModel(basis_gates=BASIS_GATES)
    nm.add_all_qubit_quantum_error(depolarizing_error(p2, 2), ["cz"])
    nm.add_all_qubit_quantum_error(depolarizing_error(p1, 1), ["sx", "x"] + (["rz"] if noisy_rz else []))
    return nm


def noisy_backend(n_qubits: int, p2: float = P_CZ, p1: float = P_1Q, seed: int | None = None,
                  method: str = "density_matrix") -> AerSimulator:
    return AerSimulator(noise_model=noise_model(p2, p1), coupling_map=CouplingMap.from_line(n_qubits),
                        basis_gates=BASIS_GATES, method=method, seed_simulator=seed)


def noisy_sampler(seed: int | None = None, p2: float = P_CZ, p1: float = P_1Q) -> AerSampler:
    """SamplerV2 de Aer con nuestro modelo de ruido (matriz densidad: exacto y rápido hasta ~10 qubits)."""
    return AerSampler(seed=seed, options={"backend_options": {
        "noise_model": noise_model(p2, p1), "method": "density_matrix"}})


def to_line(circuit: QuantumCircuit, optimization_level: int = 1) -> QuantumCircuit:
    """Circuito ISA para qubits en línea, con layout trivial (el ansatz compacto no necesita SWAP)."""
    n = circuit.num_qubits
    return transpile(circuit, basis_gates=BASIS_GATES, coupling_map=CouplingMap.from_line(n),
                     initial_layout=list(range(n)), optimization_level=optimization_level,
                     seed_transpiler=SEED)


def _inverse(isa: QuantumCircuit) -> QuantumCircuit:
    """C† con las MISMAS compuertas que C (nivel 1 fusiona rotaciones). Con nivel 0 sale con ~5x más rz
    (368 contra 72) y el plegado deja de escalar el ruido exactamente por λ: la extrapolación se pasa
    ~4 mHa. Por eso se verifica que los conteos coincidan."""
    inv = to_line(isa.inverse(), optimization_level=1)
    a, b = isa.count_ops(), inv.count_ops()
    assert all(a.get(g, 0) == b.get(g, 0) for g in ("cz", "sx")), (dict(a), dict(b))
    return inv


def fold_global(isa: QuantumCircuit, k: int) -> QuantumCircuit:
    """C (C† C)^k: misma unitaria, factor de ruido λ = 2k + 1 en TODAS las compuertas.
    Las barreras impiden que el transpilador cancele C† C."""
    inv = _inverse(isa)
    out = isa.copy()
    for _ in range(k):
        out.barrier()
        out.compose(inv, inplace=True)
        out.barrier()
        out.compose(isa, inplace=True)
    return out


def _measurement_circuits(isa: QuantumCircuit, groups) -> list[QuantumCircuit]:
    """Agrega la rotación de base y la medición SIN reoptimizar (no tocar el plegado)."""
    n = isa.num_qubits
    circs = []
    for g in groups:
        b = group_basis(g)
        rot = QuantumCircuit(n)
        for q in range(n):
            ch = b[-1 - q]
            if ch == "X":
                rot.h(q)
            elif ch == "Y":
                rot.sdg(q)
                rot.h(q)
        c = isa.copy()
        c.barrier()
        c.compose(transpile(rot, basis_gates=BASIS_GATES, optimization_level=1), inplace=True)
        c.measure_all()
        circs.append(c)
    return circs


def _sample(circs, shots, seed, make_sampler):
    seeds = np.random.SeedSequence(seed).generate_state(len(circs))
    hists = []
    for c, ns, sd in zip(circs, shots, seeds):
        r = make_sampler(int(sd)).run([(c, None, int(ns))]).result()[0]
        hists.append(_outcome_hist(r.data.meas, c.num_clbits))
    return np.array(hists)


def measure_energy_noisy(isa: QuantumCircuit, H: SparsePauliOp, total_shots: int, seed: int = SEED,
                         pilot_fraction: float = 0.1, p2: float = P_CZ, p1: float = P_1Q,
                         make_sampler=None, to_device=None, readout=None) -> dict:
    """Igual que estimator.measure_energy (QWC + piloto + Neyman, σ de las muestras) pero con ruido.
    isa: circuito YA en compuertas nativas sobre la línea, con parámetros asignados (p. ej. fold_global).
    make_sampler(seed) -> SamplerV2 (por defecto, nuestro modelo con p2, p1); to_device(circs) ->
    circuitos ubicados en un dispositivo (p. ej. hardware.py para FakeFez).
    readout: matrices de confusión por qubit (readout_calibration) para corregir el error de lectura."""
    make_sampler = make_sampler or (lambda sd: noisy_sampler(sd, p2, p1))
    const, _ = split_constant(H)
    groups = measurement_groups(H)
    values = np.array([group_outcome_values(g) for g in groups])
    if readout is not None:
        values = np.array([readout_corrected_values(v, readout) for v in values])
    circs = _measurement_circuits(isa, groups)
    if to_device is not None:
        circs = to_device(circs)
    hists = _sample(circs, allocate_shots(groups, int(total_shots * pilot_fraction)), seed, make_sampler)
    rest = total_shots - int(hists.sum())
    if rest > 0:
        _, _, sg = _combine(const, hists, values)
        hists = hists + _sample(circs, allocate_shots(groups, rest, weights=sg + 1e-12), seed + 1,
                                make_sampler)
    E, sigma, _ = _combine(const, hists, values)
    return {"energy": float(E), "sigma": float(sigma), "shots": int(hists.sum()), "n_groups": len(groups),
            "cz": int(isa.count_ops().get("cz", 0))}


def readout_calibration(n: int, shots: int = 200_000, seed: int = SEED, make_sampler=None,
                        to_device=None, p2: float = P_CZ, p1: float = P_1Q) -> list[np.ndarray]:
    """Matriz de confusión de lectura de cada qubit, A_q[i, j] = P(leer i | preparado j), midiendo
    todo en |0...0> y todo en |1...1> (modelo de lectura independiente por qubit)."""
    make_sampler = make_sampler or (lambda sd: noisy_sampler(sd, p2, p1))
    circs = []
    for prep_one in (False, True):
        c = QuantumCircuit(n)
        if prep_one:
            c.x(range(n))
        c.measure_all()
        circs.append(c)
    if to_device is not None:
        circs = to_device(circs)
    h0, h1 = _sample(circs, [shots, shots], seed, make_sampler)
    bits = (np.arange(2 ** n)[:, None] >> np.arange(n)) & 1          # bit q de cada resultado
    p1_given0 = (h0 @ bits) / shots                                   # P(leer 1 | preparado 0), por qubit
    p0_given1 = (h1 @ (1 - bits)) / shots                             # P(leer 0 | preparado 1), por qubit
    return [np.array([[1 - a, b], [a, 1 - b]]) for a, b in zip(p1_given0, p0_given1)]


def readout_corrected_values(values: np.ndarray, readout: list[np.ndarray]) -> np.ndarray:
    """v' = (A⁻¹)ᵀ v con A = ⊗_q A_q: así Σ_x p_medida(x) v'(x) = Σ_x p_real(x) v(x).
    Corregir los valores (y no los histogramas) deja la σ de las muestras como está."""
    n = len(readout)
    t = values.reshape((2,) * n)                                      # eje k = bit n-1-k
    for q, A in enumerate(readout):
        t = np.moveaxis(np.tensordot(np.linalg.inv(A).T, t, axes=([1], [n - 1 - q])), 0, n - 1 - q)
    return t.reshape(-1)


def echo_fidelity(isa: QuantumCircuit, shots: int = 200_000, seed: int = SEED,
                  p2: float = P_CZ, p1: float = P_1Q) -> dict:
    """Fidelidad por pasada estimada con el eco C C†|0> (ideal: |0...0> con probabilidad 1).
    Modelo depolarizante global: P(0) = F^2 + (1 - F^2)/2^n."""
    n = isa.num_qubits
    eco = isa.copy()
    eco.barrier()
    eco.compose(_inverse(isa), inplace=True)
    eco.measure_all()
    r = noisy_sampler(seed, p2, p1).run([(eco, None, shots)]).result()[0]
    p0 = _outcome_hist(r.data.meas, n)[0] / shots
    d = 2 ** n
    f2 = (p0 - 1 / d) / (1 - 1 / d)
    sp0 = np.sqrt(p0 * (1 - p0) / shots)
    F = float(np.sqrt(max(f2, 1e-12)))
    return {"F": F, "F_sigma": float(sp0 / (1 - 1 / d) / (2 * F)), "p0": float(p0), "shots": shots}


def richardson(lams, energies, sigmas) -> tuple[float, float]:
    """Extrapolación polinomial exacta a λ = 0 por los puntos dados; σ propagada (shots independientes)."""
    lams = np.asarray(lams, float)
    gam = np.array([np.prod([lj / (lj - li) for j, lj in enumerate(lams) if j != i])
                    for i, li in enumerate(lams)])
    return float(gam @ np.asarray(energies)), float(np.sqrt(np.sum((gam * np.asarray(sigmas)) ** 2)))


def identity_coefficient(H: SparsePauliOp) -> float:
    """Tr(H)/2^n = energía del estado totalmente mezclado (no requiere ningún cálculo clásico)."""
    return split_constant(H)[0]


def zne_shot_split(lams, total_shots: int) -> list[int]:
    """Reparto óptimo de shots entre niveles para la extrapolación: N_i ∝ |γ_i| (σ por shot ≈ igual).
    Con λ = (1, 3): 75 % / 25 %. Costo total ≈ (Σ|γ_i|)^2 veces el de una medición sin mitigar
    (4x para λ = 1,3; 12x para λ = 1,3,5)."""
    lams = np.asarray(lams, float)
    gam = np.abs([np.prod([lj / (lj - li) for j, lj in enumerate(lams) if j != i]) for i, li in enumerate(lams)])
    return [int(total_shots * g / gam.sum()) for g in gam]


def zne(isa: QuantumCircuit, H: SparsePauliOp, shots, lams=(1, 3, 5),
        rescale: bool = False, echo_shots: int = 2_000_000, seed: int = SEED, n_mc: int = 4000,
        **meas_kw) -> dict:
    """ZNE con plegado global. shots: total (se reparte con zne_shot_split) o lista por nivel.
    rescale=True: primero reescala cada nivel con la fidelidad del eco,
    E_λ' = (E_λ - (1 - F^λ) E_mix) / F^λ, y después extrapola (con λ = (1, 3) alcanza).
    meas_kw: se pasa a measure_energy_noisy (p. ej. make_sampler/to_device para otro dispositivo;
    el eco usa siempre nuestro modelo, así que rescale=True solo vale con el modelo propio)."""
    per_level = zne_shot_split(lams, shots) if np.isscalar(shots) else list(shots)
    meas = [measure_energy_noisy(fold_global(isa, (l - 1) // 2), H, n, seed=seed + 7 * l, **meas_kw)
            for l, n in zip(lams, per_level)]
    E = np.array([r["energy"] for r in meas])
    S = np.array([r["sigma"] for r in meas])
    out = {"lams": list(lams), "E_levels": E.tolist(), "sigma_levels": S.tolist(),
           "shots": int(sum(r["shots"] for r in meas)), "cz_levels": [r["cz"] for r in meas]}
    if not rescale:
        out["E"], out["sigma"] = richardson(lams, E, S)
        return out
    echo = echo_fidelity(isa, echo_shots, seed=seed + 99)
    emix = identity_coefficient(H)
    lam = np.asarray(lams, float)
    resc = lambda e, F: (e - (1 - F ** lam) * emix) / F ** lam
    out["E"], _ = richardson(lams, resc(E, echo["F"]), S)
    rng = np.random.default_rng(seed)                     # σ por Monte Carlo (E_λ y F con su ruido)
    mc = [richardson(lams, resc(E + rng.normal(0, S), echo["F"] + rng.normal(0, echo["F_sigma"])), S)[0]
          for _ in range(n_mc)]
    out["sigma"] = float(np.std(mc))
    out["F"], out["shots"] = echo["F"], out["shots"] + echo_shots
    return out
