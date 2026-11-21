"""Mitigación de errores sobre el modelo de ruido de noise.py (Octavio).

Dueño: Santino.

Qué agrega (noise.py ya tiene el plegado global, el eco, Richardson y la calibración de lectura):
  * Post-selección por simetría: en los grupos de medición que miden en Z todos los qubits del
    soporte de N_α (o N_β), cada shot da N_α. Los shots con N_α ≠ 2 (o N_β ≠ 2) son errores seguros
    y se descartan. Funciona en 8 qubits (JW) y en 5 qubits (reducción de Octavio: N_α y N_β
    siguen siendo diagonales). En los demás grupos no se puede chequear y no se descarta nada.
  * Corrección de lectura COMBINABLE con la post-selección, y error de lectura en el modelo
    (el modelo de noise.py no tiene, así que ahí la corrección no hace nada).
  * ZNE con las dos cosas adentro, y evaluación EXACTA (sin shots, matriz densidad) del sesgo y
    del costo en shots de cada combinación, para elegir la técnica antes de gastar shots.
  * ZNE plegando TAMBIÉN la rotación de base de cada grupo: (R C)((R C)†(R C))^k. noise.zne
    pliega solo C y agrega R después: el ruido de R (2-3 compuertas de 1 qubit por qubit) no escala
    con λ y deja un sesgo que la extrapolación no saca. A 1,8 Å (16 CZ): Richardson (1,3,5)
    +0,163 mHa plegando solo C contra +0,0009 mHa plegando R·C (scripts/mitigation_study.py).

Estimador de cociente por grupo. Para cada resultado x del grupo g hay dos tablas:
    u(x) = keep(x)·v(x),  w(x) = keep(x)        (v = valor de O_g, keep = cumple las simetrías)
    E_g = media(u) / media(w)   = media de v sobre los shots que se quedan.
Con corrección de lectura se reemplazan u y w por A⁻ᵀu y A⁻ᵀw (A = ⊗ matrices de confusión):
así Σ p_medida·(A⁻ᵀu) = Σ p_real·u, y el cociente es la post-selección sobre la distribución
corregida. Sin post-selección ni lectura, u = v y w = 1: el estimador de estimator.py.
σ por el método delta: Var(E_g) = Var(u − E_g·w) / (N_g · media(w)²), con la varianza de las
muestras (sin suponer nada del ruido).
"""
from __future__ import annotations

import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.quantum_info import SparsePauliOp
from qiskit_aer import AerSimulator
from qiskit_aer.noise import ReadoutError
from qiskit_aer.primitives import SamplerV2 as AerSampler

from .config import BASIS_GATES, SEED
from .estimator import allocate_shots, group_basis, measurement_groups, split_constant, value_table
from .noise import (P_1Q, P_CZ, _measurement_circuits, _sample, fold_global, noise_model, richardson,
                    zne_shot_split)


# --------------------------------------------------------------------------- modelo con lectura
def noise_model_ro(p2: float = P_CZ, p1: float = P_1Q, p_ro: float = 0.0):
    """El modelo de noise.py más un error de lectura simétrico p_ro por qubit (bit-flip al medir)."""
    nm = noise_model(p2, p1)
    if p_ro > 0:
        nm.add_all_qubit_readout_error(ReadoutError([[1 - p_ro, p_ro], [p_ro, 1 - p_ro]]))
    return nm


def sampler(seed=None, p2: float = P_CZ, p1: float = P_1Q, p_ro: float = 0.0) -> AerSampler:
    return AerSampler(seed=seed, options={"backend_options": {
        "noise_model": noise_model_ro(p2, p1, p_ro), "method": "density_matrix"}})


def readout_matrices(n: int, p_ro: float) -> list[np.ndarray]:
    """Matrices de confusión EXACTAS del modelo (para la evaluación sin shots). Con shots se usa
    noise.readout_calibration(n, make_sampler=lambda s: sampler(s, p_ro=...))."""
    A = np.array([[1 - p_ro, p_ro], [p_ro, 1 - p_ro]])
    return [A] * n


def apply_kron(vec: np.ndarray, mats: list[np.ndarray]) -> np.ndarray:
    """(⊗_q M_q) vec, con vec indexado por bitstring (bit q = qubit q, mats[q] actúa en el qubit q)."""
    n = len(mats)
    t = np.asarray(vec, float).reshape((2,) * n)  # eje k = bit n-1-k
    for q, M in enumerate(mats):
        t = np.moveaxis(np.tensordot(M, t, axes=([1], [n - 1 - q])), 0, n - 1 - q)
    return t.reshape(-1)


def measurement_circuits(isa: QuantumCircuit, groups, lam: int = 1, fold_rotations: bool = True):
    """Circuito de medición de cada grupo con factor de ruido λ (impar).
    fold_rotations=True: se pliega R·C entero (C = ansatz, R = rotación a la base del grupo).
    fold_rotations=False: como noise.zne (pliega solo C y agrega R después; deja sesgo)."""
    k = (lam - 1) // 2
    if not fold_rotations:
        return _measurement_circuits(fold_global(isa, k) if k else isa, groups)
    n = isa.num_qubits
    out = []
    for g in groups:
        b = group_basis(g)
        rot = QuantumCircuit(n)
        for q in range(n):
            if b[-1 - q] == "X":
                rot.h(q)
            elif b[-1 - q] == "Y":
                rot.sdg(q)
                rot.h(q)
        c = isa.copy()
        c.barrier()
        c.compose(transpile(rot, basis_gates=BASIS_GATES, optimization_level=1), inplace=True)
        if k:
            c = fold_global(c, k)
        c.measure_all()
        out.append(c)
    return out


# --------------------------------------------------------------------------- simetrías
def symmetry_masks(groups: list[SparsePauliOp], symmetries) -> np.ndarray:
    """keep[g, x] = True si el resultado x del grupo g es compatible con las simetrías que ese
    grupo puede chequear. symmetries: lista de (operador diagonal, autovalor esperado), p. ej.
    [(N_alpha, 2), (N_beta, 2)]. Un grupo chequea S si mide en Z todos los qubits del soporte de S
    (los qubits con I en la base del grupo no se rotan: también se miden en Z)."""
    n = groups[0].num_qubits
    keep = np.ones((len(groups), 2 ** n), bool)
    for S, target in symmetries or []:
        S = S.simplify()
        assert not np.any(S.paulis.x), "la simetría tiene que ser diagonal (solo I/Z)"
        support = np.flatnonzero(np.any(S.paulis.z, axis=0))
        ok = np.abs(value_table(S) - target) < 1e-6
        for i, g in enumerate(groups):
            b = group_basis(g)
            if all(b[-1 - q] in "ZI" for q in support):
                keep[i] &= ok
    return keep


def mitigation_tables(groups, symmetries=None, readout=None) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Tablas (U, W) del estimador de cociente y la máscara keep, shape (G, 2^n)."""
    V = np.array([value_table(g) for g in groups])
    keep = symmetry_masks(groups, symmetries)
    U, W = keep * V, keep.astype(float)
    if readout is not None:
        Ainv_T = [np.linalg.inv(A).T for A in readout]
        U = np.array([apply_kron(u, Ainv_T) for u in U])
        W = np.array([apply_kron(w, Ainv_T) for w in W])
    return U, W, keep


def combine_ratio(const: float, hists: np.ndarray, U: np.ndarray, W: np.ndarray):
    """E = const + Σ_g media(u)/media(w); σ por el método delta. hists pueden ser conteos o
    probabilidades exactas (entonces σ_g es la desviación POR SHOT). Devuelve (E, σ, σ_g por shot)."""
    hists = np.asarray(hists, float)
    n = hists.sum(1)
    mw = (hists * W).sum(1) / n
    R = (hists * U).sum(1) / n / mw
    resid = U - R[:, None] * W
    var = (hists * resid ** 2).sum(1) / np.maximum(n - 1, 1e-300)
    sg_shot = np.sqrt(var) / np.abs(mw)
    return const + R.sum(), float(np.sqrt(np.sum(sg_shot ** 2 / n))), sg_shot


# --------------------------------------------------------------------------- con shots
def measure_energy_mitigated(isa: QuantumCircuit, H: SparsePauliOp, total_shots: int,
                             seed: int = SEED, pilot_fraction: float = 0.1, symmetries=None,
                             readout=None, p2: float = P_CZ, p1: float = P_1Q, p_ro: float = 0.0,
                             make_sampler=None, to_device=None, lam: int = 1,
                             fold_rotations: bool = True) -> dict:
    """Como noise.measure_energy_noisy (QWC + piloto + Neyman, σ de las muestras) pero con
    post-selección por simetría y/o corrección de lectura. isa: circuito nativo en la línea,
    con layout trivial y parámetros asignados (noise.to_line). lam: factor de plegado (ZNE)."""
    make_sampler = make_sampler or (lambda sd: sampler(sd, p2, p1, p_ro))
    const, _ = split_constant(H)
    groups = measurement_groups(H)
    U, W, keep = mitigation_tables(groups, symmetries, readout)
    circs = measurement_circuits(isa, groups, lam, fold_rotations)
    if to_device is not None:
        circs = to_device(circs)
    s1, s2 = (int(s.generate_state(1)[0]) for s in np.random.SeedSequence(seed).spawn(2))
    hists = _sample(circs, allocate_shots(groups, int(total_shots * pilot_fraction)), s1, make_sampler)
    rest = total_shots - int(hists.sum())
    if rest > 0 and pilot_fraction < 1:
        _, _, sg = combine_ratio(const, hists, U, W)
        hists = hists + _sample(circs, allocate_shots(groups, rest, weights=sg + 1e-12), s2, make_sampler)
    E, sigma, _ = combine_ratio(const, hists, U, W)
    checked = ~keep.all(axis=1)
    kept = (hists * keep)[checked].sum() / max(hists[checked].sum(), 1)
    return {"energy": float(E), "sigma": float(sigma), "shots": int(hists.sum()),
            "n_groups": len(groups), "keep_fraction": float(kept),
            "cz": int(circs[0].count_ops().get("cz", 0))}


def zne_mitigated(isa: QuantumCircuit, H: SparsePauliOp, total_shots: int, lams=(1, 3),
                  seed: int = SEED, **kw) -> dict:
    """ZNE con plegado global de R·C (ver measurement_circuits) y, en cada nivel, post-selección
    y/o lectura (kw va a measure_energy_mitigated). Shots por nivel ∝ |γ_λ| (noise.zne_shot_split)."""
    per_level = zne_shot_split(lams, total_shots)
    meas = [measure_energy_mitigated(isa, H, n, seed=seed + 7 * l, lam=l, **kw)
            for l, n in zip(lams, per_level)]
    E, S = [r["energy"] for r in meas], [r["sigma"] for r in meas]
    e, s = richardson(lams, E, S)
    return {"energy": e, "sigma": s, "lams": list(lams), "E_levels": E, "sigma_levels": S,
            "keep_levels": [r["keep_fraction"] for r in meas],
            "shots": int(sum(r["shots"] for r in meas)), "cz_levels": [r["cz"] for r in meas]}


# --------------------------------------------------------------------------- exacto (sin shots)
def exact_probabilities(isa: QuantumCircuit, groups, p2: float = P_CZ, p1: float = P_1Q,
                        p_ro: float = 0.0, lam: int = 1, fold_rotations: bool = True) -> np.ndarray:
    """Probabilidades EXACTAS de cada resultado de cada circuito de medición con ruido (matriz
    densidad). Solo para elegir la técnica: las energías reportadas salen de shots."""
    sim = AerSimulator(method="density_matrix", noise_model=noise_model(p2, p1))
    circs = []
    for c in measurement_circuits(isa, groups, lam, fold_rotations):
        c = c.remove_final_measurements(inplace=False)
        c.save_probabilities()
        circs.append(c)
    res = sim.run(circs).result()
    P = np.array([np.asarray(res.data(i)["probabilities"]) for i in range(len(circs))])
    if p_ro > 0:
        P = np.array([apply_kron(p, readout_matrices(isa.num_qubits, p_ro)) for p in P])
    return P


def exact_estimate(P: np.ndarray, H: SparsePauliOp, U, W) -> dict:
    """Energía a la que converge el estimador con infinitos shots, y su costo: K = σ²·N con
    Neyman (σ_g por shot de combine_ratio). Shots para σ objetivo = K / σ²."""
    const, _ = split_constant(H)
    E, _, sg = combine_ratio(const, P * 1e12, U, W)   # "1e12 shots": n−1 ≈ n
    return {"E": float(E), "K": float(np.sum(sg) ** 2)}


def echo_fidelity_exact(isa: QuantumCircuit, p2: float = P_CZ, p1: float = P_1Q) -> float:
    """F por pasada a partir de P(0) EXACTA del eco C C†|0> (versión sin shots de noise.echo_fidelity)."""
    from .noise import _inverse
    n = isa.num_qubits
    eco = isa.copy()
    eco.barrier()
    eco.compose(_inverse(isa), inplace=True)
    eco.save_probabilities()
    p0 = AerSimulator(method="density_matrix", noise_model=noise_model(p2, p1)).run(eco).result() \
        .data()["probabilities"][0]
    d = 2 ** n
    return float(np.sqrt(max((p0 - 1 / d) / (1 - 1 / d), 1e-12)))


# --------------------------------------------------------------------------- eco + ZNE con shots
def echo_fidelity_ro(isa: QuantumCircuit, shots: int, seed: int = SEED, readout=None,
                     p2: float = P_CZ, p1: float = P_1Q, p_ro: float = 0.0, make_sampler=None) -> dict:
    """Como noise.echo_fidelity (eco C C†|0>, P(0) = F² + (1 − F²)/2ⁿ) pero midiendo con el
    sampler que tenga error de lectura y corrigiéndolo: P(0) = media de c(x), c = fila 0 de A⁻¹.
    Sin la corrección, el error de lectura se confunde con ruido de compuertas y F sale baja."""
    from .noise import _inverse, _outcome_hist
    make_sampler = make_sampler or (lambda sd: sampler(sd, p2, p1, p_ro))
    n = isa.num_qubits
    eco = isa.copy()
    eco.barrier()
    eco.compose(_inverse(isa), inplace=True)
    eco.measure_all()
    h = _outcome_hist(make_sampler(seed).run([(eco, None, shots)]).result()[0].data.meas, n)
    c = np.zeros(2 ** n)
    c[0] = 1.0
    if readout is not None:  # c(x) = (A⁻¹)[0, x] = ((A⁻¹)ᵀ e_0)(x)
        c = apply_kron(c, [np.linalg.inv(A).T for A in readout])
    p0 = float(h @ c / shots)
    s_p0 = float(np.sqrt(max(h @ c ** 2 / shots - p0 ** 2, 0) / shots))
    d = 2 ** n
    F = float(np.sqrt(max((p0 - 1 / d) / (1 - 1 / d), 1e-12)))
    return {"F": F, "F_sigma": s_p0 / (1 - 1 / d) / (2 * F), "p0": p0, "shots": shots}


def zne_echo(isa: QuantumCircuit, H: SparsePauliOp, total_shots: int, lams=(1, 3),
             echo_shots: int = 2_000_000, seed: int = SEED, n_mc: int = 4000, readout=None, **kw) -> dict:
    """Eco + ZNE (lo que recomienda Octavio) con el plegado de R·C y corrección de lectura:
    cada nivel se reescala E' = (E − (1 − F^λ)·Tr(H)/2ⁿ) / F^λ y después se extrapola.
    Sin post-selección: el reescalado supone que el ruido lleva al estado totalmente mezclado,
    y post-seleccionar cambia esa referencia. σ por Monte Carlo (E_λ y F con su error)."""
    per_level = zne_shot_split(lams, total_shots)
    meas = [measure_energy_mitigated(isa, H, n, seed=seed + 7 * l, lam=l, readout=readout, **kw)
            for l, n in zip(lams, per_level)]
    E = np.array([r["energy"] for r in meas])
    S = np.array([r["sigma"] for r in meas])
    ekw = {k: v for k, v in kw.items() if k in ("p2", "p1", "p_ro", "make_sampler")}
    echo = echo_fidelity_ro(isa, echo_shots, seed=seed + 99, readout=readout, **ekw)
    emix = split_constant(H)[0]
    lam = np.asarray(lams, float)
    resc = lambda e, F: (e - (1 - F ** lam) * emix) / F ** lam
    e, _ = richardson(lams, resc(E, echo["F"]), S)
    rng = np.random.default_rng(seed)
    mc = [richardson(lams, resc(E + rng.normal(0, S), echo["F"] + rng.normal(0, echo["F_sigma"])), S)[0]
          for _ in range(n_mc)]
    return {"energy": float(e), "sigma": float(np.std(mc)), "lams": list(lams), "E_levels": E.tolist(),
            "sigma_levels": S.tolist(), "F": echo["F"], "F_sigma": echo["F_sigma"],
            "shots": int(sum(r["shots"] for r in meas)) + echo_shots}
