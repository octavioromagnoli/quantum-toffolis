"""Extra: ruido de hardware real con FakeFez (IBM Heron r2, 156 qubits, calibración real de IBM).

Dueño: Octavio.

El ansatz compacto (5 qubits, 16 CZ, reduction.py) solo usa CZ entre vecinos de una línea, así que se
ubica en una línea de 5 qubits físicos de Fez sin agregar SWAP. Elegimos la línea con menos errores
esperados según la calibración del backend (best_line). Se simula con el modelo de ruido completo de
Fez (errores de compuerta, T1/T2 y lectura) usando matriz densidad. La medición y la ZNE son las de
noise.py: solo cambian el sampler y la ubicación en el chip.

    line, info = best_line(isa)
    measure_energy_noisy(isa, H5, shots, **fez_kwargs(line))
    zne(isa, H5, shots, lams=(1, 3), **fez_kwargs(line))
"""
from __future__ import annotations

from functools import lru_cache

import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit_aer.noise import NoiseModel
from qiskit_aer.noise.device import basic_device_gate_errors, basic_device_readout_errors
from qiskit_aer.primitives import SamplerV2 as AerSampler
from qiskit_ibm_runtime.fake_provider import FakeFez

from .ansatz import build_ansatz
from .config import ANSATZ_REPS, SEED
from .hamiltonian import molecule, qubit_hamiltonian
from .noise import measure_energy_noisy, readout_calibration, to_line, zne
from .reduction import build_compact_ansatz, optimize_compact, reduced_number_operators, tapering


@lru_cache(maxsize=1)
def fez():
    """(backend, modelo de ruido) de FakeFez, construidos una sola vez."""
    backend = FakeFez()
    return backend, NoiseModel.from_backend(backend)


def _error(target, name: str, qubits) -> float:
    props = target[name].get(tuple(qubits)) or target[name].get(tuple(reversed(qubits)))
    return props.error if props is not None and props.error is not None else 1.0


def _lines(coupling, n: int) -> list[tuple[int, ...]]:
    """Todos los caminos simples de n qubits físicos (en los dos sentidos)."""
    adj: dict[int, set[int]] = {}
    for a, b in coupling.get_edges():
        adj.setdefault(a, set()).add(b)
        adj.setdefault(b, set()).add(a)
    out = []

    def grow(path):
        if len(path) == n:
            out.append(tuple(path))
            return
        for q in adj[path[-1]]:
            if q not in path:
                grow(path + [q])

    for q in adj:
        grow([q])
    return out


def expected_errors(isa: QuantumCircuit, line, backend=None) -> float:
    """Número esperado de errores (a primer orden): Σ error de cada cz/sx/x en su qubit físico
    + error de lectura de cada qubit, con la calibración del backend."""
    backend = backend or fez()[0]
    t = backend.target
    tot = sum(_error(t, "measure", (q,)) for q in line)
    for inst in isa.data:
        name = inst.operation.name
        if name in ("cz", "sx", "x"):
            tot += _error(t, name, [line[isa.find_bit(q).index] for q in inst.qubits])
    return tot


def best_line(isa: QuantumCircuit) -> tuple[tuple[int, ...], dict]:
    """Línea de qubits físicos de Fez con menos errores esperados para este circuito."""
    backend = fez()[0]
    lines = _lines(backend.coupling_map, isa.num_qubits)
    costs = np.array([expected_errors(isa, ln, backend) for ln in lines])
    k = int(np.argmin(costs))
    return lines[k], {"expected_errors": float(costs[k]), "median_expected_errors": float(np.median(costs)),
                      "n_lines": len(lines)}


@lru_cache(maxsize=8)
def line_noise_model(line: tuple[int, ...]) -> NoiseModel:
    """Ruido de Fez solo en los qubits físicos `line`, renumerados 0..n-1 (qubit i = line[i]).
    Para circuitos que solo usan esos qubits es el mismo modelo que el completo (mismas funciones de
    Aer, mismos errores), pero Aer no tiene que cargar los 156 qubits en cada corrida: ~1 s menos
    por circuito medido."""
    backend = fez()[0]
    idx = {q: i for i, q in enumerate(line)}
    nm = NoiseModel(basis_gates=["cz", "id", "rz", "sx", "x"])
    for name, qubits, err in basic_device_gate_errors(target=backend.target):
        if all(q in idx for q in qubits):
            nm.add_quantum_error(err, name, [idx[q] for q in qubits])
    for qubits, err in basic_device_readout_errors(target=backend.target):
        if all(q in idx for q in qubits):
            nm.add_readout_error(err, [idx[q] for q in qubits])
    return nm


def to_fez(line):
    """Ubica circuitos ya nativos de la línea (cz, rz, sx, x) en los qubits físicos `line` de Fez.
    optimization_level=0: no reoptimiza, así se respeta el plegado de la ZNE."""
    backend = fez()[0]
    return lambda circs: transpile(circs, backend, initial_layout=list(line), optimization_level=0)


def fez_kwargs(line, full: bool = False) -> dict:
    """Argumentos para noise.measure_energy_noisy / noise.zne que miden en la línea `line` de Fez.
    full=True: circuito ubicado en el chip de 156 qubits con el modelo completo (lento; para validar)."""
    if full:
        return {"to_device": to_fez(line), "make_sampler": lambda sd: AerSampler(seed=sd, options={
            "backend_options": {"noise_model": fez()[1], "method": "density_matrix"}})}
    nm = line_noise_model(tuple(line))
    return {"make_sampler": lambda sd: AerSampler(seed=sd, options={
        "backend_options": {"noise_model": nm, "method": "density_matrix"}})}


def run_fez(d: float = 2.4, shots: int = 2_000_000, layers: int = 8, seed: int = SEED) -> dict:
    """Un punto estirado en FakeFez con el ansatz compacto: sin mitigar, corrigiendo la lectura y con
    ZNE global (λ = 1, 3, 5, Richardson) + lectura. La ZNE usa 12× shots para una σ parecida."""
    m = molecule(d)
    H8 = qubit_hamiltonian(m)
    H5 = tapering(m, H8).apply(H8)
    ans = build_compact_ansatz(m, layers=layers)
    opt = optimize_compact(ans, H5, reduced_number_operators(m)["S2"], n_random=3, seed=seed)
    isa = to_line(ans.circuit.assign_parameters(opt["x"]))
    line, info = best_line(isa)
    kw = fez_kwargs(line)
    ro = readout_calibration(isa.num_qubits, 500_000, seed=seed + 1, **kw)
    raw = measure_energy_noisy(isa, H5, shots, seed=seed + 2, **kw)
    rc = measure_energy_noisy(isa, H5, shots, seed=seed + 3, readout=ro, **kw)
    z = zne(isa, H5, 12 * shots, lams=(1, 3, 5), readout=ro, seed=seed + 4, **kw)
    uccsd = build_ansatz(m, reps=ANSATZ_REPS).circuit       # θ ≠ 0: con θ = 0 se borran compuertas
    theta = np.random.default_rng(seed).normal(0, 0.1, uccsd.num_parameters)
    uccsd = transpile(uccsd.assign_parameters(theta), fez()[0], optimization_level=1, seed_transpiler=seed)
    f = m["e_fci"]
    rows = [{"método": name, "E": e, "sigma_mHa": 1e3 * s if s is not None else None,
             "E − FCI (mHa)": 1e3 * (e - f), "shots": n}
            for name, e, s, n in [
                ("VQE en Fez, sin mitigar", raw["energy"], raw["sigma"], raw["shots"]),
                ("VQE en Fez, corrigiendo la lectura", rc["energy"], rc["sigma"], rc["shots"]),
                ("VQE en Fez, ZNE (1,3,5) + lectura", z["E"], z["sigma"], z["shots"]),
                ("CCSD", m["e_ccsd"], None, None),
                ("HF", m["e_hf"], None, None)]]
    return {"d": d, "line": line, **info, "cz": int(isa.count_ops().get("cz", 0)),
            "cz_uccsd_fez": int(uccsd.count_ops().get("cz", 0)),
            "err_exact_mHa": 1e3 * (opt["energy"] - f),
            "readout": [(float(A[1, 0]), float(A[0, 1])) for A in ro], "rows": rows}
