"""Extra "comparar ansätze" (issue #2): todos los ansätze en las 18 distancias, SIN ruido.

Para cada ansatz: qubits, parámetros, CZ (todos con todos y en línea), grupos de medición,
error máximo contra FCI (en toda la curva y en la zona estirada) y S² máximo.
Los θ se optimizan como en el pipeline: warm start + θ = 0 (optimize_params) o, para el
compacto, warm start + arranques al azar con penalización de S² (optimize_compact).
FCI se usa SOLO para reportar el error.

    python scripts/compare_ansatze.py               # ~25 min en una laptop (UCCSD ×2 tarda ~7)
    python scripts/compare_ansatze.py --rapido      # sin UCCSD ×2 ni el barrido de capas
Dueño: Gonzalo.
"""
from __future__ import annotations

import json
import os
import sys
import time

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.circuit import ParameterVector
from qiskit.circuit.library import PauliEvolutionGate

from h4vqe.adapt import ADAPT_FERMIONIC, ADAPT_FERMIONIC_13, ADAPT_QEB
from h4vqe.ansatz import Ansatz, FastAnsatz, build_ansatz, cz_count, excitation_generator
from h4vqe.config import BASIS_GATES, DISTANCES, SEED, STRETCHED_FROM
from h4vqe.estimator import measurement_groups
from h4vqe.hamiltonian import molecule, number_operators, qubit_hamiltonian
from h4vqe.optimize import optimize_params
from h4vqe.reduction import (build_compact_ansatz, optimize_compact, reduced_number_operators,
                             tapering)
from h4vqe.routing import LINE_ORDER, line_cz_count, routed_ansatz

D_CZ = 2.4   # geometría donde se cuentan CZ y grupos (el circuito es el mismo en toda la curva)


def tapered_adapt(m, H, excitations, n_qubits):
    """ADAPT con los generadores reducidos a 6 o 5 qubits (tapering de reduction.py)."""
    tp = tapering(m, H, n_qubits=n_qubits)
    gens = [tp.apply(excitation_generator(e, 2 * m["norb"])) for e in excitations]
    th = ParameterVector("θ", len(gens))
    qc = QuantumCircuit(tp.n_qubits)
    for q in range(tp.n_qubits):
        if (tp.hf_index >> q) & 1:
            qc.x(q)
    for t, G in zip(th, gens):
        qc.append(PauliEvolutionGate(G, time=t), range(tp.n_qubits))
    native = transpile(qc, basis_gates=BASIS_GATES, optimization_level=1, seed_transpiler=SEED)
    ans = Ansatz(native, list(excitations), FastAnsatz(gens, tp.hf_index, tp.n_qubits), th,
                 {"n_qubits": tp.n_qubits})
    return ans, tp


def summarize(name, errs, s2, row):
    errs = np.array(errs)
    st = np.array(DISTANCES) >= STRETCHED_FROM - 1e-9
    row |= {"ansatz": name, "err_max_mHa": float(errs.max()),
            "err_max_estirada_mHa": float(errs[st].max()), "S2_max": float(max(s2))}
    print(f"{name:32s} q={row['qubits']} θ={row['params']:2d} CZ={row['cz_all']:4d}/{row['cz_line']:4d} "
          f"grupos={row['groups']:2d} err máx={row['err_max_mHa']:8.3f} mHa "
          f"(estirada {row['err_max_estirada_mHa']:.3f})  S²≤{row['S2_max']:.1e}", flush=True)
    return row


def run_full(name, kw=None, n_qubits=8, line=None, ms=None, Hs=None):
    """Ansatz de excitaciones en 8 qubits (kw para build_ansatz) o reducido (n_qubits = 6, 5)."""
    x0, errs, s2 = None, [], []
    for m, H in zip(ms, Hs):
        if n_qubits == 8:
            a, Hr = build_ansatz(m, **kw), H
            S2 = number_operators(m["norb"])["S2"]
        else:
            a, tp = tapered_adapt(m, H, ADAPT_FERMIONIC, n_qubits)
            Hr, S2 = tp.apply(H), reduced_number_operators(m, n_qubits)["S2"]
        o = optimize_params(a, Hr, x0=x0, n_restarts=2)
        x0 = o["x"]
        psi = a.fast.state(o["x"])
        errs.append(1e3 * (o["energy"] - m["e_fci"]))
        s2.append(float(np.real(np.vdot(psi, S2.to_matrix(sparse=True) @ psi))))
    k = DISTANCES.index(D_CZ)
    m, H = ms[k], Hs[k]
    if n_qubits == 8:
        a, Hr = build_ansatz(m, **kw), H
    else:
        a, tp = tapered_adapt(m, H, ADAPT_FERMIONIC, n_qubits)
        Hr = tp.apply(H)
    row = {"qubits": n_qubits, "params": a.num_parameters, "cz_all": cz_count(a.circuit),
           "cz_line": line if line is not None else min(cz_count(a.circuit, "line", seed=s) for s in range(3)),
           "groups": len(measurement_groups(Hr))}
    return summarize(name, errs, s2, row)


def run_compact(layers, n_qubits, ms, Hs):
    x0, errs, s2 = None, [], []
    for m, H in zip(ms, Hs):
        tp = tapering(m, H, n_qubits=n_qubits)
        Hr, S2 = tp.apply(H), reduced_number_operators(m, n_qubits)["S2"]
        a = build_compact_ansatz(m, layers=layers, n_qubits=n_qubits)
        o = optimize_compact(a, Hr, S2, x0=x0, n_random=4)
        x0 = o["x"]
        psi = a.fast.state(o["x"])
        errs.append(1e3 * (o["energy"] - m["e_fci"]))
        s2.append(float(np.real(np.vdot(psi, S2.to_matrix(sparse=True) @ psi))))
    k = DISTANCES.index(D_CZ)
    a = build_compact_ansatz(ms[k], layers=layers, n_qubits=n_qubits)
    row = {"qubits": n_qubits, "params": a.num_parameters, "cz_all": cz_count(a.circuit),
           "cz_line": cz_count(a.circuit, "line"),
           "groups": len(measurement_groups(tapering(ms[k], Hs[k], n_qubits).apply(Hs[k])))}
    return summarize(f"compacto {n_qubits}q L={layers}", errs, s2, row)


def main(rapido: bool = False):
    t0 = time.time()
    ms = [molecule(d) for d in DISTANCES]
    Hs = [qubit_hamiltonian(m) for m in ms]
    k = DISTANCES.index(D_CZ)
    l15 = line_cz_count(routed_ansatz(ms[k], ADAPT_FERMIONIC, LINE_ORDER)[0])
    l13 = line_cz_count(routed_ansatz(ms[k], ADAPT_FERMIONIC_13, LINE_ORDER)[0])
    kw = {"ms": ms, "Hs": Hs}
    rows = [run_full("UCCSD ×1", {"reps": 1}, **kw)]
    if not rapido:
        rows.append(run_full("UCCSD ×2", {"reps": 2}, **kw))
    rows += [run_full("ADAPT 13", {"excitations": ADAPT_FERMIONIC_13}, line=l13, **kw),
             run_full("ADAPT 15 (pares de espín) + fSWAP", {"excitations": ADAPT_FERMIONIC},
                      line=l15, **kw),
             run_full("ADAPT 15 reducido a 6 qubits", n_qubits=6, **kw),
             run_full("ADAPT 15 reducido a 5 qubits", n_qubits=5, **kw),
             run_full("ADAPT-QEB 15", {"excitations": ADAPT_QEB, "kind": "qeb"}, **kw)]
    for L in ([7] if rapido else [5, 6, 7, 8]):
        rows.append(run_compact(L, 5, ms, Hs))
    out = os.path.join(ROOT, "scripts", "compare_ansatze.json")
    with open(out, "w", encoding="utf-8") as f:
        json.dump(rows, f, ensure_ascii=False, indent=1)
    print(f"escrito {out}  ({time.time() - t0:.0f} s)")


if __name__ == "__main__":
    main(rapido="--rapido" in sys.argv)
