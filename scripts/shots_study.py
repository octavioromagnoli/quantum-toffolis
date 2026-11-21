"""Cuántos shots hacen falta para σ = 0,5 mHa según agrupamiento y asignación (issue #3).

Cuenta ANALÍTICA, sin simular shots: para un estado |ψ> y grupos medidos con
fracciones f_g del total N, σ² = Σ_g V_g / (f_g N), V_g = <O_g²> − <O_g>².
Definimos K = σ² N (Ha²·shot); shots para σ objetivo = K / σ².
Neyman con V_g exactas es el óptimo: K = (Σ_g √V_g)².

|ψ> = estado FCI (sector N_α = N_β = 2): el VQE llega a < 1 µHa de FCI, así que
las varianzas son las mismas. Solo para planificar; la energía reportada sale de shots.

    python scripts/shots_study.py
"""
from __future__ import annotations

import os
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import numpy as np
from qiskit.quantum_info import SparsePauliOp

from h4vqe.estimator import measurement_groups, split_constant
from h4vqe.hamiltonian import molecule, qubit_hamiltonian

SIGMA_TARGET = 0.5e-3


def sector_ground_state(H: SparsePauliOp, norb: int, na=2, nb=2):
    nq = H.num_qubits
    idx = [x for x in range(2 ** nq)
           if bin(x & ((1 << norb) - 1)).count("1") == na and bin(x >> norb).count("1") == nb]
    Hm = H.to_matrix(sparse=True)[idx][:, idx].toarray()
    w, v = np.linalg.eigh(Hm)
    psi = np.zeros(2 ** nq, complex)
    psi[idx] = v[:, 0]
    return w[0], psi


def hf_state(m):
    norb = m["norb"]
    psi = np.zeros(2 ** (2 * norb), complex)
    psi[sum(1 << q for q in list(range(2)) + [norb + i for i in range(2)])] = 1
    return psi


class Terms:
    """Paulis del hamiltoniano (sin la identidad) con P_k|ψ> precalculado."""

    def __init__(self, H: SparsePauliOp):
        _, rest = split_constant(H)
        self.labels = [p.to_label() for p in rest.paulis]
        self.paulis = list(rest.paulis)
        self.c = np.real(rest.coeffs)
        self.n = len(self.c)
        self._mats = [SparsePauliOp(p).to_matrix(sparse=True) for p in self.paulis]

    def apply(self, psi):
        return np.array([M @ psi for M in self._mats])  # (n_terms, dim)


def group_var(Pv, psi, idx, coef):
    """Varianza de O = Σ coef_k P_k (k en idx) en |ψ>, con Pv[k] = P_k|ψ>."""
    Opsi = coef @ Pv[idx]
    mean = np.real(np.vdot(psi, Opsi))
    return max(np.real(np.vdot(Opsi, Opsi)) - mean ** 2, 0.0)


def qwc(a: str, b: str) -> bool:
    return all(x == "I" or y == "I" or x == y for x, y in zip(a, b))


def sorted_insertion(T: Terms, qubit_wise: bool):
    """Sorted insertion (Crawford et al. 2021): por |c| decreciente, cada término al
    primer grupo con el que conmuta. Suele dar menos varianza que el coloreo."""
    order = np.argsort(-np.abs(T.c))
    groups = []
    for k in order:
        for g in groups:
            if all(qwc(T.labels[k], T.labels[j]) if qubit_wise else T.paulis[k].commutes(T.paulis[j])
                   for j in g):
                g.append(k)
                break
        else:
            groups.append([k])
    return groups


def estimator_groups(T: Terms, H, method: str, grouping: str):
    """Los MISMOS grupos que usa h4vqe/estimator.py, como listas de índices de términos."""
    lab2k = {l: k for k, l in enumerate(T.labels)}
    return [[lab2k[p.to_label()] for p in g.paulis]
            for g in measurement_groups(H, method=method, grouping=grouping)]


def K_alloc(V, f):
    return float(np.sum(V / f))


def study_groups(T, groups, Pv, psi, Pv_hf, psi_hf):
    V = np.array([group_var(Pv, psi, g, T.c[g]) for g in groups])
    Vhf = np.array([group_var(Pv_hf, psi_hf, g, T.c[g]) for g in groups])
    sumc = np.array([np.abs(T.c[g]).sum() for g in groups])
    f_sumc = sumc / sumc.sum()
    f_ney = np.sqrt(V) / np.sqrt(V).sum()
    f_hf = np.sqrt(Vhf) + 1e-3 * sumc  # HF tiene V_g = 0 en grupos diagonales: piso chico
    f_hf = f_hf / f_hf.sum()
    f_2st = 0.1 * f_sumc + 0.9 * f_ney  # lo que hace measure_energy (piloto 10 % + Neyman)
    return {"n_groups": len(groups), "sumc": K_alloc(V, f_sumc), "neyman_hf": K_alloc(V, f_hf),
            "2etapas": K_alloc(V, f_2st), "neyman": K_alloc(V, f_ney)}


def eff_basis(T, g):
    """Base efectiva de un grupo QWC: letra por qubit (los qubits libres se miden en Z)."""
    b = ["Z"] * len(T.labels[0])
    for k in g:
        for q, ch in enumerate(T.labels[k]):
            if ch != "I":
                b[q] = ch
    return "".join(b)


def shot_reuse(T, groups, Pv, psi, iters=30):
    """Cada base QWC medida también estima todos los Paulis compatibles de OTROS grupos.
    Estimador: P_k = promedio ponderado (por shots) sobre las bases compatibles.
    E = Σ_g media_g(O'_g), O'_g = Σ_{k compat. g} w_kg c_k P_k, Σ_g w_kg = 1.
    Asignación: punto fijo de Neyman sobre las V'_g."""
    bases = [eff_basis(T, g) for g in groups]
    C = np.array([[all(ch == "I" or ch == b[q] for q, ch in enumerate(T.labels[k])) for b in bases]
                  for k in range(T.n)], float)  # (terms, groups)
    V0 = np.array([group_var(Pv, psi, g, T.c[g]) for g in groups])
    f = np.sqrt(V0) / np.sqrt(V0).sum()
    best = np.inf
    for _ in range(iters):
        W = C * f[None, :]
        W = W / W.sum(1, keepdims=True)
        V = np.array([group_var(Pv, psi, np.arange(T.n), W[:, g] * T.c) for g in range(len(groups))])
        best = min(best, K_alloc(V, np.maximum(f, 1e-12)))
        f = np.sqrt(V) + 1e-12
        f = f / f.sum()
    return best


def main():
    print(f"K = σ²·N en mHa²·Mshots ; shots para σ = {1e3*SIGMA_TARGET} mHa = K / 0,25 (en M)\n")
    for d in [0.9, 1.8, 2.0, 2.4, 3.0]:
        m = molecule(d)
        H = qubit_hamiltonian(m)
        e0, psi = sector_ground_state(H, m["norb"])
        psi_hf = hf_state(m)
        T = Terms(H)
        Pv, Pv_hf = T.apply(psi), T.apply(psi_hf)
        res = {}
        for name, groups in [("QWC qiskit", estimator_groups(T, H, "coloring", "qwc")),
                             ("QWC sorted-ins.", estimator_groups(T, H, "sorted", "qwc")),
                             ("GC qiskit", estimator_groups(T, H, "coloring", "gc")),
                             ("GC sorted-ins.", estimator_groups(T, H, "sorted", "gc"))]:
            res[name] = study_groups(T, groups, Pv, psi, Pv_hf, psi_hf)
            if name.startswith("QWC"):
                res[name]["+reuso"] = shot_reuse(T, groups, Pv, psi)
        print(f"d = {d} Å   (E_FCI − e0 = {1e3*(m['e_fci']-e0):.1e} mHa)")
        conv = 1e6 / 1e6  # Ha²·shot -> mHa²·Mshot
        ref = res["QWC qiskit"]["2etapas"]
        print(f"  {'agrupamiento':16s} {'grupos':>6s}  " + "  ".join(f"{k:>10s}" for k in
              ["sumc", "neyman_hf", "2etapas", "neyman", "+reuso"]) + "   shots σ=0,5 (mejor)   ganancia vs hoy")
        for name, r in res.items():
            vals = [r.get(k) for k in ["sumc", "neyman_hf", "2etapas", "neyman", "+reuso"]]
            best = min(v for v in vals if v is not None)
            print(f"  {name:16s} {r['n_groups']:6d}  " +
                  "  ".join(f"{v*conv:10.3f}" if v is not None else f"{'-':>10s}" for v in vals) +
                  f"   {best/SIGMA_TARGET**2/1e6:8.2f} M          ×{ref/best:.2f}")
        print()


if __name__ == "__main__":
    main()
