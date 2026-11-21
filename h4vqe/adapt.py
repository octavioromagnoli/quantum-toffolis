"""ADAPT-VQE: arma el ansatz agregando de a una excitación (la de gradiente máximo).

Dueño: Gonzalo.

Idea (Grimsley et al., Nat. Commun. 10, 3007, 2019):
    1. Arrancar de |HF> con el ansatz vacío.
    2. Para cada operador G_k del pool, el gradiente de agregar exp(-i θ G_k) al final
       (con θ = 0) es  g_k = <ψ| i[G_k, H] |ψ> = 2 Re <Hψ| -i G_k ψ>.
    3. Agregar el de |g_k| máximo y reoptimizar TODOS los parámetros (warm start).
    4. Repetir hasta que max|g_k| < grad_tol (criterio que NO usa FCI).
       Ojo: con 1e-4 se traba en mesetas (a 1,6, 2,8 y 3,0 Å para con 0,6 a 2 mHa de
       error y después el gradiente vuelve a crecer). Con 1e-5 llega a FCI en toda la curva.

Un mismo operador se puede elegir más de una vez.

La búsqueda es sin ruido con `FastAnsatz` (permitido por la consigna). La lista que
devuelve se pasa a `build_ansatz(m, excitations=lista)`, que arma el circuito
que después se mide con shots.

Funciones:
    adapt_vqe(m, H, ...)         -> lista propia para UNA geometría
    adapt_vqe_multi(ms, Hs, ...) -> UNA lista común para varias geometrías
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import minimize

from .ansatz import (FastAnsatz, excitation_generator, hf_occupation, qubit_annihilation_ops,
                     uccsd_excitations)
from .hamiltonian import annihilation_ops


def fermionic_pool(m: dict) -> list:
    """Las 26 excitaciones simples y dobles ocupado -> vacío de UCCSD (conservan S_z)."""
    return uccsd_excitations(m)


class _Pool:
    """Matrices G_k y G_k² del pool y estado HF, calculadas una sola vez."""

    def __init__(self, m: dict, excitations: list, kind: str = "fermionic"):
        nq = 2 * m["norb"]
        a = annihilation_ops(nq) if kind == "fermionic" else qubit_annihilation_ops(nq)
        self.excitations = list(excitations)
        self.G = [excitation_generator(e, nq, a, kind).to_matrix(sparse=False)
                  for e in excitations]
        self.G2 = [G @ G for G in self.G]
        self.psi0 = np.zeros(2 ** nq, dtype=complex)
        self.psi0[sum(1 << q for q in hf_occupation(m))] = 1.0

    def fast(self, idx: list[int]) -> FastAnsatz:
        return FastAnsatz.from_matrices([self.G[i] for i in idx], [self.G2[i] for i in idx],
                                        self.psi0)

    def gradients(self, psi: np.ndarray, Hpsi: np.ndarray) -> np.ndarray:
        return np.array([2.0 * np.real(np.vdot(Hpsi, -1j * (G @ psi))) for G in self.G])


def _reoptimize(fast: FastAnsatz, Hm, x0: np.ndarray) -> tuple[np.ndarray, float]:
    r = minimize(lambda x: fast.energy_and_grad(x, Hm), x0, jac=True, method="L-BFGS-B",
                 options={"maxiter": 5000, "ftol": 1e-14, "gtol": 1e-9})
    return r.x, float(r.fun)


def adapt_vqe(m: dict, H, pool: list | None = None, grad_tol: float = 1e-5,
              max_ops: int = 40, e_ref: float | None = None, err_tol: float | None = None,
              kind: str = "fermionic", verbose: bool = False) -> dict:
    """ADAPT-VQE para una geometría.

    Criterio de parada principal: max|g_k| < grad_tol (no usa FCI).
    `e_ref`/`err_tol` (por ejemplo FCI y 0,5 mHa) son opcionales y solo para diagnóstico.
    """
    pool = _Pool(m, pool or fermionic_pool(m), kind)
    Hm = H.to_matrix(sparse=True)
    idx: list[int] = []
    x = np.zeros(0)
    psi = pool.psi0
    energy = float(np.real(np.vdot(psi, Hm @ psi)))
    history = [{"n_ops": 0, "energy": energy, "max_grad": None, "added": None}]
    reason = "max_ops"
    for _ in range(max_ops):
        g = pool.gradients(psi, Hm @ psi)
        k = int(np.argmax(np.abs(g)))
        history[-1]["max_grad"] = float(abs(g[k]))
        if abs(g[k]) < grad_tol:
            reason = "grad_tol"
            break
        if e_ref is not None and err_tol is not None and energy - e_ref < err_tol:
            reason = "err_tol"
            break
        idx.append(k)
        fast = pool.fast(idx)
        x, energy = _reoptimize(fast, Hm, np.append(x, 0.0))
        psi = fast.state(x)
        history.append({"n_ops": len(idx), "energy": energy, "max_grad": None,
                        "added": pool.excitations[k]})
        if verbose:
            err = f"  err={1e3 * (energy - e_ref):.4f} mHa" if e_ref is not None else ""
            print(f"  op {len(idx):2d}: {pool.excitations[k]}  |g|={abs(g[k]):.2e}  "
                  f"E={energy:.8f}{err}")
    return {"excitations": [pool.excitations[i] for i in idx], "x": x, "energy": energy,
            "history": history, "stop": reason}


def line_cost(order: list[int], exc) -> tuple[float, list[int]]:
    """Costo en CZ de agregar `exc` con qubits en línea (ruteo con fSWAP de routing.py).

    Devuelve (CZ aproximadas, nuevo orden de modos). Doble en línea ≈ 18 CZ, simple ≈ 2,
    cada fSWAP ≈ 2.
    """
    from qiskit import QuantumCircuit
    from .routing import _gather
    new = list(order)
    n_sw = _gather(QuantumCircuit(len(order)), new, list(exc[0]) + list(exc[1]))
    return (18.0 if len(exc[0]) == 2 else 2.0) + 2.0 * n_sw, new


def adapt_vqe_multi(ms: list[dict], Hs: list, pool: list | None = None,
                    grad_tol: float = 1e-5, max_ops: int = 40, kind: str = "fermionic",
                    init_order: list[int] | None = None, cost_power: float = 0.0,
                    verbose: bool = False) -> dict:
    """UNA lista de excitaciones común a varias geometrías.

    En cada paso se calcula el gradiente de cada operador en TODAS las geometrías y se
    agrega el de mayor max_d |g_k(d)|; después se reoptimiza cada geometría por separado
    (sus θ son distintos, pero el circuito y la cantidad de CZ son los mismos).
    Para cuando el gradiente máximo es < grad_tol en todas las geometrías.

    Co-ADAPT (Ramôa et al. 2026): con cost_power > 0 el puntaje es
    max_d|g_k| / costo_k^cost_power, donde costo_k son las CZ en línea de agregar el
    operador dado el orden actual de modos (seguido como en routing.routed_ansatz,
    arrancando en init_order). Así la lista sale barata de rutear EN ESE ORDEN.
    """
    pools = [_Pool(m, pool or fermionic_pool(m), kind) for m in ms]
    Hms = [H.to_matrix(sparse=True) for H in Hs]
    idx: list[int] = []
    xs = [np.zeros(0) for _ in ms]
    psis = [p.psi0 for p in pools]
    energies = [float(np.real(np.vdot(s, Hm @ s))) for s, Hm in zip(psis, Hms)]
    history = []
    reason = "max_ops"
    order = list(init_order or range(2 * ms[0]["norb"]))
    for _ in range(max_ops):
        G = np.array([p.gradients(s, Hm @ s) for p, s, Hm in zip(pools, psis, Hms)])
        gmax = np.max(np.abs(G), axis=0)
        if gmax.max() < grad_tol:
            reason = "grad_tol"
            break
        costs = [line_cost(order, e) for e in pools[0].excitations]
        score = gmax / np.array([c for c, _ in costs]) ** cost_power
        score[gmax < grad_tol] = -1.0
        k = int(np.argmax(score))
        order = costs[k][1]
        idx.append(k)
        for j, (p, Hm) in enumerate(zip(pools, Hms)):
            fast = p.fast(idx)
            xs[j], energies[j] = _reoptimize(fast, Hm, np.append(xs[j], 0.0))
            psis[j] = fast.state(xs[j])
        history.append({"n_ops": len(idx), "added": pools[0].excitations[k],
                        "max_grad": float(gmax[k]), "line_cost": costs[k][0],
                        "energies": list(energies)})
        if verbose:
            print(f"  op {len(idx):2d}: {pools[0].excitations[k]}  max|g|={gmax[k]:.2e}  "
                  f"costo línea≈{costs[k][0]:.0f}")
    return {"excitations": [pools[0].excitations[i] for i in idx], "xs": xs,
            "energies": energies, "history": history, "stop": reason}


# --------------------------------------------------------------------------------------
# Listas ya calculadas (para no correr ADAPT en Colab). Se regeneran con:
#   adapt_vqe_multi([molecule(d) for d in DISTANCES], [...], kind=..., grad_tol=1e-5)
# y truncando a los primeros k operadores. Excitaciones: ((ocupados), (vacíos)),
# qubits 0-3 = alfa, 4-7 = beta (convención de hamiltonian.py).
# --------------------------------------------------------------------------------------

# Recomendado. Pool fermiónico, lista común para las 18 distancias, primeros 13 operadores.
# Sin ruido y con el warm start del pipeline: error máx. 0,84 mHa (1,4 Å), ≤ 0,6 mHa
# para d ≥ 1,8 Å; 185 CZ todos con todos, ~515 en línea (UCCSD×2: ~584 / ~1600).
# Corte elegido por la varianza ⟨H²⟩ − ⟨H⟩² (cae de 3e-3 a 2,6e-4 Ha² en k = 13),
# que se puede medir sin conocer FCI.
ADAPT_FERMIONIC = [
    ((1, 5), (3, 7)), ((0, 5), (3, 6)), ((0, 4), (2, 6)), ((1, 4), (2, 7)),
    ((1, 5), (2, 6)), ((0, 5), (2, 7)), ((0, 4), (3, 7)), ((1, 4), (3, 6)),
    ((0, 1), (2, 3)), ((4, 5), (6, 7)), ((0,), (2,)), ((5,), (7,)), ((0, 5), (2, 7)),
]

# Extra. Pool QEB (excitaciones de qubit), primeros 15 operadores: 159 CZ / ~451 en línea.
# OJO: con el warm start del pipeline cae en un mínimo malo a 2,0 Å (9,8 mHa, S² = 0,78).
# Usar con build_ansatz(m, excitations=ADAPT_QEB, kind="qeb") y varios reinicios.
ADAPT_QEB = [
    ((1, 5), (3, 7)), ((0, 5), (3, 6)), ((0, 4), (2, 6)), ((1, 4), (2, 7)),
    ((1, 5), (2, 6)), ((0, 5), (2, 7)), ((0, 4), (3, 7)), ((1, 4), (3, 6)),
    ((0, 1), (2, 3)), ((4, 5), (6, 7)), ((0,), (2,)), ((5,), (7,)), ((4,), (6,)),
    ((1,), (3,)), ((0, 5), (2, 7)),
]
