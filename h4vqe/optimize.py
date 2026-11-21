"""Búsqueda de parámetros SIN ruido (permitido por la consigna).

Dueño: Gonzalo.

Los θ* se buscan con el simulador numpy exacto del ansatz. La energía que se
REPORTA no sale de acá: sale de medir el circuito con shots (estimator.py).
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import minimize

from .ansatz import Ansatz


def optimize_params(ans: Ansatz, H, x0=None, tol: float = 1e-10, maxiter: int = 2000,
                    n_restarts: int = 0, seed: int = 0) -> dict:
    """L-BFGS-B con gradiente analítico. Devuelve dict con x, energy, nfev.

    Siempre prueba también θ = 0 (|HF⟩), además del x0 dado (warm start de la distancia
    anterior), y se queda con la energía más baja (es variacional, así que el mínimo es el
    mejor). Con la lista ADAPT de 15 operadores, el warm start solo cae en mínimos malos
    a 2,6-3,0 Å (hasta 2 mHa), y desde θ = 0 da 0,05-0,15 mHa.
    """
    Hm = H.to_matrix(sparse=True)
    fun = lambda x: ans.fast.energy_and_grad(x, Hm)
    zero = np.zeros(ans.num_parameters)
    bases = [zero] if x0 is None else [np.asarray(x0, float), zero]
    rng = np.random.default_rng(seed)
    starts = list(bases)
    starts += [b + rng.normal(0, 0.1, ans.num_parameters) for b in bases for _ in range(n_restarts)]
    best = None
    for s in starts:
        r = minimize(fun, s, jac=True, method="L-BFGS-B",
                     options={"maxiter": maxiter, "ftol": tol, "gtol": 1e-8})
        if best is None or r.fun < best.fun:
            best = r
    return {"x": best.x, "energy": float(best.fun), "nfev": int(best.nfev),
            "success": bool(best.success)}


def exact_state(ans: Ansatz, x) -> np.ndarray:
    """Solo para diagnóstico (fidelidad, varianza). NUNCA para la energía reportada."""
    return ans.fast.state(x)
