"""Reducción por simetrías (8 -> 5 qubits) y ansatz compacto nativo de la línea.

Dueño: Octavio.

Simetrías Z2 del H4 lineal en Jordan-Wigner (todas de tipo Z; conmutan con H, N, S_z y S^2):
    P_alfa = Z0 Z1 Z2 Z3                  paridad de N_alfa
    P_beta = Z4 Z5 Z6 Z7                  paridad de N_beta
    P_u    = Π Z sobre los orbitales impares ante la inversión de la cadena (σu), α y β
Tapering (Bravyi, Gambetta, Mezzacapo, Temme 2017): se lleva el conjunto de generadores a forma
escalonada sobre GF(2), así cada generador g_i tiene un qubit "exclusivo" q_i. Con
U_i = (X_{q_i} + g_i)/√2 se cumple U_i g_i U_i = X_{q_i}; en U H U el qubit q_i solo aparece como
I o X, y X_{q_i} se reemplaza por el autovalor de g_i en el sector del estado HF. Quedan 5 qubits.
El estado HF reducido sigue siendo un bitstring (los bits HF de los qubits que quedan).

Ansatz compacto ("brickwork"): X en HF, una capa de RY y L capas de [CZ en ladrillo + RY].
Cada capa tiene 2 CZ entre vecinos (qubits en línea: 0 SWAP). Con L = 8 son 16 CZ y 45 parámetros,
y llega a FCI (< 0,001 mHa) en toda la curva. OJO: la versión con TODAS las CZ en cada capa
(cadena 0-1,1-2,2-3,3-4 simultánea) tiene un techo de fidelidad 0,961 que no mejora con más capas.

Interfaz (la misma que ansatz.Ansatz, así optimize/pipeline/estimator funcionan igual):
    tapering(m) -> Tapering           .apply(op8) -> op5, .hf_index, .removed, .kept
    reduced_hamiltonian(m) -> SparsePauliOp (5 qubits)
    reduced_number_operators(m) -> {"N_alpha","N_beta","S2"} en 5 qubits
    build_compact_ansatz(m, layers=8) -> Ansatz (5 qubits)
    optimize_compact(ans, H5, S2_5=None, mu=0.5, x0=None, n_random=3) -> dict como optimize_params
"""
from __future__ import annotations

import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.circuit import ParameterVector
from qiskit.quantum_info import SparsePauliOp
from scipy.optimize import minimize

from .ansatz import Ansatz, build_ansatz, cz_count, hf_occupation
from .config import ANSATZ_REPS, BASIS_GATES, SEED
from .estimator import measurement_groups
from .hamiltonian import molecule, number_operators, qubit_hamiltonian
from .optimize import optimize_params


# --------------------------------------------------------------------------- simetrías
def inversion_parities(m: dict) -> np.ndarray:
    """+1 / -1 para cada orbital molecular ante la inversión de la cadena (átomo i <-> átomo 3-i).
    En STO-3G cada H aporta un solo orbital 1s, así que invertir = dar vuelta los coeficientes."""
    C = m["hf"].mo_coeff
    S = m["mol"].intor("int1e_ovlp")
    par = np.einsum("ip,ij,jp->p", C[::-1, :], S, C)        # <φ_p | Inv | φ_p>
    out = np.sign(np.round(par, 6))
    if not np.allclose(np.abs(par), 1, atol=1e-5):
        raise ValueError(f"orbitales sin paridad definida ante inversión: {par}")
    return out


def _zstring(qubits, n):
    return SparsePauliOp.from_sparse_list([("Z" * len(qubits), sorted(qubits), 1.0)], n)


def z2_symmetries(m: dict) -> list[SparsePauliOp]:
    norb = m["norb"]
    nq = 2 * norb
    u = [p for p, s in enumerate(inversion_parities(m)) if s < 0]
    return [_zstring(range(norb), nq), _zstring(range(norb, nq), nq),
            _zstring(u + [norb + p for p in u], nq)]


def _commutes(A: SparsePauliOp, B: SparsePauliOp) -> bool:
    c = (A.dot(B) - B.dot(A)).simplify(atol=1e-10)
    return bool(np.allclose(c.coeffs, 0, atol=1e-10))


def _rref_gf2(rows: np.ndarray):
    """Forma escalonada reducida sobre GF(2). Devuelve (filas, columnas pivote)."""
    R = rows.copy() % 2
    pivots, r = [], 0
    for c in range(R.shape[1]):
        hit = [i for i in range(r, R.shape[0]) if R[i, c]]
        if not hit:
            continue
        R[[r, hit[0]]] = R[[hit[0], r]]
        for i in range(R.shape[0]):
            if i != r and R[i, c]:
                R[i] ^= R[r]
        pivots.append(c)
        r += 1
        if r == R.shape[0]:
            break
    return R[:r], pivots


PREFERRED_REMOVED = (0, 4, 7)   # la elección validada en el análisis (ver docs/ideas/diferencial.md)


class Tapering:
    def __init__(self, generators: list[SparsePauliOp], hf_index: int,
                 prefer_removed=PREFERRED_REMOVED):
        nq = generators[0].num_qubits
        Z = np.array([g.paulis[0].z.astype(np.uint8) for g in generators])
        # Las columnas preferidas van primero: si son independientes, quedan como pivotes (qubits que
        # se eliminan). La elección no cambia la física (todas son unitariamente equivalentes), pero sí
        # la forma de los términos: cuántos grupos de medición hay y qué tan bien entrena el ansatz.
        prefer = [q for q in prefer_removed if q < nq]
        order = prefer + [q for q in range(nq) if q not in prefer]
        R, piv = _rref_gf2(Z[:, order])
        self.removed = [order[p] for p in piv]
        R_cols = np.zeros_like(R)
        R_cols[:, order] = R
        self.generators = [_zstring([q for q in range(nq) if row[q]], nq) for row in R_cols]
        self.sector = [(-1) ** sum((hf_index >> q) & 1 for q in range(nq) if row[q]) for row in R_cols]
        self.kept = [q for q in range(nq) if q not in self.removed]
        self.n_qubits = len(self.kept)
        self.U = []
        for g, q in zip(self.generators, self.removed):
            X = SparsePauliOp.from_sparse_list([("X", [q], 1.0)], nq)
            self.U.append(((X + g) / np.sqrt(2)).simplify())
        self.hf_index = sum(((hf_index >> q) & 1) << k for k, q in enumerate(self.kept))

    def apply(self, op: SparsePauliOp) -> SparsePauliOp:
        """Operador de 8 qubits que conmuta con las simetrías -> operador de 5 qubits en el sector HF."""
        for U in self.U:
            op = U.dot(op).dot(U).simplify(atol=1e-12)
        labels, coefs = [], []
        for p, c in zip(op.paulis, op.coeffs):
            for q, s in zip(self.removed, self.sector):
                if p.z[q]:
                    raise ValueError("el operador no conmuta con las simetrías")
                if p.x[q]:
                    c = c * s
            lab = p.to_label()[::-1]                      # lab[k] = letra del qubit k
            labels.append("".join(lab[k] for k in self.kept)[::-1])
            coefs.append(c)
        return SparsePauliOp(labels, coefs).simplify(atol=1e-12).chop(1e-12)


def tapering(m: dict, H: SparsePauliOp | None = None) -> Tapering:
    gens = z2_symmetries(m)
    if H is not None:
        for g in gens:
            assert _commutes(H, g), "una simetría no conmuta con H"
    hf_index = sum(1 << q for q in hf_occupation(m))
    return Tapering(gens, hf_index)


def reduced_hamiltonian(m: dict, H: SparsePauliOp | None = None) -> SparsePauliOp:
    H = qubit_hamiltonian(m) if H is None else H
    return tapering(m, H).apply(H)


def reduced_number_operators(m: dict) -> dict[str, SparsePauliOp]:
    tp = tapering(m)
    return {k: tp.apply(v) for k, v in number_operators(m["norb"]).items()}


# --------------------------------------------------------------------------- ansatz compacto
def brickwork_ops(n: int, layers: int, hf_index: int):
    ops = [("x", q) for q in range(n) if (hf_index >> q) & 1]
    k = 0
    for q in range(n):
        ops.append(("ry", q, k)); k += 1
    for layer in range(layers):
        for q in range(layer % 2, n - 1, 2):
            ops.append(("cz", q, q + 1))
        for q in range(n):
            ops.append(("ry", q, k)); k += 1
    return ops, k


class FastBrickwork:
    """Simulador numpy de estados reales (RY, CZ, X) con gradiente adjunto. Solo para OPTIMIZAR."""

    def __init__(self, n: int, ops, n_params: int, hf_index: int):
        self.n, self.ops, self.n_params = n, ops, n_params
        self.idx = np.arange(2 ** n)
        self.bit = [(self.idx >> q) & 1 for q in range(n)]
        self.cz_sign = {}
        self.psi0 = np.zeros(2 ** n, dtype=complex)
        self.psi0[hf_index] = 1.0

    def _apply(self, v, op, x, inv=False, deriv=False):
        if op[0] == "x":
            return v[self.idx ^ (1 << op[1])]
        if op[0] == "cz":
            key = (op[1], op[2])
            if key not in self.cz_sign:
                self.cz_sign[key] = np.where(self.bit[op[1]] & self.bit[op[2]], -1.0, 1.0)
            return v * self.cz_sign[key]
        q, t = op[1], (-x[op[2]] if inv else x[op[2]])
        c, s = np.cos(t / 2), np.sin(t / 2)
        if deriv:                                        # d/dt RY(t) = RY(t + π) / 2
            c, s = -0.5 * s, 0.5 * c
        w = v.reshape(2 ** (self.n - q - 1), 2, 2 ** q)
        a0, a1 = w[:, 0, :], w[:, 1, :]
        return np.stack([c * a0 - s * a1, s * a0 + c * a1], axis=1).reshape(-1)

    def _state_real(self, x):
        v = np.zeros(2 ** self.n)
        v[0] = 1.0
        for op in self.ops:
            v = self._apply(v, op, x)
        return v

    def state(self, x) -> np.ndarray:
        return self._state_real(np.asarray(x, float)).astype(complex)

    def energy_and_grad(self, x, Hm):
        """Hm: matriz (densa o sparse) del observable de 5 qubits. Gradiente exacto en 3 pasadas."""
        x = np.asarray(x, float)
        phi = self._state_real(x)
        lam = np.real(Hm @ phi)
        e = float(phi @ lam)
        g = np.zeros(self.n_params)
        for op in reversed(self.ops):
            phi = self._apply(phi, op, x, inv=True)
            if op[0] == "ry":
                g[op[2]] += 2.0 * float(lam @ self._apply(phi, op, x, deriv=True))
            lam = self._apply(lam, op, x, inv=True)
        return e, g


def build_compact_ansatz(m: dict, layers: int = 8) -> Ansatz:
    tp = tapering(m)
    n = tp.n_qubits
    ops, n_params = brickwork_ops(n, layers, tp.hf_index)
    theta = ParameterVector("θ", n_params)
    qc = QuantumCircuit(n)
    for op in ops:
        if op[0] == "x":
            qc.x(op[1])
        elif op[0] == "ry":
            qc.ry(theta[op[2]], op[1])
        else:
            qc.cz(op[1], op[2])
    native = transpile(qc, basis_gates=BASIS_GATES, optimization_level=1, seed_transpiler=SEED)
    return Ansatz(native, [], FastBrickwork(n, ops, n_params, tp.hf_index), theta,
                  {"kind": "brickwork-5q", "layers": layers, "n_qubits": n, "logical": qc,
                   "removed_qubits": tp.removed, "sector": tp.sector})


def optimize_compact(ans: Ansatz, H5: SparsePauliOp, S2_5: SparsePauliOp | None = None,
                     mu: float = 0.5, x0=None, n_random: int = 3, seed: int = SEED,
                     maxiter: int = 5000) -> dict:
    """Busca θ* SIN ruido (lo permite la consigna). Minimiza <H> + mu <S^2>: la penalización solo
    guía la búsqueda hacia el singlete (a 3 Å el triplete está a 0,95 mHa); la energía que se
    devuelve es <H> sin penalización. Arranca del θ anterior (warm start) y de n_random puntos al azar."""
    Hm = H5.to_matrix(sparse=True)
    C = Hm if (S2_5 is None or mu == 0) else Hm + mu * S2_5.to_matrix(sparse=True)
    rng = np.random.default_rng(seed)
    n = ans.num_parameters
    starts = [np.zeros(n) if x0 is None else np.asarray(x0, float)]
    starts += [rng.uniform(-np.pi, np.pi, n) for _ in range(n_random)]
    best = None
    for s in starts:
        r = minimize(lambda x: ans.fast.energy_and_grad(x, C), s, jac=True, method="L-BFGS-B",
                     options={"maxiter": maxiter, "gtol": 1e-10, "ftol": 1e-15})
        if best is None or r.fun < best.fun:
            best = r
    e = ans.fast.energy_and_grad(best.x, Hm)[0]
    return {"x": best.x, "energy": float(e), "nfev": int(best.nfev), "success": bool(best.success)}


# --------------------------------------------------------------------------- tablas para el notebook
def reduction_table(distances) -> list[dict]:
    """Por distancia: términos de Pauli y grupos QWC en 8 y 5 qubits, y chequeo de H₅ contra PySCF
    (fundamental de H₅ = E_FCI y <HF|H₅|HF> = E_HF)."""
    rows = []
    for d in distances:
        m = molecule(d)
        H8 = qubit_hamiltonian(m)
        tp = tapering(m, H8)
        H5 = tp.apply(H8)
        M5 = H5.to_matrix(sparse=True)
        rows.append({"d": d, "paulis_8q": len(H8), "paulis_5q": len(H5),
                     "groups_8q": len(measurement_groups(H8)), "groups_5q": len(measurement_groups(H5)),
                     "err_fci": float(np.linalg.eigvalsh(M5.toarray())[0] - m["e_fci"]),
                     "err_hf": float(M5[tp.hf_index, tp.hf_index].real - m["e_hf"])})
    return rows


def compare_ansatze(d: float = 2.0, layers: int = 8) -> list[dict]:
    """UCCSD en 8 qubits contra el ansatz compacto en 5: qubits, parámetros, CZ, grupos y error sin ruido."""
    m = molecule(d)
    H8 = qubit_hamiltonian(m)
    H5 = tapering(m, H8).apply(H8)
    u = build_ansatz(m, reps=ANSATZ_REPS)
    eu = optimize_params(u, H8, n_restarts=2)["energy"]
    c = build_compact_ansatz(m, layers=layers)
    ec = optimize_compact(c, H5, reduced_number_operators(m)["S2"])["energy"]
    return [{"ansatz": name, "qubits": a.circuit.num_qubits, "parámetros": a.num_parameters,
             "CZ (todos con todos)": cz_count(a.circuit), "CZ (línea)": cz_count(a.circuit, "line"),
             "grupos QWC": len(measurement_groups(H)), "error sin ruido (mHa)": 1e3 * (e - m["e_fci"])}
            for name, a, H, e in [(f"UCCSD ×{ANSATZ_REPS}", u, H8, eu), (f"compacto L={layers}", c, H5, ec)]]
