"""Reducción por simetrías (8 -> 6 o 5 qubits) y ansatz compacto nativo de la línea.

Dueño: Octavio.

Simetrías Z2 del H4 lineal en Jordan-Wigner (todas de tipo Z; conmutan con H, N, S_z y S^2):
    P_alfa = Z0 Z1 Z2 Z3                  paridad de N_alfa
    P_beta = Z4 Z5 Z6 Z7                  paridad de N_beta
    P_u    = Π Z sobre los orbitales impares ante la inversión de la cadena (σu), α y β
Tapering (Bravyi, Gambetta, Mezzacapo, Temme 2017): se lleva el conjunto de generadores a forma
escalonada sobre GF(2), así cada generador g_i tiene un qubit "exclusivo" q_i. Con
U_i = (X_{q_i} + g_i)/√2 se cumple U_i g_i U_i = X_{q_i}; en U H U el qubit q_i solo aparece como
I o X, y X_{q_i} se reemplaza por el autovalor de g_i en el sector del estado HF.
El estado HF reducido sigue siendo un bitstring (los bits HF de los qubits que quedan).

Dos niveles (`n_qubits`):
    6  solo P_alfa y P_beta (equivale al mapeo de paridad con la reducción de 2 qubits). Toda
       excitación conserva N_alfa y N_beta, así que los 26 generadores del pool de UCCSD/ADAPT se
       reducen con .apply(G). No usa la inversión: sirve también con orbitales rotados o localizados.
       Es la versión para el ansatz de excitaciones (ADAPT, Gonzalo).
    5  además P_u (default). Solo se reducen las excitaciones simétricas ante la inversión (14 de 26)
       y exige orbitales canónicos (con paridad definida). Es la del ansatz compacto de 16 CZ.

Ansatz compacto ("brickwork"): X en HF, una capa de RY y L capas de [CZ en ladrillo + RY].
CZ solo entre vecinos (qubits en línea: 0 SWAP). En 5 qubits cada capa tiene 2 CZ: con L = 8 son
16 CZ y 45 parámetros, y llega a FCI (< 0,001 mHa) en toda la curva. OJO: la versión con TODAS las
CZ en cada capa (cadena 0-1,1-2,2-3,3-4 simultánea) tiene un techo de fidelidad 0,961 que no mejora
con más capas.

Interfaz (la misma que ansatz.Ansatz, así optimize/pipeline/estimator funcionan igual):
    tapering(m, n_qubits=5) -> Tapering    .apply(op8) -> op reducido, .hf_index, .removed, .kept
    reduced_hamiltonian(m, n_qubits=5) -> SparsePauliOp
    reduced_number_operators(m, n_qubits=5) -> {"N_alpha","N_beta","S2"} reducidos
    build_compact_ansatz(m, layers=None, n_qubits=5) -> Ansatz  (layers=None: COMPACT_LAYERS)
    optimize_compact(ans, H, S2=None, mu=0.5, x0=None, n_random=3) -> dict como optimize_params
    neyman_cost(H, psi) -> (Σ_g √Var_g)²: shots · σ² con los grupos y el reparto de estimator.py
"""
from __future__ import annotations

import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.circuit import ParameterVector
from qiskit.quantum_info import SparsePauliOp
from scipy.optimize import minimize

from .ansatz import (Ansatz, build_ansatz, cz_count, excitation_generator, hf_occupation,
                     uccsd_excitations)
from .config import ANSATZ_REPS, BASIS_GATES, SEED
from .estimator import measurement_groups
from .hamiltonian import molecule, number_operators, qubit_hamiltonian, sector_ground_state
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


def z2_symmetries(m: dict, n_qubits: int = 5) -> list[SparsePauliOp]:
    """Generadores que se eliminan: las 2 paridades (-> 6 qubits) o también la inversión (-> 5)."""
    norb = m["norb"]
    nq = 2 * norb
    gens = [_zstring(range(norb), nq), _zstring(range(norb, nq), nq)]
    if n_qubits == nq - 2:
        return gens
    if n_qubits != nq - 3:
        raise ValueError(f"n_qubits tiene que ser {nq - 2} o {nq - 3}, no {n_qubits}")
    u = [p for p, s in enumerate(inversion_parities(m)) if s < 0]
    return gens + [_zstring(u + [norb + p for p in u], nq)]


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


# Qubits que se eliminan en cada nivel. 5: la elección validada con el ansatz compacto (ver
# docs/ideas/diferencial.md). 6: la de menos shots para la misma σ con el estimador de main
# (sorted insertion): ×1,02–1,07 de los de 8 qubits a 2,0–2,4 Å, contra ×1,21–1,27 con (0, 4).
PREFERRED_REMOVED = {5: (0, 4, 7), 6: (3, 5)}


class Tapering:
    def __init__(self, generators: list[SparsePauliOp], hf_index: int,
                 prefer_removed=PREFERRED_REMOVED[5]):
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
        """Operador de 8 qubits que conmuta con las simetrías -> operador reducido en el sector HF.
        Si no conmuta (por ejemplo, una excitación σg -> σu en 5 qubits) tira ValueError."""
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


def tapering(m: dict, H: SparsePauliOp | None = None, n_qubits: int = 5) -> Tapering:
    """Tapering a 6 (paridades) o 5 qubits (+ inversión). Orden de qubits: el de qubit_hamiltonian."""
    gens = z2_symmetries(m, n_qubits)
    if H is not None:
        for g in gens:
            assert _commutes(H, g), "una simetría no conmuta con H"
    hf_index = sum(1 << q for q in hf_occupation(m))
    return Tapering(gens, hf_index, prefer_removed=PREFERRED_REMOVED[n_qubits])


def reduced_hamiltonian(m: dict, H: SparsePauliOp | None = None, n_qubits: int = 5) -> SparsePauliOp:
    H = qubit_hamiltonian(m) if H is None else H
    return tapering(m, H, n_qubits).apply(H)


def reduced_number_operators(m: dict, n_qubits: int = 5) -> dict[str, SparsePauliOp]:
    tp = tapering(m, n_qubits=n_qubits)
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


COMPACT_LAYERS = {5: 8, 6: 10}   # capas que llegan a FCI en toda la curva (5q: 16 CZ; 6q: 25 CZ)


def build_compact_ansatz(m: dict, layers: int | None = None, n_qubits: int = 5) -> Ansatz:
    tp = tapering(m, n_qubits=n_qubits)
    n = tp.n_qubits
    layers = COMPACT_LAYERS[n_qubits] if layers is None else layers
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
                  {"kind": f"brickwork-{n}q", "layers": layers, "n_qubits": n, "logical": qc,
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
def neyman_cost(H: SparsePauliOp, psi: np.ndarray) -> float:
    """(Σ_g √Var_g)² en el estado psi, con los grupos de estimator.measurement_groups. Con el reparto
    de Neyman, σ² = costo / shots: para la misma σ, los shots son proporcionales a este número."""
    total = 0.0
    for g in measurement_groups(H):
        v = g.to_matrix(sparse=True) @ psi
        total += np.sqrt(max(float(np.vdot(v, v).real - np.vdot(psi, v).real ** 2), 0.0))
    return total ** 2


def reduction_table(distances) -> list[dict]:
    """Por distancia: términos de Pauli y grupos QWC en 8, 6 y 5 qubits, y chequeo de Ĥ₆ y Ĥ₅ contra
    PySCF (fundamental = E_FCI y <HF|Ĥ|HF> = E_HF). err_fci y err_hf: el peor de los dos."""
    rows = []
    for d in distances:
        m = molecule(d)
        H8 = qubit_hamiltonian(m)
        row = {"d": d, "paulis_8q": len(H8), "groups_8q": len(measurement_groups(H8)),
               "err_fci": 0.0, "err_hf": 0.0}
        for n in (6, 5):
            tp = tapering(m, H8, n)
            Hr = tp.apply(H8)
            M = Hr.to_matrix(sparse=True)
            e_fci = float(np.linalg.eigvalsh(M.toarray())[0] - m["e_fci"])
            e_hf = float(M[tp.hf_index, tp.hf_index].real - m["e_hf"])
            row[f"paulis_{n}q"], row[f"groups_{n}q"] = len(Hr), len(measurement_groups(Hr))
            row["err_fci"] = max(row["err_fci"], e_fci, key=abs)
            row["err_hf"] = max(row["err_hf"], e_hf, key=abs)
        rows.append(row)
    return rows


def levels_table(d: float = 2.0) -> list[dict]:
    """8 -> 6 -> 5 qubits en una distancia: términos, grupos, shots para la misma σ (relativos a 8
    qubits, en el estado exacto) y cuántas excitaciones del pool de UCCSD/ADAPT se pueden reducir."""
    m = molecule(d)
    H8 = qubit_hamiltonian(m)
    pool = [excitation_generator(e, H8.num_qubits) for e in uccsd_excitations(m)]
    _, psi8 = sector_ground_state(H8, m["norb"], m["nelec"])
    base = neyman_cost(H8, psi8)
    rows = [{"qubits": H8.num_qubits, "simetrías": "—", "términos de Pauli": len(H8),
             "grupos QWC": len(measurement_groups(H8)), "shots (misma σ)": 1.0,
             "excitaciones del pool": f"{len(pool)}/{len(pool)}"}]
    for n, sym in ((6, "N_α, N_β (paridad)"), (5, "+ inversión")):
        tp = tapering(m, H8, n)
        Hr = tp.apply(H8)
        psi = np.linalg.eigh(Hr.to_matrix())[1][:, 0]
        ok = 0
        for G in pool:
            try:
                tp.apply(G)
                ok += 1
            except ValueError:
                pass
        rows.append({"qubits": n, "simetrías": sym, "términos de Pauli": len(Hr),
                     "grupos QWC": len(measurement_groups(Hr)),
                     "shots (misma σ)": round(neyman_cost(Hr, psi) / base, 2),
                     "excitaciones del pool": f"{ok}/{len(pool)}"})
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
