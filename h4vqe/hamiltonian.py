"""Molécula H4 lineal -> hamiltoniano de qubits + referencias clásicas.

Dueño: Octavio.

Convención de qubits (Jordan-Wigner, little-endian de Qiskit):
    qubits 0..norb-1      -> espín-orbitales alfa (orbital 0 = el de menor energía)
    qubits norb..2norb-1  -> espín-orbitales beta
El estado HF ocupa los qubits 0,1 (alfa) y 4,5 (beta): bitstring '00110011'.
"""
from __future__ import annotations

import numpy as np
from pyscf import gto, scf, mcscf, ao2mo, cc
from qiskit.quantum_info import SparsePauliOp

from .config import BASIS


ZIGZAG_ANGLE = 120.0   # ángulo H-H-H de la cadena en zigzag (grados)


def h4_geometry(d: float, shape: str = "linear") -> str:
    """4 H con vecinos a d Å: cadena recta (la del desafío), cadena en zigzag o cuadrado de lado d."""
    if shape == "linear":
        xyz = [(0, 0, i * d) for i in range(4)]
    elif shape == "zigzag":
        a = np.radians(ZIGZAG_ANGLE / 2)
        xyz = [(i * d * np.sin(a), (i % 2) * d * np.cos(a), 0) for i in range(4)]
    elif shape == "square":
        xyz = [(0, 0, 0), (d, 0, 0), (d, d, 0), (0, d, 0)]
    else:
        raise ValueError(f"forma desconocida: {shape!r}")
    return "; ".join(f"H {x:.6f} {y:.6f} {z:.6f}" for x, y, z in xyz)


def stable_rhf(mol, max_iter: int = 10):
    """RHF que sigue las inestabilidades internas hasta llegar a una solución estable.

    En el H4 cuadrado hay orbitales degenerados y el SCF puede terminar en un punto silla, según el
    entorno (hilos, versión de Python): a 1,2 Å queda 76 mHa por encima del RHF estable, y con esa
    referencia CCSD daba +180 mHa en vez de −6,5 mHa. En la cadena recta y en el zigzag RHF ya es
    estable en las 18 distancias, así que ahí esto no cambia nada.
    """
    hf = scf.RHF(mol)
    hf.conv_tol = 1e-10
    hf.kernel()
    for _ in range(max_iter):
        mo, _, stable, _ = hf.stability(return_status=True)
        if stable:
            break
        hf.kernel(dm0=hf.make_rdm1(mo, hf.mo_occ))
    return hf


def molecule(d: float, basis: str = BASIS, shape: str = "linear") -> dict:
    """Integrales en base de orbitales HF + energías HF, FCI y CCSD."""
    mol = gto.M(atom=h4_geometry(d, shape), basis=basis, unit="Angstrom", verbose=0)
    hf = stable_rhf(mol)
    norb = mol.nao_nr()
    nelec = mol.nelectron
    cas = mcscf.CASCI(hf, norb, nelec)
    h1, ecore = cas.get_h1cas()
    g = ao2mo.restore(1, cas.get_h2cas(), norb)
    e_fci = cas.kernel()[0]

    ccsd = cc.CCSD(hf)
    ccsd.conv_tol = 1e-9
    ccsd.max_cycle = 500
    ccsd.kernel()

    return {
        "d": d, "shape": shape, "mol": mol, "hf": hf, "norb": norb,
        "nelec": (nelec // 2, nelec // 2),
        "h": h1, "g": g, "constante": ecore,
        "e_hf": hf.e_tot, "e_fci": e_fci,
        "e_ccsd": ccsd.e_tot, "ccsd_converged": bool(ccsd.converged),
    }


def _annihilation_jw(p: int, n_qubits: int) -> SparsePauliOp:
    return SparsePauliOp.from_sparse_list(
        [("Z" * p + "X", list(range(p + 1)), 0.5),
         ("Z" * p + "Y", list(range(p + 1)), 0.5j)], n_qubits)


def annihilation_ops(n_qubits: int) -> list[SparsePauliOp]:
    return [_annihilation_jw(p, n_qubits) for p in range(n_qubits)]


def qubit_hamiltonian(m: dict) -> SparsePauliOp:
    """Ĥ = E_core + Σ h_pq a†_p a_q + ½ Σ g_pqrs a†_p a†_r a_s a_q  (Jordan-Wigner)."""
    norb = m["norb"]
    nq = 2 * norb
    h, g = m["h"], m["g"]
    a = annihilation_ops(nq)
    ad = [op.adjoint() for op in a]
    hop = [[[ad[s * norb + p].dot(a[s * norb + q]).simplify() for q in range(norb)]
            for p in range(norb)] for s in range(2)]
    terms = [SparsePauliOp("I" * nq, m["constante"])]
    for s in range(2):
        for p in range(norb):
            for q in range(norb):
                terms.append(h[p, q] * hop[s][p][q])
    for s in range(2):
        for t in range(2):
            for p in range(norb):
                for q in range(norb):
                    second = SparsePauliOp.sum([g[p, q, r, u] * hop[t][r][u]
                                                for r in range(norb) for u in range(norb)])
                    terms.append(0.5 * hop[s][p][q].dot(second))
                    if s == t:
                        for u in range(norb):
                            terms.append(-0.5 * g[p, q, q, u] * hop[s][p][u])
    return SparsePauliOp.sum(terms).simplify(atol=1e-10).chop(1e-10)


def number_operators(norb: int) -> dict[str, SparsePauliOp]:
    """N_alfa, N_beta y S^2 en qubits, para chequear simetrías del estado."""
    nq = 2 * norb
    a = annihilation_ops(nq)
    n = [a[p].adjoint().dot(a[p]).simplify() for p in range(nq)]
    na = SparsePauliOp.sum(n[:norb]).simplify()
    nb = SparsePauliOp.sum(n[norb:]).simplify()
    sz = 0.5 * (na - nb)
    sp = SparsePauliOp.sum([a[p].adjoint().dot(a[norb + p]) for p in range(norb)]).simplify()
    s2 = (sp.dot(sp.adjoint()) + sz.dot(sz) - sz).simplify(atol=1e-12)
    return {"N_alpha": na, "N_beta": nb, "S2": s2}


def exact_ground_energy(H: SparsePauliOp) -> float:
    """Diagonalización completa (256x256). Solo como chequeo contra FCI."""
    return float(np.linalg.eigvalsh(H.to_matrix())[0])


def sector_indices(norb: int, nelec: tuple[int, int]) -> np.ndarray:
    """Índices de la base computacional con N_alfa y N_beta fijos (bit q = ocupación del qubit q)."""
    i = np.arange(2 ** (2 * norb))
    bits = (i[:, None] >> np.arange(2 * norb)) & 1
    return i[(bits[:, :norb].sum(1) == nelec[0]) & (bits[:, norb:].sum(1) == nelec[1])]


def hf_index(norb: int, nelec: tuple[int, int]) -> int:
    """Índice del determinante HF: los n_alfa y n_beta orbitales más bajos ocupados."""
    return sum(1 << q for q in range(nelec[0])) + sum(1 << (norb + q) for q in range(nelec[1]))


def sector_ground_state(H: SparsePauliOp, norb: int, nelec: tuple[int, int]) -> tuple[float, np.ndarray]:
    """Estado fundamental de Ĥ solo en el sector (N_alfa, N_beta): 36x36 para H4 en vez de 256x256."""
    idx = sector_indices(norb, nelec)
    w, v = np.linalg.eigh(H.to_matrix(sparse=True)[idx][:, idx].toarray())
    psi = np.zeros(2 ** (2 * norb), dtype=complex)
    psi[idx] = v[:, 0]
    return float(w[0]), psi


def validation_table(distances) -> list[dict]:
    """Por distancia: chequeo de Ĥ contra PySCF (HF y FCI), peso de HF en FCI y error de CCSD."""
    rows = []
    s2_op = None
    for d in distances:
        m = molecule(d)
        H = qubit_hamiltonian(m)
        k = hf_index(m["norb"], m["nelec"])
        e_hf_q = H.to_matrix(sparse=True)[k, k].real
        e_sec, psi = sector_ground_state(H, m["norb"], m["nelec"])
        if s2_op is None:
            s2_op = number_operators(m["norb"])["S2"].to_matrix(sparse=True)
        rows.append({
            "d": d, "n_paulis": len(H),
            "err_hf": e_hf_q - m["e_hf"],             # <HF|Ĥ|HF> − E_HF(PySCF)
            "err_fci": e_sec - m["e_fci"],            # E_0(Ĥ, sector 2α2β) − E_FCI(PySCF)
            "s2": float(np.vdot(psi, s2_op @ psi).real),  # 0 = singulete
            "w_hf": abs(psi[k]) ** 2,                 # peso del determinante HF en el estado FCI
            "ccsd_minus_fci_mHa": 1e3 * (m["e_ccsd"] - m["e_fci"]),
            "ccsd_converged": m["ccsd_converged"],
        })
    return rows
