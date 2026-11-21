"""Extra: otra geometría de H4 (cuadrado y zigzag) con el mismo flujo que la cadena recta.

Dueño: Octavio.

Solo cambia la geometría (hamiltonian.h4_geometry). El hamiltoniano, el ansatz, la búsqueda de θ*
y la energía medida con shots son los mismos de la curva principal.
"""
from __future__ import annotations

from .ansatz import build_ansatz, cz_count
from .config import ANSATZ_REPS, SEED
from .estimator import measure_energy
from .hamiltonian import hf_index, molecule, qubit_hamiltonian, sector_ground_state
from .optimize import optimize_params


def scan_geometry(shape: str, distances, shots: int, reps: int = ANSATZ_REPS,
                  seed: int = SEED, verbose: bool = True) -> list[dict]:
    """VQE (energía medida ± σ) contra HF, CCSD y FCI para una geometría, en varias distancias."""
    rows, x_prev, cz = [], None, None
    for d in distances:
        m = molecule(d, shape=shape)
        H = qubit_hamiltonian(m)
        ans = build_ansatz(m, reps=reps)
        if cz is None:   # el circuito tiene la misma estructura en todas las distancias
            cz = {"cz_all": cz_count(ans.circuit), "cz_line": cz_count(ans.circuit, "line")}
        if x_prev is None:   # primer punto: desde HF con reinicios al azar
            opt = optimize_params(ans, H, n_restarts=2)
        else:                # resto: desde el θ* del punto anterior y desde HF; el mejor
            opt = min(optimize_params(ans, H, x0=x_prev), optimize_params(ans, H),
                      key=lambda o: o["energy"])
        x_prev = opt["x"]
        meas = measure_energy(ans.circuit.assign_parameters(opt["x"]), H, shots, seed=seed)
        _, psi = sector_ground_state(H, m["norb"], m["nelec"])
        row = {"shape": shape, "d": d, "E_vqe": meas["energy"], "sigma": meas["sigma"],
               "shots": meas["shots"], "E_opt_exact": opt["energy"],
               "E_hf": m["e_hf"], "E_fci": m["e_fci"], "E_ccsd": m["e_ccsd"],
               "ccsd_converged": m["ccsd_converged"],
               "w_hf": abs(psi[hf_index(m["norb"], m["nelec"])]) ** 2, **cz}
        row["err_vqe_mHa"] = 1e3 * (row["E_vqe"] - row["E_fci"])
        row["err_ccsd_mHa"] = 1e3 * (row["E_ccsd"] - row["E_fci"])
        rows.append(row)
        if verbose:
            print(f"{shape:6s} d={d:4.2f} Å  VQE−FCI = {row['err_vqe_mHa']:+6.2f} ± {1e3*row['sigma']:.2f} mHa  "
                  f"CCSD−FCI = {row['err_ccsd_mHa']:+7.2f} mHa{'' if row['ccsd_converged'] else ' (CCSD no convergió)'}")
    return rows
