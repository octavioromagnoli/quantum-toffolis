"""Genera notebook/H4_VQE.ipynb a partir de h4vqe/*.py y notebook/md/*.md.

    python notebook/build_notebook.py            # genera
    python notebook/build_notebook.py --run      # genera y ejecuta (verifica de punta a punta)

No editar el .ipynb a mano: se pisa. Editar los módulos o los .md y regenerar.
"""
from __future__ import annotations

import pathlib
import subprocess
import sys

import nbformat as nbf

ROOT = pathlib.Path(__file__).resolve().parent.parent
MD = ROOT / "notebook" / "md"
PKG = ROOT / "h4vqe"
OUT = ROOT / "notebook" / "H4_VQE.ipynb"

MODULES = ["__init__", "config", "hamiltonian", "ansatz", "optimize", "estimator",
           "analysis", "pipeline", "geometry", "reduction", "noise", "hardware"]

INSTALL = '''\
# Instala dependencias si faltan (Colab arranca sin Qiskit ni PySCF).
import importlib.util, os, subprocess, sys
REQS = {"qiskit": "qiskit==2.6.2", "qiskit_aer": "qiskit-aer==0.18.4", "pyscf": "pyscf==2.15.0",
        "qiskit_ibm_runtime": "qiskit-ibm-runtime==0.52.1"}   # runtime: solo para FakeFez
missing = [pip for mod, pip in REQS.items() if importlib.util.find_spec(mod) is None]
if missing:
    subprocess.check_call([sys.executable, "-m", "pip", "install", "-q", *missing])
os.makedirs("h4vqe", exist_ok=True)
print("python", sys.version.split()[0], "- instalado:", missing or "nada (ya estaba)")'''

HAMILTONIAN = '''\
import pandas as pd, matplotlib.pyplot as plt
from h4vqe import config
from h4vqe.hamiltonian import validation_table
val = pd.DataFrame(validation_table(config.DISTANCES))
print(f"{len(val)} distancias | máx |<HF|H|HF> − E_HF| = {val.err_hf.abs().max():.1e} Ha | "
      f"máx |E_0(H, sector 2α2β) − E_FCI| = {val.err_fci.abs().max():.1e} Ha | máx |S²| = {val.s2.abs().max():.1e}")
worst = val.loc[val.ccsd_minus_fci_mHa.idxmin()]
print(f"CCSD convergió en todas: {val.ccsd_converged.all()} | peor CCSD − FCI = {worst.ccsd_minus_fci_mHa:.1f} mHa a {worst.d} Å")
fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 3.6))
for a in (a1, a2):
    a.axvspan(config.STRETCHED_FROM, val.d.max() + 0.1, color="0.93", zorder=0)
    a.set_xlabel("d (Å)"); a.set_xlim(val.d.min() - 0.1, val.d.max() + 0.1)
a1.axhspan(-1e3 * config.CHEM_ACC, 1e3 * config.CHEM_ACC, color="tab:green", alpha=0.25, label="±1,6 mHa (precisión química)")
a1.axhline(0, color="k", lw=0.8)
a1.plot(val.d, val.ccsd_minus_fci_mHa, "o-", label="CCSD − FCI")
a1.set_ylabel("E − E_FCI (mHa)"); a1.legend(loc="lower left")
a2.plot(val.d, val.w_hf, "o-", color="tab:red")
a2.set_ylabel(r"peso de HF en FCI  $|\\langle \\mathrm{HF}|\\Psi_\\mathrm{FCI}\\rangle|^2$"); a2.set_ylim(0, 1)
plt.tight_layout(); plt.show()'''

GEOMETRY = '''\
import pandas as pd, matplotlib.pyplot as plt
from h4vqe import config
from h4vqe.geometry import scan_geometry
from h4vqe.hamiltonian import ZIGZAG_ANGLE
GEO_D = [0.7, 1.2, 1.8, 2.4, 3.0]
GEO_SHOTS = 2_000_000   # σ ≈ 1,4 mHa por punto
geo = pd.DataFrame(scan_geometry("square", GEO_D, GEO_SHOTS) + scan_geometry("zigzag", GEO_D, GEO_SHOTS))
geo["beats_ccsd"] = geo.err_vqe_mHa.abs() + 1e3 * geo.sigma < geo.err_ccsd_mHa.abs()
fig, axs = plt.subplots(1, 2, figsize=(11, 3.6), sharey=True)
for a, shape, title in zip(axs, ["square", "zigzag"], ["H$_4$ cuadrado (lado d)", f"H$_4$ en zigzag ({ZIGZAG_ANGLE:.0f}°)"]):
    g = geo[geo["shape"] == shape]
    a.axhspan(-1e3 * config.CHEM_ACC, 1e3 * config.CHEM_ACC, color="tab:green", alpha=0.25, label="±1,6 mHa")
    a.axhline(0, color="k", lw=0.8)
    a.errorbar(g.d, g.err_vqe_mHa, yerr=1e3 * g.sigma, fmt="o-", capsize=3, label="VQE − FCI (medido, ±σ)")
    a.plot(g.d, g.err_ccsd_mHa, "s--", color="tab:orange", label="CCSD − FCI")
    bad = g[~g.ccsd_converged]
    a.plot(bad.d, bad.err_ccsd_mHa, "x", color="k", ms=10, mew=2, label="CCSD no convergió")
    a.set_title(title); a.set_xlabel("d (Å)")
axs[0].set_ylabel("E − E_FCI (mHa)"); axs[0].legend(loc="lower left", fontsize=8)
plt.tight_layout(); plt.show()
geo[["shape", "d", "E_vqe", "sigma", "err_vqe_mHa", "err_ccsd_mHa", "ccsd_converged", "beats_ccsd", "cz_all", "cz_line"]]'''

HARDWARE = '''\
import pandas as pd
from h4vqe.hardware import run_fez
fz = run_fez(2.4)
print(f"línea de Fez: qubits {fz['line']} | errores esperados {fz['expected_errors']:.3f} "
      f"(mediana de las {fz['n_lines']} líneas posibles: {fz['median_expected_errors']:.3f})")
print(f"CZ: compacto {fz['cz']} (sin SWAP) | UCCSD ×2 ubicado en Fez: {fz['cz_uccsd_fez']} | "
      f"error del compacto sin ruido: {fz['err_exact_mHa']:.1e} mHa")
print("lectura por qubit, P(1|0) / P(0|1): " + ", ".join(f"{a:.2%} / {b:.2%}" for a, b in fz["readout"]))
pd.DataFrame(fz["rows"]).astype({"shots": "Int64"})'''

CHECKS = '''\
import numpy as np, importlib, time
import h4vqe
for mod in ["config","hamiltonian","ansatz","optimize","estimator","analysis","pipeline"]:
    importlib.reload(importlib.import_module("h4vqe." + mod))
from qiskit.quantum_info import Statevector
from h4vqe import config
from h4vqe.hamiltonian import molecule, qubit_hamiltonian
from h4vqe.ansatz import build_ansatz, cz_count

m = molecule(2.0); H = qubit_hamiltonian(m)
ans = build_ansatz(m, reps=config.ANSATZ_REPS)
e_hf = Statevector(ans.fast.psi0).expectation_value(H).real
print(f"términos de Pauli: {len(H)} | qubits: {H.num_qubits}")
print(f"E_HF (PySCF) = {m['e_hf']:.8f}   <HF|H|HF> = {e_hf:.8f}")
print(f"E_FCI = {m['e_fci']:.8f}   E_CCSD = {m['e_ccsd']:.8f}  (CCSD − FCI = {1e3*(m['e_ccsd']-m['e_fci']):+.2f} mHa)")
x = np.random.default_rng(0).normal(0, .2, ans.num_parameters)
fid = abs(np.vdot(Statevector(ans.circuit.assign_parameters(x)).data, ans.fast.state(x)))**2
print(f"fidelidad circuito Qiskit vs simulador numpy = {fid:.12f}")
print(f"parámetros: {ans.num_parameters} | CZ: {cz_count(ans.circuit)} (todos con todos), {cz_count(ans.circuit, 'line')} (línea)")'''

RUN = '''\
from h4vqe.pipeline import run_curve
QUICK = False   # True: 4 distancias y menos shots, para probar rápido
dists = [0.7, 0.9, 1.8, 3.0] if QUICK else config.DISTANCES
shots = 400_000 if QUICK else config.SHOTS_PER_POINT
t0 = time.time()
rows, ctxs = run_curve(dists, shots=shots)
print(f"total: {time.time()-t0:.0f} s")'''

TABLE = '''\
import pandas as pd
from h4vqe.analysis import results_table
df = results_table(rows)
cols = ["d","E_vqe","sigma_mHa","shots","cz_all","cz_line","E_fci","E_ccsd","err_vqe_mHa","err_ccsd_mHa","beats_ccsd","chem_acc"]
pd.set_option("display.float_format", lambda v: f"{v:.6f}")
df[cols]'''

PLOT = '''\
from h4vqe.analysis import plot_curve
fig = plot_curve(df, "curva_h4.png")'''

RE_DE = '''\
from h4vqe.analysis import fit_minimum, dissociation_energy, fit_minimum as _fm
fit = fit_minimum(df.d, df.E_vqe, df.sigma)
far = df[df.d == df.d.max()].iloc[0]
de = dissociation_energy(far.E_vqe, far.sigma, fit)
fit_fci = _fm(df.d, df.E_fci, np.full(len(df), 1e-6))
print(f"R_e (VQE) = {fit['R_e']:.4f} ± {fit['R_e_sigma']:.4f} Å   | R_e (FCI, mismo ajuste) = {fit_fci['R_e']:.4f} Å")
print(f"D_e (VQE) = {1e3*de['D_e']:.2f} ± {1e3*de['D_e_sigma']:.2f} mHa | D_e (FCI) = {1e3*(far.E_fci - fit_fci['E_min']):.2f} mHa")
s = df[df.d >= config.STRETCHED_FROM]
print(f"zona estirada: le gana a CCSD en {int(s.beats_ccsd.sum())}/{len(s)} puntos; precisión química en {int(s.chem_acc.sum())}/{len(s)}")'''

RELIABILITY = '''\
from h4vqe.analysis import variance
from h4vqe.hamiltonian import number_operators
from h4vqe.optimize import optimize_params
diag = []
for d, c in ctxs.items():
    if d < config.STRETCHED_FROM: continue
    Hm = c["H"].to_matrix(sparse=True); psi = c["ansatz"].fast.state(c["opt"]["x"])
    ops = number_operators(c["m"]["norb"]); sv = Statevector(psi)
    a1 = build_ansatz(c["m"], reps=1); o1 = optimize_params(a1, c["H"], n_restarts=2)
    diag.append({"d": d, "var_H": variance(psi, Hm), "N_alpha": sv.expectation_value(ops["N_alpha"]).real,
                 "N_beta": sv.expectation_value(ops["N_beta"]).real, "S2": sv.expectation_value(ops["S2"]).real,
                 "E(reps=1) - E(reps=2) [mHa]": 1e3*(o1["energy"] - c["opt"]["energy"])})
pd.DataFrame(diag)'''


def md(name: str):
    return nbf.v4.new_markdown_cell((MD / name).read_text())


def build() -> nbf.NotebookNode:
    nb = nbf.v4.new_notebook()
    cells = [md("00_intro.md"), md("01_setup.md"), nbf.v4.new_code_cell(INSTALL)]
    for mod in MODULES:
        src = (PKG / f"{mod}.py").read_text()
        cells.append(nbf.v4.new_code_cell(f"%%writefile h4vqe/{mod}.py\n{src}"))
    cells += [md("02_hamiltonian.md"), nbf.v4.new_code_cell(HAMILTONIAN),
              md("03_ansatz.md"), md("04_measurement.md"),
              nbf.v4.new_code_cell(CHECKS),
              md("05_results.md"), nbf.v4.new_code_cell(RUN), nbf.v4.new_code_cell(TABLE),
              nbf.v4.new_code_cell(PLOT), nbf.v4.new_code_cell(RE_DE),
              md("06_reliability.md"), nbf.v4.new_code_cell(RELIABILITY),
              md("07_extras.md"), md("07_geometry.md"), nbf.v4.new_code_cell(GEOMETRY),
              md("07_hardware.md"), nbf.v4.new_code_cell(HARDWARE),
              md("08_decisions.md")]
    nb.cells = cells
    nb.metadata["kernelspec"] = {"name": "python3", "display_name": "Python 3", "language": "python"}
    nb.metadata["colab"] = {"provenance": []}
    return nb


if __name__ == "__main__":
    nbf.write(build(), OUT)
    print("escrito", OUT)
    if "--run" in sys.argv:
        subprocess.check_call([sys.executable, "-m", "jupyter", "nbconvert", "--to", "notebook",
                               "--execute", "--inplace", "--ExecutePreprocessor.timeout=-1",
                               "--ExecutePreprocessor.kernel_name=h4vqe-local", str(OUT)])
        print("ejecutado OK")
