"""Genera notebook/H4_VQE.ipynb a partir de h4vqe/*.py y notebook/md/*.md.

    python notebook/build_notebook.py            # genera
    python notebook/build_notebook.py --run      # genera y ejecuta (verifica de punta a punta)
    python notebook/build_notebook.py --quick --run   # igual, con QUICK = True (antes de cada PR)

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
           "analysis", "pipeline"]

INSTALL = '''\
# Instala dependencias si faltan (Colab arranca sin Qiskit ni PySCF).
import importlib.util, os, subprocess, sys
REQS = {"qiskit": "qiskit==2.6.2", "qiskit_aer": "qiskit-aer==0.18.4", "pyscf": "pyscf==2.15.0"}
missing = [pip for mod, pip in REQS.items() if importlib.util.find_spec(mod) is None]
if missing:
    subprocess.check_call([sys.executable, "-m", "pip", "install", "-q", *missing])
os.makedirs("h4vqe", exist_ok=True)
print("python", sys.version.split()[0], "- instalado:", missing or "nada (ya estaba)")'''

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
from h4vqe.analysis import fit_minimum, dissociation_energy, fci_minimum, compare_fits, DEFAULT_FIT
# Referencia: mínimo FCI verdadero (minimizando PySCF en d), NO un ajuste sobre la grilla.
ref = fci_minimum()
far = df[df.d == df.d.max()].iloc[0]
De_fci = far.E_fci - ref["E_min"]
# Comparación de modelos. "sesgo" = modelo aplicado a FCI en la misma grilla − mínimo verdadero:
# mide el error del método de ajuste, sin mirar el VQE. Con eso elegimos el modelo.
cmp = compare_fits(df, ref)
display(cmp.round(3))
fit = fit_minimum(df.d, df.E_vqe, df.sigma, DEFAULT_FIT)
de = dissociation_energy(far.E_vqe, far.sigma, fit)
print(f"modelo: {fit['model']} sobre d = {fit['fit_points']}")
print(f"R_e (VQE) = {fit['R_e']:.4f} ± {fit['R_e_sigma']:.4f} Å   | R_e (FCI) = {ref['R_e']:.4f} Å   "
      f"→ diferencia {1e3*(fit['R_e']-ref['R_e']):+.1f} mÅ = {(fit['R_e']-ref['R_e'])/fit['R_e_sigma']:+.1f} σ")
print(f"D_e (VQE) = {1e3*de['D_e']:.2f} ± {1e3*de['D_e_sigma']:.2f} mHa | D_e (FCI) = {1e3*De_fci:.2f} mHa  "
      f"→ diferencia {1e3*(de['D_e']-De_fci):+.2f} mHa = {(de['D_e']-De_fci)/de['D_e_sigma']:+.1f} σ")
s = df[df.d >= config.STRETCHED_FROM]
print(f"zona estirada: le gana a CCSD en {int(s.beats_ccsd.sum())}/{len(s)} puntos; precisión química en {int(s.chem_acc.sum())}/{len(s)}")'''

RELIABILITY = '''\
from h4vqe.pipeline import run_reliability
# Tres puntos de la zona estirada; en cada uno, HF vs UCCSD ×1 vs el ansatz de la curva.
rel_d = [1.8, 3.0] if QUICK else [1.8, 2.4, 3.0]
t0 = time.time()
rel = pd.DataFrame(run_reliability(ctxs, rel_d, shots=shots,
                                   var_shots=500_000 if QUICK else 2_000_000,
                                   echo_shots=200_000 if QUICK else 1_000_000,
                                   sym_shots=100_000 if QUICK else 200_000))
print(f"total: {time.time()-t0:.0f} s")
show = rel[["d","estado","cz","E","sigma","var","var_sigma","var_echo","var_echo_up95","echo_events",
            "N_alpha","N_alpha_sigma","N_beta","N_beta_sigma","S2","S2_sigma","var_exacta"]]
with pd.option_context("display.float_format", lambda v: f"{v:.3g}"):
    display(show)'''

DEPTH = '''\
# Convergencia con la profundidad (sin FCI): cuánto baja la energía MEDIDA al agregar
# parámetros, y cuánto cae la varianza (eco) en el mismo paso.
piv = rel.pivot(index="d", columns="estado", values="E")
lbl = [c for c in piv.columns if c.startswith("curva")][0]
lbl1 = [c for c in piv.columns if c.startswith("UCCSD")][0]
sig = rel.pivot(index="d", columns="estado", values="sigma")
conv = pd.DataFrame({
    "E(×1) − E(curva) [mHa]": 1e3*(piv[lbl1] - piv[lbl]),
    "σ [mHa]": 1e3*np.hypot(sig[lbl1], sig[lbl]),
    "Var eco ×1 [Ha²]": rel[rel.estado == lbl1].set_index("d").var_echo,
    "Var eco curva ≤ [Ha²] (95 %)": rel[rel.estado == lbl].set_index("d").var_echo_up95,
})
with pd.option_context("display.float_format", lambda v: f"{v:.3g}"):
    display(conv)
# Recién ahora, como comprobación: el error real contra FCI.
chk = df.set_index("d").loc[rel_d, ["err_vqe_mHa", "sigma_mHa", "err_ccsd_mHa"]]
chk'''


def md(name: str):
    return nbf.v4.new_markdown_cell((MD / name).read_text())


def build(quick: bool = False) -> nbf.NotebookNode:
    nb = nbf.v4.new_notebook()
    run = RUN.replace("QUICK = False", "QUICK = True ") if quick else RUN
    cells = [md("00_intro.md"), md("01_setup.md"), nbf.v4.new_code_cell(INSTALL)]
    for mod in MODULES:
        src = (PKG / f"{mod}.py").read_text()
        cells.append(nbf.v4.new_code_cell(f"%%writefile h4vqe/{mod}.py\n{src}"))
    cells += [md("02_hamiltonian.md"), md("03_ansatz.md"), md("04_measurement.md"),
              nbf.v4.new_code_cell(CHECKS),
              md("05_results.md"), nbf.v4.new_code_cell(run), nbf.v4.new_code_cell(TABLE),
              nbf.v4.new_code_cell(PLOT), nbf.v4.new_code_cell(RE_DE),
              md("06_reliability.md"), nbf.v4.new_code_cell(RELIABILITY), nbf.v4.new_code_cell(DEPTH),
              md("07_extras.md"), md("08_decisions.md")]
    nb.cells = cells
    nb.metadata["kernelspec"] = {"name": "python3", "display_name": "Python 3", "language": "python"}
    nb.metadata["colab"] = {"provenance": []}
    return nb


if __name__ == "__main__":
    nbf.write(build(quick="--quick" in sys.argv), OUT)
    print("escrito", OUT)
    if "--run" in sys.argv:
        subprocess.check_call([sys.executable, "-m", "jupyter", "nbconvert", "--to", "notebook",
                               "--execute", "--inplace", "--ExecutePreprocessor.timeout=-1",
                               "--ExecutePreprocessor.kernel_name=h4vqe-local", str(OUT)])
        print("ejecutado OK")
