"""Genera el notebook de entrega a partir de h4vqe/*.py y notebook/md/*.md.

    notebook/H4_VQE.ipynb   todo en un solo notebook, contado como historia: de lo mínimo que pide
                            la consigna (UCCSD) a la versión final (ADAPT, mediciones calibradas,
                            confiabilidad medida) y los extras (ruido, mitigación, FakeFez, otra
                            geometría, comparación de ansätze y escalado a FeMoco)

    python notebook/build_notebook.py                  # genera el notebook
    python notebook/build_notebook.py --run            # genera y ejecuta en una carpeta vacía
    python notebook/build_notebook.py --quick --run    # QUICK = True en la curva (antes de cada PR)
    python notebook/build_notebook.py --kernel NOMBRE  # kernel para ejecutar (default h4vqe-local)
    python notebook/build_notebook.py --text-only      # actualiza solo el texto del ya ejecutado

No editar el .ipynb a mano: se pisa. Editar los módulos o los .md y regenerar.
La ejecución en una carpeta temporal vacía (idea de Santino, build_minimo.py) reproduce un Colab
nuevo: %%writefile crea h4vqe/ ahí y no se usa el paquete del repo.
"""
from __future__ import annotations

import pathlib
import re
import shutil
import sys
import tempfile
import time

import nbformat as nbf

ROOT = pathlib.Path(__file__).resolve().parent.parent
MD = ROOT / "notebook" / "md"
PKG = ROOT / "h4vqe"
OUT = ROOT / "notebook" / "H4_VQE.ipynb"
SCALING_CSV = ROOT / "results" / "scaling" / "recursos.csv"

MODULES = ["__init__", "config", "hamiltonian", "ansatz", "optimize", "estimator", "adapt", "routing",
           "analysis", "pipeline", "reduction", "noise", "mitigation", "hardware", "geometry",
           "scaling", "historia"]

INSTALL = '''\
# Instala dependencias si faltan (Colab arranca sin Qiskit ni PySCF).
import importlib.util, os, subprocess, sys
REQS = {"qiskit": "qiskit==2.6.2", "qiskit_aer": "qiskit-aer==0.18.4", "pyscf": "pyscf==2.15.0",
        "qiskit_ibm_runtime": "qiskit-ibm-runtime==0.52.1"}   # este último, solo para FakeFez (§9)
missing = [pip for mod, pip in REQS.items() if importlib.util.find_spec(mod) is None]
if missing:
    subprocess.check_call([sys.executable, "-m", "pip", "install", "-q", *missing])
os.makedirs("h4vqe", exist_ok=True)
import time; T_INICIO = time.time()   # el resumen final dice cuánto tardó todo
print("python", sys.version.split()[0], "- instalado:", missing or "nada (ya estaba)")'''

# §1 ------------------------------------------------------------------------------------------
HAM_STRUCT = '''\
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from IPython.display import display, Math
from h4vqe.hamiltonian import molecule, qubit_hamiltonian
from h4vqe.historia import hamiltonian_structure, hamiltonian_latex
ham = {d: molecule(d) for d in (0.9, 2.4)}
Hd = {d: qubit_hamiltonian(m) for d, m in ham.items()}
for d, m in ham.items():
    e = m["hf"].mo_energy
    print(f"d = {d} Å: {len(Hd[d])} términos de Pauli en {Hd[d].num_qubits} qubits | E_nuc = {m['constante']:.4f} Ha | "
          f"energías de los orbitales: {np.round(e, 3)} Ha | brecha HOMO-LUMO = {e[2] - e[1]:.3f} Ha")
    display(Math(hamiltonian_latex(Hd[d], 8)))
with pd.option_context("display.float_format", lambda v: f"{v:.3f}", "display.max_colwidth", 90):
    display(hamiltonian_structure(Hd))'''

HAMILTONIAN = '''\
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

# §2 ------------------------------------------------------------------------------------------
MINIMO = '''\
# Lo mínimo: UCCSD ×1 y ×2. Error SIN ruido (con el estado: solo para elegir el ansatz), S² y costo.
import time
from h4vqe.historia import minimum_table
t0 = time.time()
minimo = minimum_table((1.8, 2.0, 2.4, 3.0))
with pd.option_context("display.float_format", lambda v: f"{v:.3g}"):
    display(minimo)
for a, g in minimo.groupby("ansatz"):
    w = g.loc[g["error sin ruido (mHa)"].idxmax()]
    print(f"{a}: error máx. {w['error sin ruido (mHa)']:.3g} mHa (a {w.d} Å), S² máx. {g['S²'].max():.2f}, "
          f"{w['CZ en línea']} CZ en línea -> P(ninguna CZ falla) = {w['P(ninguna CZ falla), línea']:.2f}")
print(f"({time.time() - t0:.0f} s)")'''

# §3 ------------------------------------------------------------------------------------------
ADAPT_CONV = '''\
# ADAPT al ir agregando operadores (en el orden en que los eligió), en las 18 distancias.
# Var(H) no usa FCI: es el criterio de corte. El error contra FCI es solo comprobación.
from h4vqe.adapt import ADAPT_FERMIONIC
from h4vqe.historia import adapt_convergence, exc_label
t0 = time.time()
conv_adapt = adapt_convergence(config.DISTANCES, ks=(8, 10, 12, 13, 15))
by_k, by_k_est = conv_adapt.groupby("operadores"), conv_adapt[conv_adapt.d >= config.STRETCHED_FROM].groupby("operadores")
res_adapt = pd.DataFrame({
    "CZ (todos con todos)": by_k.CZ.first(),
    "Var(H) máx., d ≥ 1,8 Å (Ha²)": by_k_est["Var(H) (Ha²)"].max(),
    "Var(H) máx., toda la curva (Ha²)": by_k["Var(H) (Ha²)"].max(),
    "error máx., toda la curva (mHa)": by_k["error (mHa)"].max(),
    "S² máx.": by_k["S²"].max()})
with pd.option_context("display.float_format", lambda v: f"{v:.3g}"):
    display(res_adapt)
print("lista final (15): " + " | ".join(f"{i + 1}: {exc_label(e)}" for i, e in enumerate(ADAPT_FERMIONIC)))
print(f"({time.time() - t0:.0f} s)")'''

CHECKS = '''\
import importlib
import h4vqe
for mod in ["config", "hamiltonian", "ansatz", "optimize", "estimator", "adapt", "routing", "analysis",
            "pipeline", "reduction", "noise", "mitigation", "hardware", "geometry", "scaling", "historia"]:
    importlib.reload(importlib.import_module("h4vqe." + mod))
from qiskit.quantum_info import Statevector
from h4vqe import config
from h4vqe.hamiltonian import molecule, qubit_hamiltonian
from h4vqe.ansatz import build_ansatz, cz_count
from h4vqe.pipeline import curve_ansatz_kwargs, line_cz_fswap

m = molecule(2.0); H = qubit_hamiltonian(m)
ans = build_ansatz(m, **curve_ansatz_kwargs())   # config.ANSATZ
e_hf = Statevector(ans.fast.psi0).expectation_value(H).real
print(f"términos de Pauli: {len(H)} | qubits: {H.num_qubits}")
print(f"E_HF (PySCF) = {m['e_hf']:.8f}   <HF|H|HF> = {e_hf:.8f}")
print(f"E_FCI = {m['e_fci']:.8f}   E_CCSD = {m['e_ccsd']:.8f}  (CCSD − FCI = {1e3*(m['e_ccsd']-m['e_fci']):+.2f} mHa)")
x = np.random.default_rng(0).normal(0, .2, ans.num_parameters)
fid = abs(np.vdot(Statevector(ans.circuit.assign_parameters(x)).data, ans.fast.state(x)))**2
print(f"fidelidad circuito Qiskit vs simulador numpy = {fid:.12f}")
print(f"ansatz: {config.ANSATZ} | parámetros: {ans.num_parameters} | CZ: {cz_count(ans.circuit)} (todos con todos), "
      f"{cz_count(ans.circuit, 'line')} (línea, transpiler), {line_cz_fswap(m, ans)} (línea, ruteo fSWAP)")'''

# §4 ------------------------------------------------------------------------------------------
SHOTS = '''\
# Shots para σ = 0,5 mHa con cada forma de agrupar (estado del ansatz de la curva; Neyman ideal).
# Solo para PLANIFICAR: las σ reportadas salen siempre de las muestras.
from h4vqe.historia import shots_table, GROUPINGS
t0 = time.time()
shots_tab = shots_table((0.9, 1.8, 2.4, 3.0))
tab = shots_tab.pivot_table(index="agrupamiento", columns="d", values="shots para σ = 0,5 mHa (M)")
tab = tab.loc[[name for _, _, name in GROUPINGS]].round(1)
tab.columns = [f"shots (M) a {d:g} Å" for d in tab.columns]
tab.insert(0, "grupos", shots_tab.groupby("agrupamiento").grupos.max())
display(tab)
s8 = shots_tab.pivot_table(index="agrupamiento", columns="d", values="σ con 8 M de shots (mHa)")
r = tab.iloc[0, 1:] / tab.iloc[1, 1:]
print(f"coloreo / sorted insertion (QWC): {r.min():.2f}–{r.max():.2f}× más shots | σ con 8 M de shots en la zona estirada: "
      f"{s8.loc[GROUPINGS[0][2], [1.8, 2.4, 3.0]].min():.2f}–{s8.loc[GROUPINGS[0][2], [1.8, 2.4, 3.0]].max():.2f} mHa (lo mínimo) -> "
      f"{s8.loc[GROUPINGS[1][2], [1.8, 2.4, 3.0]].min():.2f}–{s8.loc[GROUPINGS[1][2], [1.8, 2.4, 3.0]].max():.2f} mHa (la curva)")
print(f"({time.time() - t0:.0f} s)")'''

# §5 ------------------------------------------------------------------------------------------
RUN = '''\
from h4vqe.pipeline import run_curve
QUICK = False   # True: 4 distancias y menos shots, para probar rápido
dists = [0.7, 0.9, 1.8, 3.0] if QUICK else config.DISTANCES
shots = 400_000 if QUICK else config.SHOTS_PER_POINT
t0 = time.time()
rows, ctxs = run_curve(dists, shots=shots)
print(f"total: {time.time()-t0:.0f} s")'''

TABLE = '''\
from h4vqe.analysis import results_table
df = results_table(rows)
cols = ["d","E_vqe","sigma_mHa","shots","cz_all","cz_line","cz_line_fswap","E_fci","E_ccsd","err_vqe_mHa","err_ccsd_mHa","beats_ccsd","chem_acc"]
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

# §6 ------------------------------------------------------------------------------------------
LIMIT = '''\
# Límite físico sin FCI: 4 átomos de H separados (con 1 electrón, HF es exacto en la base).
from h4vqe.historia import isolated_atoms_energy
E_atoms = isolated_atoms_energy()
lejos = df[df.d == df.d.max()].iloc[0]
print(f"4 átomos de H separados: 4·E_H = {E_atoms:.6f} Ha")
for nombre, E in [("VQE (medido)", lejos.E_vqe), ("CCSD", lejos.E_ccsd)]:
    print(f"  a {lejos.d:.1f} Å, {nombre:12s}: E − 4·E_H = {1e3 * (E - E_atoms):+7.2f} mHa")'''

RELIABILITY = '''\
from h4vqe.pipeline import run_reliability
# Tres puntos de la zona estirada; en cada uno, HF, un ansatz menos profundo (ADAPT 12 ops),
# el de la curva (ADAPT 15 ops) y uno más profundo (UCCSD ×2). Ver pipeline.depth_references.
rel_d = [1.8, 3.0] if QUICK else [1.8, 2.4, 3.0]
t0 = time.time()
rel = pd.DataFrame(run_reliability(ctxs, rel_d, shots=shots,
                                   var_shots=500_000 if QUICK else 2_000_000,
                                   echo_shots=200_000 if QUICK else 1_000_000,
                                   sym_shots=100_000 if QUICK else 200_000))
print(f"total: {time.time()-t0:.0f} s")
# Var(grupos) puede salir negativa por ruido (σ ≈ 2e-3 Ha²): se muestra recortada en 0, con su σ.
# Solo separa HF del resto; la que resuelve es la del eco.
rel["var_grupos (≥0)"] = rel["var"].clip(lower=0)
# echo_cz: CZ del circuito de eco (U, e^{−iHt}, U†). Es un diagnóstico de simulador, no de hardware.
show = rel[["d","rol","estado","cz","E","sigma","var_grupos (≥0)","var_sigma","var_echo","var_echo_up95","echo_events",
            "echo_cz","N_alpha","N_alpha_sigma","N_beta","N_beta_sigma","S2","S2_sigma","var_exacta"]]
with pd.option_context("display.float_format", lambda v: f"{v:.3g}"):
    display(show)'''

DEPTH = '''\
# Convergencia con la profundidad (sin FCI). Si el ansatz de la curva convergió:
#   - el menos profundo da una energía MEDIDA más alta y una varianza (eco) mayor;
#   - el más profundo (muchas más CZ) no baja la energía medida más allá de σ.
E = rel.pivot(index="d", columns="rol", values="E")
S = rel.pivot(index="d", columns="rol", values="sigma")
by = lambda rol, col: rel[rel.rol == rol].set_index("d")[col]
conv = pd.DataFrame({
    "E(menos) − E(curva) [mHa]": 1e3*(E["menos profundo"] - E["curva"]),
    "σ₁ [mHa]": 1e3*np.hypot(S["menos profundo"], S["curva"]),
})
if "más profundo" in E:
    conv["E(curva) − E(más) [mHa]"] = 1e3*(E["curva"] - E["más profundo"])
    conv["σ₂ [mHa]"] = 1e3*np.hypot(S["curva"], S["más profundo"])
conv["Var eco menos [Ha²]"] = by("menos profundo", "var_echo")
conv["Var eco curva ≤ [Ha²] (95 %)"] = by("curva", "var_echo_up95")
with pd.option_context("display.float_format", lambda v: f"{v:.3g}"):
    display(conv)
# Recién ahora, como comprobación: el error real contra FCI.
chk = df.set_index("d").loc[rel_d, ["err_vqe_mHa", "sigma_mHa", "err_ccsd_mHa"]]
chk'''

# §7 ------------------------------------------------------------------------------------------
NOISE = '''\
# Los tres candidatos con el modelo de ruido (qubits en línea), a 2,4 Å: valor esperado EXACTO con
# matriz densidad (sin shots), o sea el sesgo que deja el ruido y que más shots no sacan.
from h4vqe.historia import noisy_ansatz_bias
t0 = time.time()
ruido_tab, ruido_ctx = noisy_ansatz_bias(2.4)
with pd.option_context("display.float_format", lambda v: f"{v:.3g}"):
    display(ruido_tab)
print(f"CCSD − FCI a 2,4 Å = {1e3 * (ruido_ctx['m']['e_ccsd'] - ruido_ctx['m']['e_fci']):+.1f} mHa  ({time.time() - t0:.0f} s)")'''

REDUCTION = '''\
from h4vqe.reduction import reduction_table, levels_table, compare_ansatze
red = pd.DataFrame(reduction_table(config.DISTANCES))
print(f"{len(red)} distancias | máx |E_0(Ĥ₆, Ĥ₅) − E_FCI| = {red.err_fci.abs().max():.1e} Ha | "
      f"máx |<HF|Ĥ₆, Ĥ₅|HF> − E_HF| = {red.err_hf.abs().max():.1e} Ha")
print(f"8 → 6 → 5 qubits | términos de Pauli: {red.paulis_8q.max()} → {red.paulis_6q.max()} → {red.paulis_5q.max()} | "
      f"grupos QWC (máx. en la curva): {red.groups_8q.max()} → {red.groups_6q.max()} → {red.groups_5q.max()}")
display(pd.DataFrame(levels_table(2.4)))   # zona estirada
pd.DataFrame(compare_ansatze(2.0))'''

# §8 ------------------------------------------------------------------------------------------
MITIGATION = '''\
# Sesgo EXACTO de cada técnica (sin shots: probabilidades con matriz densidad), a 2,4 Å, con el
# ansatz compacto (5 qubits, 16 CZ), el modelo de noise.py y 2 % de error de lectura.
# Molécula, Ĥ₅ y θ* son los de §7.1 (MISMO proceso: ver "qué no funcionó").
from qiskit.quantum_info import Statevector
from h4vqe import mitigation as mt
from h4vqe.estimator import measurement_groups, split_constant
from h4vqe.noise import to_line
from h4vqe.reduction import reduced_number_operators
t0 = time.time()
d, m, H5 = ruido_ctx["d"], ruido_ctx["m"], ruido_ctx["H5"]
P_RO, SIG = 0.02, 0.5e-3
ops = reduced_number_operators(m)
bound = ruido_ctx["compact"].circuit.assign_parameters(ruido_ctx["x_compact"]); isa = to_line(bound)
E0 = Statevector(bound).expectation_value(H5).real
groups = measurement_groups(H5); sym = [(ops["N_alpha"], 2), (ops["N_beta"], 2)]
ro = mt.readout_matrices(H5.num_qubits, P_RO)
P = {l: mt.exact_probabilities(isa, groups, p_ro=P_RO, lam=l) for l in (1, 3, 5)}
tabs = {"sin mitigar": mt.mitigation_tables(groups), "RO": mt.mitigation_tables(groups, readout=ro),
        "PS + RO": mt.mitigation_tables(groups, symmetries=sym, readout=ro)}
est = {k: {l: mt.exact_estimate(P[l], H5, U, W) for l in (1, 3, 5)} for k, (U, W, _) in tabs.items()}
F, emix = mt.echo_fidelity_exact(isa), split_constant(H5)[0]
def gam(lams):
    lams = np.asarray(lams, float)
    return np.array([np.prod([b / (b - a) for j, b in enumerate(lams) if j != i]) for i, a in enumerate(lams)])
mitig_rows = []
def add(name, base, lams, rescale=False):
    E = np.array([est[base][l]["E"] for l in lams]); K = np.array([est[base][l]["K"] for l in lams])
    if rescale:   # eco: E' = (E − (1 − F^λ) E_mix) / F^λ
        fl = F ** np.asarray(lams, float); E, K = (E - (1 - fl) * emix) / fl, K / fl ** 2
    g = gam(lams)
    mitig_rows.append({"técnica": name, "sesgo (mHa)": 1e3 * (g @ E - E0),
                       "shots para σ = 0,5 mHa (M)": float(np.sum(np.abs(g) * np.sqrt(K)) ** 2) / SIG ** 2 / 1e6})
for k in tabs: add(k, k, (1,))
add("ZNE lineal (1,3) + PS + RO", "PS + RO", (1, 3))
add("eco + ZNE lineal (1,3) + RO", "RO", (1, 3), rescale=True)
add("ZNE Richardson (1,3,5) + PS + RO", "PS + RO", (1, 3, 5))
mitig = pd.DataFrame(mitig_rows).round(3)
print(f"d = {d} Å | {int(isa.count_ops()['cz'])} CZ | F(eco) = {F:.4f} | E0 − FCI = {1e3*(E0 - m['e_fci']):.1e} mHa | "
      f"CCSD − FCI = {1e3*(m['e_ccsd']-m['e_fci']):.1f} mHa | {time.time()-t0:.0f} s")
mitig'''

# §9 ------------------------------------------------------------------------------------------
HARDWARE = '''\
from h4vqe.hardware import run_fez
fz = run_fez(2.4)
print(f"línea de Fez: qubits {fz['line']} | errores esperados {fz['expected_errors']:.3f} "
      f"(mediana de las {fz['n_lines']} líneas posibles: {fz['median_expected_errors']:.3f})")
print(f"CZ: compacto {fz['cz']} (sin SWAP) | UCCSD ×2 ubicado en Fez: {fz['cz_uccsd_fez']} | "
      f"error del compacto sin ruido: {fz['err_exact_mHa']:.1e} mHa")
print("lectura por qubit, P(1|0) / P(0|1): " + ", ".join(f"{a:.2%} / {b:.2%}" for a, b in fz["readout"]))
pd.DataFrame(fz["rows"]).astype({"shots": "Int64"})'''

# §10 -----------------------------------------------------------------------------------------
GEOMETRY = '''\
from h4vqe.geometry import scan_geometry
from h4vqe.hamiltonian import ZIGZAG_ANGLE
GEO_D = [0.7, 1.2, 1.8, 2.4, 3.0]
GEO_SHOTS = 2_000_000
geo = pd.DataFrame(scan_geometry("square", GEO_D, GEO_SHOTS) + scan_geometry("zigzag", GEO_D, GEO_SHOTS))
geo["beats_ccsd"] = geo.err_vqe_mHa.abs() + 1e3 * geo.sigma < geo.err_ccsd_mHa.abs()
geo["err_sin_ruido_mHa"] = 1e3 * (geo.E_opt_exact - geo.E_fci)   # θ* con el estado exacto (sin shots)
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
geo[["shape", "d", "E_vqe", "sigma", "err_vqe_mHa", "err_sin_ruido_mHa", "err_ccsd_mHa", "ccsd_converged", "w_hf",
     "beats_ccsd", "cz_all", "cz_line"]]'''

# §12 -----------------------------------------------------------------------------------------
SCALING = '''\
from h4vqe.scaling import D_EQ, D_STRETCH, chain_resources, extrapolate, plot_scaling, uccsd_cz_transpiled
t0 = time.time()
esc = pd.read_csv("h4vqe/escalado_recursos.csv")    # precalculado (scripts/scaling.py): N = 4…14
# Reproducibilidad: recalculamos en vivo N = 4, 6 y 8 y comparamos con lo guardado.
vivo = pd.DataFrame([chain_resources(n, dd, verbose=False) for n in (4, 6, 8) for dd in (D_EQ, D_STRETCH)])
jj = vivo.merge(esc, on=["N", "d"], suffixes=("", "_guardado"))
dif = max(float(np.max(np.abs(jj[c] / jj[c + "_guardado"] - 1))) for c in ["terminos", "grupos_qwc", "params", "cz_cota", "shots_fci"])
print(f"recalculado en vivo N = 4, 6, 8: máxima diferencia relativa con lo guardado = {dif:.1e} "
      f"| |ΔCCSD| máx. = {np.max(np.abs(jj.ccsd_menos_fci_mHa - jj.ccsd_menos_fci_mHa_guardado)):.1e} mHa")
print(f"CZ de UCCSD ×1 para H₄ transpilado acá: {uccsd_cz_transpiled(4)} (guardado: "
      f"{int(esc.cz_transpilado[esc.N == 4].iloc[0])}; §2: {int(minimo.loc[minimo.ansatz == 'UCCSD ×1', 'CZ todos con todos'].iloc[0])})")
ext, lineas = extrapolate(esc)
print("\\n".join(lineas))
fig = plot_scaling(esc); plt.show()
filas = []
for N in (4, 10, 14):
    r = esc[(esc.N == N) & np.isclose(esc.d, D_STRETCH)].iloc[0]
    filas.append({"sistema": f"H{N}", "qubits": 2 * N, "términos": r.terminos, "grupos QWC": r.grupos_qwc,
                  "parámetros UCCSD": r.params, "CZ UCCSD ×1": r.cz_transpilado, "log10 P(ninguna CZ falla)": r.log10_p_sin_error,
                  "shots por energía (σ = 1,6 mHa)": r.shots_fci, "|CCSD − FCI| (mHa)": abs(r.ccsd_menos_fci_mHa)})
for _, r in ext[np.isclose(ext.d, D_STRETCH)].iterrows():
    filas.append({"sistema": f"FeMoco, {r.modelo} (extrapolado)", "qubits": r.qubits, "términos": r.terminos,
                  "grupos QWC": r.grupos_qwc, "parámetros UCCSD": r.params, "CZ UCCSD ×1": r.cz_transpilado,
                  "log10 P(ninguna CZ falla)": r.log10_p_sin_error, "shots por energía (σ = 1,6 mHa)": r.shots_fci,
                  "|CCSD − FCI| (mHa)": np.nan})
escala = pd.DataFrame(filas)
with pd.option_context("display.float_format", lambda v: f"{v:.3g}"):
    display(escala)
print(f"(a {D_STRETCH} Å; {time.time() - t0:.0f} s)")'''

# §13 -----------------------------------------------------------------------------------------
SUMMARY = '''\
# Resumen: todo sale de las celdas de arriba (acá no se recalcula nada).
from h4vqe.historia import plot_evolution
st = df[df.d >= config.STRETCHED_FROM]
u2 = minimo[minimo.ansatz == "UCCSD ×2"]
cuad = cmp.set_index("modelo").loc["cuadrático 5"]
sh = shots_tab.set_index(["d", "agrupamiento"])
rz = ruido_tab.set_index("ansatz")
mi = mitig.set_index("técnica")
fzr = pd.DataFrame(fz["rows"]).set_index("método")
cur = rel[rel.rol == "curva"]
gs = geo[(geo.d >= config.STRETCHED_FROM) & geo.ccsd_converged]
fem = ext[np.isclose(ext.d, D_STRETCH) & (ext.modelo == "Reiher 2017")].iloc[0]
rng = lambda s, f="{:.2f}": f"{f.format(s.min())}–{f.format(s.max())}"
ccsd24 = 1e3 * abs(ruido_ctx["m"]["e_ccsd"] - ruido_ctx["m"]["e_fci"])
fez_z = fzr.loc["VQE en Fez, ZNE (1,3,5) + lectura"]

print("1) NUESTROS MEJORES RESULTADOS")
mejores = pd.DataFrame([
    ("Curva central, d ≥ 1,8 Å (8 M shots por punto)",
     f"error ≤ {st.err_vqe_mHa.max():.2f} mHa con σ = {rng(st.sigma_mHa)} mHa",
     f"{int(st.beats_ccsd.sum())}/{len(st)} le gana a CCSD ({rng(st.err_ccsd_mHa, '{:.0f}')} mHa); {int(st.chem_acc.sum())}/{len(st)} con precisión química"),
    ("Curva central, las 18 distancias", f"error ≤ {df.err_vqe_mHa.max():.2f} mHa",
     f"{int(df.chem_acc.sum())}/{len(df)} con precisión química"),
    ("R_e", f"{fit['R_e']:.4f} ± {fit['R_e_sigma']:.4f} Å", f"FCI: {ref['R_e']:.4f} Å ({(fit['R_e'] - ref['R_e']) / fit['R_e_sigma']:+.1f}σ)"),
    ("D_e", f"{1e3 * de['D_e']:.2f} ± {1e3 * de['D_e_sigma']:.2f} mHa", f"FCI: {1e3 * De_fci:.2f} mHa ({(de['D_e'] - De_fci) / de['D_e_sigma']:+.1f}σ)"),
    ("Confiabilidad medida (sin FCI), 1,8–3,0 Å",
     f"Var(eco) ≤ {cur.var_echo_up95.max():.1e} Ha²; |S²| ≤ {cur.S2.abs().max():.3f}; N_α = N_β = {cur.N_alpha.mean():.0f} en todos los shots",
     f"HF: Var = {rng(rel[rel.rol == 'HF'].var_echo)} Ha²"),
    ("Con ruido (0,1 %/CZ, 2 % de lectura), 2,4 Å: compacto 5q + eco + ZNE + RO",
     f"sesgo exacto {mi.loc['eco + ZNE lineal (1,3) + RO', 'sesgo (mHa)']:+.3f} mHa ({mi.loc['eco + ZNE lineal (1,3) + RO', 'shots para σ = 0,5 mHa (M)']:.0f} M shots para σ = 0,5 mHa)",
     "medido con shots (§8): 6/6 con precisión química en 1,8, 2,4 y 3,0 Å"),
    ("FakeFez (ruido de IBM Fez), 2,4 Å: lectura + ZNE",
     f"{fez_z['E − FCI (mHa)']:+.2f} ± {fez_z['sigma_mHa']:.2f} mHa",
     f"sin mitigar: {fzr.loc['VQE en Fez, sin mitigar', 'E − FCI (mHa)']:+.1f} mHa; CCSD: {fzr.loc['CCSD', 'E − FCI (mHa)']:+.1f} mHa"),
    ("Otras geometrías, d ≥ 1,8 Å donde CCSD converge (UCCSD ×2, 2 M shots)",
     f"|VQE − FCI| ≤ {gs.err_vqe_mHa.abs().max():.2f} mHa (σ ≈ {1e3 * gs.sigma.mean():.1f} mHa)",
     f"le gana a CCSD en {int(gs.beats_ccsd.sum())}/{len(gs)} puntos ({rng(gs.err_ccsd_mHa.abs(), '{:.0f}')} mHa)"),
    ("Escalado a FeMoco (54 orbitales, 108 qubits)",
     f"{fem.params:.2e} parámetros, {fem.cz_transpilado:.1e} CZ, {fem.shots_fci:.1e} shots por energía",
     f"P(ninguna CZ falla) = 10^{fem.log10_p_sin_error:.0f}"),
], columns=["qué", "resultado", "comparación"])
with pd.option_context("display.max_colwidth", 140):
    display(mejores)

print("2) LO MÍNIMO CONTRA LO FINAL")
uccsd_err = u2[u2.d >= config.STRETCHED_FROM]["error sin ruido (mHa)"].max()
adapt_err = 1e3 * (st.E_opt_exact - st.E_fci).max()
lvf = pd.DataFrame([
    ("ansatz", "UCCSD ×2 (Trotter)", "ADAPT-VQE, 15 operadores", "elegidos por gradiente"),
    ("parámetros", int(u2["parámetros"].iloc[0]), int(df.n_params.iloc[0]), f"{u2['parámetros'].iloc[0] / df.n_params.iloc[0]:.1f}× menos"),
    ("CZ, todos con todos", int(u2["CZ todos con todos"].iloc[0]), int(df.cz_all.iloc[0]), f"{u2['CZ todos con todos'].iloc[0] / df.cz_all.iloc[0]:.1f}× menos"),
    ("CZ, qubits en línea", int(u2["CZ en línea"].iloc[0]), int(df.cz_line_fswap.iloc[0]), f"{u2['CZ en línea'].iloc[0] / df.cz_line_fswap.iloc[0]:.1f}× menos (ruteo fSWAP)"),
    ("grupos de medición (QWC)", int(sh.loc[(2.4, GROUPINGS[0][2]), "grupos"]), int(sh.loc[(2.4, GROUPINGS[1][2]), "grupos"]), "sorted insertion"),
    ("shots para σ = 0,5 mHa a 2,4 Å", f"{sh.loc[(2.4, GROUPINGS[0][2]), 'shots para σ = 0,5 mHa (M)']:.1f} M", f"{sh.loc[(2.4, GROUPINGS[1][2]), 'shots para σ = 0,5 mHa (M)']:.1f} M",
     f"{sh.loc[(2.4, GROUPINGS[0][2]), 'shots para σ = 0,5 mHa (M)'] / sh.loc[(2.4, GROUPINGS[1][2]), 'shots para σ = 0,5 mHa (M)']:.1f}× menos"),
    ("σ con 8 M de shots, d ≥ 1,8 Å", "0,65–0,72 mHa (primera corrida)", f"{rng(st.sigma_mHa)} mHa (medida)", "misma cantidad de shots"),
    ("error sin ruido máx., d ≥ 1,8 Å", f"{uccsd_err:.1e} mHa", f"{adapt_err:.2f} mHa", "ambos muy por debajo de σ"),
    ("error medido máx., d ≥ 1,8 Å", "0,74 mHa (primera corrida)", f"{st.err_vqe_mHa.max():.2f} mHa", "7/7 contra CCSD en los dos"),
    ("R_e − R_e(FCI)", f"{1e3 * (cuad['R_e (Å)'] - ref['R_e']):+.1f} mÅ (cuadrático, mismos datos)", f"{1e3 * (fit['R_e'] - ref['R_e']):+.1f} ± {1e3 * fit['R_e_sigma']:.1f} mÅ (cuártico)", "ajuste elegido por su sesgo"),
    ("confiabilidad", "con el estado del simulador", "medida con shots (eco, S², N, profundidad)", "sin FCI y sin el estado"),
    ("error con ruido a 2,4 Å, sin mitigar", f"{rz.loc['UCCSD ×2 (lo mínimo)', 'error con ruido, sin mitigar (mHa)']:+.0f} mHa (peor que HF: {rz.loc['HF (sin circuito)', 'error con ruido, sin mitigar (mHa)']:+.0f})",
     f"{rz.loc['ADAPT 15 + ruteo fSWAP', 'error con ruido, sin mitigar (mHa)']:+.0f} mHa (ADAPT); {rz.loc['compacto 5 qubits (L = 8)', 'error con ruido, sin mitigar (mHa)']:+.1f} mHa (compacto 5q)", "§7"),
    ("error con ruido a 2,4 Å, mitigado", "—", f"{mi.loc['eco + ZNE lineal (1,3) + RO', 'sesgo (mHa)']:+.3f} mHa (exacto); −1,0 ± 0,5 y −0,1 ± 0,5 mHa (medidos)", "§8"),
], columns=["", "lo mínimo", "final", "cómo"])
with pd.option_context("display.max_colwidth", 120):
    display(lvf)

print("3) LA EVOLUCIÓN CON RUIDO (2,4 Å, sesgo exacto con matriz densidad)")
cp = rz.loc["compacto 5 qubits (L = 8)"]
evo = pd.DataFrame([
    {"etapa": "UCCSD ×2\\n(lo mínimo)", **rz.loc["UCCSD ×2 (lo mínimo)", ["CZ (línea)", "error sin ruido (mHa)", "error con ruido, sin mitigar (mHa)"]].to_dict()},
    {"etapa": "ADAPT 15 + fSWAP\\n(la curva)", **rz.loc["ADAPT 15 + ruteo fSWAP", ["CZ (línea)", "error sin ruido (mHa)", "error con ruido, sin mitigar (mHa)"]].to_dict()},
    {"etapa": "compacto, 5 qubits", **cp[["CZ (línea)", "error sin ruido (mHa)", "error con ruido, sin mitigar (mHa)"]].to_dict()},
    {"etapa": "compacto + eco + ZNE\\n+ lectura", "CZ (línea)": cp["CZ (línea)"], "error sin ruido (mHa)": cp["error sin ruido (mHa)"],
     "error con ruido, sin mitigar (mHa)": mi.loc["eco + ZNE lineal (1,3) + RO", "sesgo (mHa)"]},
]).rename(columns={"CZ (línea)": "CZ en línea", "error con ruido, sin mitigar (mHa)": "error con ruido (mHa)"})
evo["|CCSD − FCI| (mHa)"] = ccsd24
evo["CZ en línea"] = evo["CZ en línea"].astype(int)
with pd.option_context("display.float_format", lambda v: f"{v:.3g}"):
    display(evo.assign(etapa=evo.etapa.str.replace("\\n", " ")))
fig = plot_evolution(evo); plt.show()
print(f"tiempo total de ejecución del notebook: {(time.time() - T_INICIO) / 60:.0f} min")'''


def md(name: str):
    return nbf.v4.new_markdown_cell((MD / name).read_text(encoding="utf-8"))


def code(src: str):
    return nbf.v4.new_code_cell(src)


# Metadatos que colapsan el código de una celda: Colab (cellView = form, con el # @title de la
# primera línea) y Jupyter/VS Code (jupyter.source_hidden). La celda se ejecuta igual.
HIDDEN = {"cellView": "form", "jupyter": {"source_hidden": True}}

# En el notebook va solo el nombre del equipo, al principio (00_intro.md). Las líneas de "dueño"
# de los módulos quedan en el repo, pero no se copian al notebook.
_NAME_FIXES = [(" (Octavio)", ""), ("reducción de Octavio:", "reducción a 5 qubits:"),
               ("(lo que recomienda Octavio)", "(la combinación recomendada)"), ("(ADAPT, Gonzalo)", "(ADAPT)")]
_NAMES = re.compile(r"Octavio|Gonzalo|Santino|Joaquín")


def clean_module(src: str) -> str:
    lines = src.split("\n")
    if "Dueño:" in src:
        i = next(k for k, ln in enumerate(lines) if ln.startswith("Dueño:"))
        j = i + 1      # líneas de continuación: no vacías y sin cerrar el docstring
        while j < len(lines) and lines[j].strip() and not lines[j].lstrip().startswith('"""'):
            j += 1
        if j < len(lines) and not lines[j].strip() and i > 0 and not lines[i - 1].strip():
            j += 1     # no dejar dos líneas en blanco seguidas
        lines = lines[:i] + lines[j:]
    src = "\n".join(lines)
    for old, new in _NAME_FIXES:
        src = src.replace(old, new)
    assert not _NAMES.search(src), _NAMES.search(src)
    return src


def hidden(src: str):
    c = code(src)
    c.metadata.update(HIDDEN)
    return c


def _modules():
    """Una sola celda, colapsada, que escribe el paquete h4vqe/ (antes eran 17 celdas %%writefile)."""
    files = {f"{m}.py": clean_module((PKG / f"{m}.py").read_text(encoding="utf-8")) for m in MODULES}
    files["escalado_recursos.csv"] = SCALING_CSV.read_text(encoding="utf-8")
    body = []
    for name, src in files.items():
        esc = src.replace("\\", "\\\\").replace("\'\'\'", "\\\'\\\'\\\'")
        body.append(f"    {name!r}: \'\'\'{esc}\'\'\',")
    return hidden("# @title Paquete `h4vqe/`: escribe los módulos (código oculto, doble clic para verlo)\n"
                  "import os\n"
                  "os.makedirs(\"h4vqe\", exist_ok=True)\n"
                  "ARCHIVOS = {\n" + "\n".join(body) + "\n}\n"
                  "for nombre, fuente in ARCHIVOS.items():\n"
                  "    with open(os.path.join(\"h4vqe\", nombre), \"w\", encoding=\"utf-8\") as f:\n"
                  "        f.write(fuente)\n"
                  "print(f\"h4vqe/: {len(ARCHIVOS)} archivos escritos\")")


def build(quick: bool = False) -> nbf.NotebookNode:
    run = RUN.replace("QUICK = False", "QUICK = True ") if quick else RUN
    cells = [
        md("00_intro.md"),
        md("01_setup.md"), hidden("# @title Instala las dependencias (código oculto)\n" + INSTALL), _modules(),
        md("02_problema.md"), code(HAM_STRUCT),
        md("02_lectura.md"), code(HAMILTONIAN),
        md("02_ccsd.md"),
        md("03_minimo.md"), code(MINIMO),
        md("03_hipotesis.md"),
        md("04_adapt.md"), code(ADAPT_CONV),
        md("04_ruteo.md"), code(CHECKS),
        md("05_medicion.md"), code(SHOTS),
        md("05_calibracion.md"),
        md("06_curva.md"), code(run), code(TABLE), code(PLOT),
        md("06_ajuste.md"), code(RE_DE),
        md("06_resultado.md"),
        md("07_confiabilidad.md"), code(LIMIT), code(RELIABILITY), code(DEPTH),
        md("07_veredicto.md"),
        md("08_ruido.md"), code(NOISE),
        md("08_reduccion.md"), code(REDUCTION),
        md("08_reduccion_resultado.md"),
        md("09_mitigacion.md"), code(MITIGATION),
        md("09_mitigacion_shots.md"),
        md("10_fakefez.md"), code(HARDWARE),
        md("11_geometria.md"), code(GEOMETRY),
        md("11_geometria_resultado.md"),
        md("12_ansatze.md"),
        md("13_escalado.md"), code(SCALING),
        md("13_escalado_resultado.md"),
        md("14_resumen.md"), code(SUMMARY),
        md("14_conclusiones.md"),
    ]
    nb = nbf.v4.new_notebook()
    nb.cells = cells
    nb.metadata["kernelspec"] = {"name": "python3", "display_name": "Python 3", "language": "python"}
    nb.metadata["colab"] = {"provenance": [], "toc_visible": True}
    todo = [c.source[:80] for c in cells if c.cell_type == "markdown" and "TODO" in c.source]
    assert not todo, f"quedó un TODO en el texto: {todo}"
    return nb


def execute(nb: nbf.NotebookNode, kernel: str) -> nbf.NotebookNode:
    """Ejecuta en una carpeta temporal vacía, como un Colab nuevo."""
    from nbclient import NotebookClient

    workdir = pathlib.Path(tempfile.mkdtemp(prefix="h4vqe_nb_"))
    try:
        t0 = time.time()
        NotebookClient(nb, timeout=None, kernel_name=kernel,
                       resources={"metadata": {"path": str(workdir)}}).execute()
        print(f"  ejecutado OK en {time.time() - t0:.0f} s")
    finally:
        shutil.rmtree(workdir, ignore_errors=True)
    return nb


def update_text() -> None:
    """Reemplaza solo las celdas de texto del notebook ya ejecutado (las salidas quedan)."""
    old, new = nbf.read(OUT, as_version=4), build()
    if [c.cell_type for c in old.cells] != [c.cell_type for c in new.cells]:
        raise ValueError("cambió la estructura de celdas, hay que volver a ejecutar con --run")
    for c_old, c_new in zip(old.cells, new.cells):
        if c_new.cell_type == "markdown":
            c_old.source = c_new.source
        elif c_old.source != c_new.source:
            raise ValueError("cambió una celda de código, hay que volver a ejecutar con --run")
    nbf.write(old, OUT)
    print("texto actualizado en", OUT)



MD_CORTO = ROOT / "notebook" / "md_corto"
OUT_CORTO = ROOT / "notebook" / "H4_VQE_corto.ipynb"


def build_corto() -> None:
    """Copia del notebook YA EJECUTADO con el texto corto de notebook/md_corto: las mismas celdas de
    código y las mismas salidas; solo cambian las celdas de texto. No hace falta volver a correr."""
    global MD
    md_orig, MD = MD, MD_CORTO
    try:
        new = build()
    finally:
        MD = md_orig
    nb = nbf.read(OUT, as_version=4)
    if [c.cell_type for c in nb.cells] != [c.cell_type for c in new.cells]:
        raise ValueError("la estructura de celdas no coincide: regenerá y ejecutá H4_VQE.ipynb primero")
    for c_old, c_new in zip(nb.cells, new.cells):
        if c_new.cell_type == "markdown":
            c_old.source = c_new.source
        elif c_old.source != c_new.source:
            raise ValueError("cambió una celda de código: ejecutá H4_VQE.ipynb con --run primero")
    nbf.write(nb, OUT_CORTO)
    print("escrito", OUT_CORTO)

if __name__ == "__main__":
    args = sys.argv[1:]
    if "--corto" in args:
        build_corto()
        sys.exit()
    kernel = args[args.index("--kernel") + 1] if "--kernel" in args else "h4vqe-local"
    if "--text-only" in args:
        update_text()
    else:
        nb = build(quick="--quick" in args)
        nbf.write(nb, OUT)
        print("escrito", OUT)
        if "--run" in args:
            nbf.write(execute(nb, kernel), OUT)
