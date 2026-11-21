"""Genera los dos notebooks de entrega a partir de h4vqe/*.py y notebook/md/*.md.

    notebook/H4_VQE.ipynb         lo que pide la consigna (curva, R_e, D_e, confiabilidad)
    notebook/H4_VQE_extras.ipynb  extras: menos qubits, mitigación, FakeFez, otra geometría

    python notebook/build_notebook.py                  # genera los dos
    python notebook/build_notebook.py --run            # genera y ejecuta los dos en carpetas vacías
    python notebook/build_notebook.py --run core       # solo uno: core o extras
    python notebook/build_notebook.py --quick --run    # QUICK = True en la curva (antes de cada PR)
    python notebook/build_notebook.py --kernel NOMBRE  # kernel para ejecutar (default h4vqe-local)
    python notebook/build_notebook.py --text-only      # actualiza solo el texto de los ya ejecutados

No editar los .ipynb a mano: se pisan. Editar los módulos o los .md y regenerar.
La ejecución en una carpeta temporal vacía (idea de Santino, build_minimo.py) reproduce un Colab
nuevo: %%writefile crea h4vqe/ ahí y no se usa el paquete del repo.
"""
from __future__ import annotations

import pathlib
import shutil
import sys
import tempfile
import time

import nbformat as nbf

ROOT = pathlib.Path(__file__).resolve().parent.parent
MD = ROOT / "notebook" / "md"
PKG = ROOT / "h4vqe"
OUT = {"core": ROOT / "notebook" / "H4_VQE.ipynb",
       "extras": ROOT / "notebook" / "H4_VQE_extras.ipynb"}

CORE_MODULES = ["__init__", "config", "hamiltonian", "ansatz", "optimize", "estimator",
                "adapt", "routing", "analysis", "pipeline"]
EXTRA_MODULES = CORE_MODULES + ["reduction", "noise", "mitigation", "hardware", "geometry"]

INSTALL = '''\
# Instala dependencias si faltan (Colab arranca sin Qiskit ni PySCF).
import importlib.util, os, subprocess, sys
REQS = {"qiskit": "qiskit==2.6.2", "qiskit_aer": "qiskit-aer==0.18.4", "pyscf": "pyscf==2.15.0"}
missing = [pip for mod, pip in REQS.items() if importlib.util.find_spec(mod) is None]
if missing:
    subprocess.check_call([sys.executable, "-m", "pip", "install", "-q", *missing])
os.makedirs("h4vqe", exist_ok=True)
print("python", sys.version.split()[0], "- instalado:", missing or "nada (ya estaba)")'''

INSTALL_EXTRAS = INSTALL.replace(
    '"pyscf": "pyscf==2.15.0"}',
    '"pyscf": "pyscf==2.15.0",\n        "qiskit_ibm_runtime": "qiskit-ibm-runtime==0.52.1"}   # solo para FakeFez')
assert INSTALL_EXTRAS != INSTALL

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

REDUCTION = '''\
import pandas as pd
from h4vqe import config
from IPython.display import display
from h4vqe.reduction import reduction_table, levels_table, compare_ansatze
red = pd.DataFrame(reduction_table(config.DISTANCES))
print(f"{len(red)} distancias | máx |E_0(Ĥ₆, Ĥ₅) − E_FCI| = {red.err_fci.abs().max():.1e} Ha | "
      f"máx |<HF|Ĥ₆, Ĥ₅|HF> − E_HF| = {red.err_hf.abs().max():.1e} Ha")
print(f"8 → 6 → 5 qubits | términos de Pauli: {red.paulis_8q.max()} → {red.paulis_6q.max()} → {red.paulis_5q.max()} | "
      f"grupos QWC (máx. en la curva): {red.groups_8q.max()} → {red.groups_6q.max()} → {red.groups_5q.max()}")
display(pd.DataFrame(levels_table(2.4)))   # zona estirada
pd.DataFrame(compare_ansatze(2.0))'''

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

MITIGATION = '''\
# Sesgo EXACTO de cada técnica (sin shots: probabilidades con matriz densidad), a 2,4 Å, con el
# ansatz compacto (5 qubits, 16 CZ), el modelo de noise.py y 2 % de error de lectura.
# Sirve para elegir la técnica antes de gastar shots; las energías con shots están en la tabla de arriba
# (scripts/mitigation_shots.py, ~20 min por punto). θ* y Ĥ₅ salen del MISMO proceso (ver "qué no funcionó").
import time
import numpy as np, pandas as pd
from qiskit.quantum_info import Statevector
from h4vqe import mitigation as mt
from h4vqe.estimator import measurement_groups, split_constant
from h4vqe.hamiltonian import molecule, qubit_hamiltonian
from h4vqe.noise import to_line
from h4vqe.reduction import build_compact_ansatz, optimize_compact, reduced_number_operators, tapering
t0 = time.time()
d, P_RO, SIG = 2.4, 0.02, 0.5e-3
m = molecule(d); H8 = qubit_hamiltonian(m); H5 = tapering(m, H8).apply(H8); ops = reduced_number_operators(m)
ans = build_compact_ansatz(m); x = optimize_compact(ans, H5, ops["S2"])["x"]
bound = ans.circuit.assign_parameters(x); isa = to_line(bound)
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
rows = []
def add(name, base, lams, rescale=False):
    E = np.array([est[base][l]["E"] for l in lams]); K = np.array([est[base][l]["K"] for l in lams])
    if rescale:   # eco: E' = (E − (1 − F^λ) E_mix) / F^λ
        fl = F ** np.asarray(lams, float); E, K = (E - (1 - fl) * emix) / fl, K / fl ** 2
    g = gam(lams)
    rows.append({"técnica": name, "sesgo (mHa)": 1e3 * (g @ E - E0),
                 "shots para σ = 0,5 mHa (M)": float(np.sum(np.abs(g) * np.sqrt(K)) ** 2) / SIG ** 2 / 1e6})
for k in tabs: add(k, k, (1,))
add("ZNE lineal (1,3) + PS + RO", "PS + RO", (1, 3))
add("eco + ZNE lineal (1,3) + RO", "RO", (1, 3), rescale=True)
add("ZNE Richardson (1,3,5) + PS + RO", "PS + RO", (1, 3, 5))
print(f"d = {d} Å | {int(isa.count_ops()['cz'])} CZ | F(eco) = {F:.4f} | CCSD − FCI = {1e3*(m['e_ccsd']-m['e_fci']):.1f} mHa | {time.time()-t0:.0f} s")
pd.DataFrame(rows).round(3)'''

CHECKS = '''\
import numpy as np, importlib, time
import h4vqe
for mod in ["config","hamiltonian","ansatz","optimize","estimator","adapt","routing","analysis","pipeline"]:
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


def md(name: str):
    return nbf.v4.new_markdown_cell((MD / name).read_text(encoding="utf-8"))


def code(src: str):
    return nbf.v4.new_code_cell(src)


def _modules(mods) -> list:
    return [code(f"%%writefile h4vqe/{m}.py\n" + (PKG / f"{m}.py").read_text(encoding="utf-8")) for m in mods]


def _finish(cells) -> nbf.NotebookNode:
    nb = nbf.v4.new_notebook()
    nb.cells = cells
    nb.metadata["kernelspec"] = {"name": "python3", "display_name": "Python 3", "language": "python"}
    nb.metadata["colab"] = {"provenance": []}
    todo = [c.source[:80] for c in cells if c.cell_type == "markdown" and "TODO" in c.source]
    assert not todo, f"quedó un TODO en el texto: {todo}"
    return nb


def build_core(quick: bool = False) -> nbf.NotebookNode:
    """Lo que pide la consigna: curva con energías medidas, R_e, D_e y confiabilidad."""
    run = RUN.replace("QUICK = False", "QUICK = True ") if quick else RUN
    return _finish([md("00_intro.md"), md("01_setup.md"), code(INSTALL), *_modules(CORE_MODULES),
                    md("02_hamiltonian.md"), code(HAMILTONIAN),
                    md("03_ansatz.md"), md("04_measurement.md"), code(CHECKS),
                    md("05_results.md"), code(run), code(TABLE), code(PLOT), code(RE_DE),
                    md("06_reliability.md"), code(RELIABILITY), code(DEPTH),
                    md("08_decisions.md")])


def build_extras() -> nbf.NotebookNode:
    """Extras: menos qubits, mitigación, ruido de hardware (FakeFez) y otra geometría."""
    return _finish([md("00_extras.md"), md("01_setup.md"), code(INSTALL_EXTRAS), *_modules(EXTRA_MODULES),
                    md("02_reduction.md"), code(REDUCTION),
                    md("07_ansatz_comparison.md"),
                    md("07_mitigation.md"), code(MITIGATION),
                    md("07_hardware.md"), code(HARDWARE),
                    md("07_geometry.md"), code(GEOMETRY)])


BUILD = {"core": build_core, "extras": build_extras}


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


def update_text(which: str) -> None:
    """Reemplaza solo las celdas de texto de un notebook ya ejecutado (las salidas quedan)."""
    old, new = nbf.read(OUT[which], as_version=4), BUILD[which]()
    if [c.cell_type for c in old.cells] != [c.cell_type for c in new.cells]:
        raise ValueError(f"{which}: cambió la estructura de celdas, hay que volver a ejecutar con --run")
    for c_old, c_new in zip(old.cells, new.cells):
        if c_new.cell_type == "markdown":
            c_old.source = c_new.source
        elif c_old.source != c_new.source:
            raise ValueError(f"{which}: cambió una celda de código, hay que volver a ejecutar con --run")
    nbf.write(old, OUT[which])
    print("texto actualizado en", OUT[which])


if __name__ == "__main__":
    args = sys.argv[1:]
    which = [w for w in ("core", "extras") if w in args] or ["core", "extras"]
    kernel = args[args.index("--kernel") + 1] if "--kernel" in args else "h4vqe-local"
    for w in which:
        if "--text-only" in args:
            update_text(w)
            continue
        nb = build_core(quick="--quick" in args) if w == "core" else build_extras()
        nbf.write(nb, OUT[w])
        print("escrito", OUT[w])
        if "--run" in args:
            nbf.write(execute(nb, kernel), OUT[w])
