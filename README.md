# quantum-toffolis — H₄ con VQE (Quantum Jam 2026)

- Consigna: [CONSIGNA.md](CONSIGNA.md)
- **Plan y tareas por persona: [PLAN.md](PLAN.md)**
- Notebook de entrega (generado, no editar a mano): [notebook/H4_VQE.ipynb](notebook/H4_VQE.ipynb).
  Uno solo, como pide la consigna, contado como historia: de lo mínimo (VQE con UCCSD) a la versión
  final (ADAPT, mediciones calibradas, confiabilidad medida) y todos los extras (ruido, 5 qubits,
  mitigación, FakeFez, otra geometría, comparación de ansätze y escalado a FeMoco), con un resumen
  final que compara lo mínimo contra lo final.

## Estructura

```
h4vqe/                 paquete (un dueño por módulo, ver PLAN.md)
  config.py            distancias, shots, ansatz de la curva (ANSATZ = "adapt")
  hamiltonian.py       PySCF -> Jordan-Wigner, HF/FCI/CCSD, validación      (Octavio)
  ansatz.py            excitaciones (UCCSD, listas) + simulador numpy exacto (Gonzalo)
  adapt.py             ADAPT-VQE y listas guardadas (ADAPT_FERMIONIC, 15 ops) (Gonzalo)
  routing.py           ruteo fSWAP para qubits en línea (378 CZ)            (Gonzalo)
  optimize.py          θ* sin ruido, L-BFGS-B con gradiente analítico       (Gonzalo)
  estimator.py         energía medida con SamplerV2 de Aer, ±σ              (Santino)
  analysis.py          tabla, R_e, D_e, gráficos, varianza y eco            (Joaquín)
  pipeline.py          curva con warm start + diagnósticos de confiabilidad (Joaquín)
  reduction.py         tapering 8 -> 5 qubits y ansatz compacto de 16 CZ    (Octavio)
  noise.py             modelo de ruido propio, lectura y ZNE                (Octavio)
  hardware.py          extra FakeFez (necesita qiskit-ibm-runtime)          (Octavio)
  geometry.py          extra H4 cuadrado y en zigzag                        (Octavio)
  mitigation.py        post-selección, lectura combinable y ZNE plegando R·C (Santino)
  scaling.py           extra de escalado H4 -> H14 -> FeMoco                (Santino)
  historia.py          diagnósticos del relato: lo mínimo, ADAPT, shots, ruido por ansatz (Santino)
notebook/
  md/*.md              texto de cada sección del notebook
  build_notebook.py    arma H4_VQE.ipynb (módulos y datos con %%writefile + texto)
scripts/               estudios largos que el notebook no recalcula (calibración de σ, mitigación con shots,
                       comparación de ansätze, escalado); sus resultados van a results/
tests/                 test_smoke, test_routing, test_analysis, test_reduction, test_hardware, test_scaling
```

## Arranque rápido

```bash
uv venv --python 3.12 .venv
uv pip install --python .venv/bin/python --only-binary pyscf "qiskit==2.6.2" "qiskit-aer==0.18.4" "pyscf==2.15.0" "qiskit-ibm-runtime==0.52.1" scipy matplotlib pandas nbformat nbconvert ipykernel pytest
.venv/bin/python -m ipykernel install --sys-prefix --name h4vqe-local
.venv/bin/python -m pytest -q tests                  # ~30 s
.venv/bin/python notebook/build_notebook.py                # genera el notebook
.venv/bin/python notebook/build_notebook.py --run          # genera y ejecuta en una carpeta vacía (como un Colab nuevo)
.venv/bin/python notebook/build_notebook.py --quick --run  # curva corta, para probar rápido
.venv/bin/python notebook/build_notebook.py --text-only    # solo actualiza el texto del ya ejecutado
```

En macOS, `--only-binary pyscf` evita compilar PySCF, que falla. `qiskit-ibm-runtime` solo hace falta para el extra de FakeFez (`hardware.py`, `tests/test_hardware.py`).

El kernel `h4vqe-local` es uno solo por usuario: si hay varios worktrees, apunta al `.venv` del último que lo registró. Para ejecutar el notebook desde otro worktree, registrá un kernel con otro nombre (`ipykernel install --sys-prefix --name ...`) y pasalo con `build_notebook.py --run --kernel NOMBRE`.

## Probar en Colab

```bash
uv tool install google-colab-cli
colab new -s h4
colab exec -s h4 -f notebook/H4_VQE.ipynb
colab stop -s h4
```

Con nuestra cuenta, el login de `colab` (OAuth) devolvió solo los permisos de email y perfil, sin el de Colab, y no anduvo. Lo que sí funciona: subir el `.ipynb` a Drive, abrirlo con Colab y hacer "Entorno de ejecución → Ejecutar todo". La curva completa tarda ~25 min en Colab gratis.
