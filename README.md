# quantum-toffolis — H₄ con VQE (Quantum Jam 2026)

- Consigna: [CONSIGNA.md](CONSIGNA.md)
- **Plan y tareas por persona: [PLAN.md](PLAN.md)**
- Notebook de entrega (generado, no editar a mano): [notebook/H4_VQE.ipynb](notebook/H4_VQE.ipynb)

## Estructura

```
h4vqe/                 paquete (un dueño por módulo, ver PLAN.md)
  config.py            distancias, shots, ansatz por defecto
  hamiltonian.py       PySCF -> Jordan-Wigner, HF/FCI/CCSD            (Octavio)
  ansatz.py            UCCSD con Trotter + simulador numpy exacto      (Gonzalo)
  optimize.py          θ* sin ruido, L-BFGS-B con gradiente analítico  (Gonzalo)
  estimator.py         energía medida con SamplerV2 de Aer, ±σ         (Santino)
  analysis.py          tabla, R_e, D_e, gráficos, diagnósticos         (Joaquín)
  pipeline.py          barrido de la curva con warm start              (Joaquín)
notebook/
  md/*.md              texto de cada sección del notebook
  build_notebook.py    arma H4_VQE.ipynb (módulos con %%writefile + texto)
tests/test_smoke.py    chequeos de cordura (~10 s)
```

## Arranque rápido

```bash
uv venv --python 3.12 .venv
uv pip install --python .venv/bin/python --only-binary pyscf "qiskit==2.6.2" "qiskit-aer==0.18.4" "pyscf==2.15.0" scipy matplotlib pandas nbformat nbconvert ipykernel pytest
.venv/bin/python -m ipykernel install --sys-prefix --name h4vqe-local
.venv/bin/python tests/test_smoke.py
.venv/bin/python notebook/build_notebook.py          # solo genera
.venv/bin/python notebook/build_notebook.py --run    # genera y ejecuta (~10 min)
```

En macOS, `--only-binary pyscf` evita compilar PySCF, que falla.

## Probar en Colab

```bash
uv tool install google-colab-cli
colab new -s h4
colab exec -s h4 -f notebook/H4_VQE.ipynb
colab stop -s h4
```

Con nuestra cuenta, el login de `colab` (OAuth) devolvió solo los permisos de email y perfil, sin el de Colab, y no anduvo. Lo que sí funciona: subir el `.ipynb` a Drive, abrirlo con Colab y hacer "Entorno de ejecución → Ejecutar todo". La curva completa tarda ~25 min en Colab gratis.
