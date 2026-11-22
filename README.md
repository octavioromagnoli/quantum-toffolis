# H₄ con VQE: la curva de disociación con precisión química

[![Abrir en Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/octavioromagnoli/quantum-toffolis/blob/main/notebook/H4_VQE.ipynb)

Proyecto del equipo **quantum toffolis** para la Quantum Jam 2026. Con VQE (Qiskit + PySCF) calculamos la curva de energía de la cadena lineal H₄ (base STO-3G, 8 qubits) entre 0,6 y 3,0 Å. Donde los métodos clásicos aproximados empiezan a fallar, al estirar los enlaces, el resultado le gana a CCSD.

Todas las energías se **miden** con `SamplerV2` de `qiskit_aer`, con shots. No se lee el estado del simulador. El notebook corre de punta a punta en un Colab nuevo, instala sus dependencias y no usa datos precalculados.

## Resultados

| | VQE (este trabajo) | referencia |
|---|---|---|
| Error en la zona estirada (d ≥ 1,8 Å) | ≤ 0,56 mHa, σ ≈ 0,5 mHa | CCSD: 9–31 mHa |
| Precisión química (< 1,6 mHa) | 18 de 18 distancias | — |
| $R_e$ | 0,8873 ± 0,0018 Å | FCI: 0,8886 Å |
| $D_e$ | 313,00 ± 0,66 mHa | FCI: 313,21 mHa |
| Compuertas CZ (todos con todos / en línea) | 194 / 378 | UCCSD ×2: 582 / 1648 |
| Con ruido (0,1 % por CZ), 2,4 Å | 5 qubits, 16 CZ, sesgo mitigado de +0,008 mHa | UCCSD ×2 sin mitigar: 548 mHa |
| Ruido de un chip de IBM (FakeFez), 2,4 Å | −0,71 ± 1,23 mHa | CCSD: −30,6 mHa |

## Cómo llegamos

El [notebook](notebook/H4_VQE.ipynb) cuenta el recorrido. Hay también una [versión con texto corto](notebook/H4_VQE_corto.ipynb). Arranca de lo mínimo que pide la consigna, un VQE con UCCSD, y resuelve los problemas en el orden en que aparecieron:

1. **Demasiadas CZ.** ADAPT-VQE con una lista común de 15 excitaciones, más un ruteo propio con fSWAP para qubits en línea, baja de 1648 a 378 CZ.
2. **Demasiados shots y una σ dudosa.** Agrupamiento QWC con *sorted insertion* y reparto de Neyman (hasta 1,9× menos shots). La barra de error se calibró con 3600 repeticiones.
3. **Sesgo en el ajuste de $R_e$.** El método de ajuste se eligió midiendo su sesgo contra el mínimo FCI exacto.
4. **Confiabilidad sin mirar FCI.** Varianza de la energía por eco de Loschmidt, simetrías $N_\alpha$, $N_\beta$ y $S^2$ medidas, y convergencia con la profundidad del ansatz.
5. **El ruido.** Tapering por simetrías (8 → 5 qubits), un ansatz compacto de 16 CZ y mitigación (lectura, post-selección y ZNE).

Extras: FakeFez, otras geometrías (H₄ cuadrado y en zigzag) y una comparación de ansätze con y sin ruido. La consigna completa está en [CONSIGNA.md](CONSIGNA.md).

## Correrlo

**En Colab:** abrí el badge de arriba y hacé *Entorno de ejecución → Ejecutar todo*. La curva completa tarda unos 25 min en Colab gratis.

**En local** (Python 3.12, con [uv](https://docs.astral.sh/uv/)):

```bash
uv venv --python 3.12 .venv
uv pip install --python .venv/bin/python --only-binary pyscf "qiskit==2.6.2" "qiskit-aer==0.18.4" "pyscf==2.15.0" "qiskit-ibm-runtime==0.52.1" scipy matplotlib pandas nbformat nbconvert ipykernel pytest
.venv/bin/python -m ipykernel install --sys-prefix --name h4vqe-local
.venv/bin/python -m pytest -q tests                        # ~30 s
.venv/bin/python notebook/build_notebook.py --run          # genera y ejecuta el notebook en una carpeta vacía
.venv/bin/python notebook/build_notebook.py --quick --run  # curva corta, para probar rápido
```

En macOS, `--only-binary pyscf` evita compilar PySCF, que falla. `qiskit-ibm-runtime` solo hace falta para FakeFez.

El `.ipynb` se genera a partir de `h4vqe/*.py` y `notebook/md/*.md`, así que no se edita a mano: hay que cambiar los módulos o el texto y regenerarlo. Con `--text-only` se actualiza solo el texto del notebook ya ejecutado.

## Estructura

```
h4vqe/                 el paquete
  hamiltonian.py       PySCF -> Jordan-Wigner; HF, FCI y CCSD de referencia
  ansatz.py            excitaciones (UCCSD) y simulador numpy exacto
  adapt.py, routing.py ADAPT-VQE y ruteo fSWAP para qubits en línea
  optimize.py          búsqueda de θ sin ruido, con gradiente analítico
  estimator.py         energía medida con SamplerV2 de Aer, con su σ
  analysis.py          tabla, R_e, D_e, gráficos, varianza y eco
  pipeline.py          curva con warm start y diagnósticos de confiabilidad
  reduction.py         tapering 8 -> 5 qubits y ansatz compacto de 16 CZ
  noise.py             modelo de ruido, lectura y ZNE
  mitigation.py        post-selección, lectura combinable y ZNE
  hardware.py          FakeFez
  geometry.py          H4 cuadrado y en zigzag
  historia.py          diagnósticos del relato del notebook
  scaling.py           escalado H4 -> H14 -> FeMoco (exploración, fuera del notebook)
notebook/              H4_VQE.ipynb, su versión corta, el generador y el texto de cada sección
scripts/               estudios largos que el notebook no recalcula (calibración de σ, mitigación con shots, ...)
results/               resultados de esos estudios
tests/                 tests con pytest
presentacion/          diapositivas
```

## Equipo

Octavio, Gonzalo, Santino y Joaquín. Cómo nos repartimos el trabajo está en [PLAN.md](PLAN.md).
