# Plan de trabajo — H₄ con VQE (5 horas, 4 personas)

## Dónde estamos (base ya funcionando)

La base en `main` **ya cumple el desafío central sin ruido**. Corre de punta a punta y genera el notebook:

| qué | estado |
|---|---|
| Hamiltoniano H₄ (PySCF → Jordan-Wigner, 8 qubits, 185 Paulis) + HF/FCI/CCSD | ✅ |
| Ansatz UCCSD ×2 (52 parámetros): error < 0,001 mHa contra FCI en toda la curva (sin ruido) | ✅ |
| Optimización sin ruido con gradiente analítico (segundos por punto) | ✅ |
| Energía **medida** con `SamplerV2` de Aer: 46 grupos QWC, piloto + Neyman, σ desde las muestras | ✅ |
| Tabla (E, σ, shots, CZ), gráfico, R_e y D_e con incertidumbre Monte Carlo | ✅ |
| CCSD falla en la zona estirada: entre 9 y 31 mHa *por debajo* de FCI | ✅ ya lo ganamos |

**Probado en Colab limpio (Python 3.13, 2026-11-21):** corre de punta a punta sin errores y da los mismos números (R_e = 0,8936 ± 0,0019 Å, D_e = 312,95 ± 0,84 mHa, 7/7 contra CCSD). **Tarda ~25 min** (la curva sola: 1479 s, ~80 s por punto, contra ~12 s en una Mac). Después de esa corrida, el estimador se aceleró ~25 % (`get_counts` → numpy). Bajar shots (Santino) y CZ (Gonzalo) también reduce el tiempo en Colab.

**Puntos flojos** (en eso se va el tiempo):
1. **Demasiadas CZ**: UCCSD ×2 tiene ~584 CZ con conectividad total y ~1600 en línea. Con el ruido de la consigna (0,1 % por CZ) eso es inútil, y el jurado mira las CZ.
2. **Demasiados shots**: 8 M por punto para σ ≈ 0,7 mHa. Hay que bajar ese número o justificarlo.
3. **Calibración de σ**: la dispersión de z da 1,16 en 60 repeticiones (debería ser 1). Hay que confirmar que es estadística o corregirlo.
4. Extras: ruido, mitigación, 6 qubits, comparar ansätze.

## Cómo trabajamos sin pisarnos

```bash
git pull origin main
uv venv --python 3.12 .venv && uv pip install --python .venv/bin/python --only-binary pyscf "qiskit==2.6.2" "qiskit-aer==0.18.4" "pyscf==2.15.0" scipy matplotlib pandas nbformat nbconvert ipykernel pytest
.venv/bin/python tests/test_smoke.py           # tiene que imprimir OK
.venv/bin/python notebook/build_notebook.py    # genera notebook/H4_VQE.ipynb y notebook/H4_VQE_extras.ipynb
```

- **Cada uno toca solo sus archivos** (tabla de abajo). `config.py` y `pipeline.py` los integra Joaquín: si necesitan cambiarlos, avisen.
- **Nadie edita los `.ipynb` a mano.** Se genera desde `h4vqe/*.py` y `notebook/md/*.md`. Cada uno escribe su sección en su `.md`.
- Una rama por persona (`octavio/...`, `gonzalo/...`, `santino/...`, `joaquin/...`) y un PR chico a `main` cada vez que algo funciona. Mergear seguido, no al final.
- Antes de cada PR corré `tests/test_smoke.py`.
- **No romper interfaces**: `build_ansatz(m, excitations=..., reps=...) -> Ansatz` (con `.circuit`, `.fast`), `measure_energy(circ, H, shots) -> {"energy","sigma","shots",...}`, `molecule(d) -> dict`, `qubit_hamiltonian(m)`.

| persona | archivos propios | sección del notebook |
|---|---|---|
| Octavio | `h4vqe/hamiltonian.py`, `h4vqe/reduction.py` (nuevo), `h4vqe/noise.py` (nuevo) | `02_hamiltonian.md`, extras de ruido y 6 qubits |
| Gonzalo | `h4vqe/ansatz.py`, `h4vqe/optimize.py`, `h4vqe/adapt.py` (nuevo) | `03_ansatz.md`, extra de comparar ansätze |
| Santino | `h4vqe/estimator.py`, `h4vqe/mitigation.py` (nuevo) | `04_measurement.md`, extra de mitigación |
| Joaquín | `h4vqe/analysis.py`, `h4vqe/pipeline.py`, `h4vqe/config.py`, `notebook/*` | `00`, `05`, `06`, `08`, integración y Colab |

## Cronograma (T = hora desde que arrancamos)

| T | todos |
|---|---|
| 0:00–0:20 | pull, instalar, correr el test y el notebook en modo `QUICK`, leer su módulo |
| 0:20–2:00 | **bloque 1**: tareas prioritarias (P1) de cada uno → PR |
| 2:00 | **sync de 10 min**: Gonzalo pasa la nueva cantidad de CZ, Santino los nuevos shots/σ, Octavio el modelo de ruido |
| 2:10–3:30 | **bloque 2**: P2 y extras (el ruido ya con el ansatz corto de Gonzalo) |
| 3:30 | **congelamiento de código**: solo arreglos. Joaquín integra y corre la curva completa |
| 3:30–4:15 | cada uno termina su `.md` (decisiones, qué no funcionó, limitaciones) con los números finales |
| 4:15–4:45 | Joaquín corre el notebook completo en Colab limpio; los demás revisan el notebook ejecutado |
| 4:45–5:00 | subir a la carpeta de Drive compartida con el jurado y verificar los permisos |

## Tareas

Cada uno tiene su issue asignado con checklist: Octavio [#1](https://github.com/octavioromagnoli/quantum-toffolis/issues/1), Gonzalo [#2](https://github.com/octavioromagnoli/quantum-toffolis/issues/2), Santino [#3](https://github.com/octavioromagnoli/quantum-toffolis/issues/3), Joaquín [#4](https://github.com/octavioromagnoli/quantum-toffolis/issues/4).

### Octavio: hamiltoniano, simetrías y ruido
- **P1** Validar el hamiltoniano: diagonalizar Ĥ solo en el sector N_α = N_β = 2 y comparar con FCI en todas las distancias. Explicar en `02_hamiltonian.md` por qué CCSD falla (no variacional; correlación fuerte al estirar), con un gráfico de E_CCSD − E_FCI.
- **P1** `h4vqe/noise.py`: modelo de ruido de la consigna (`ruido.py`: cz 0,1 %, 1 qubit 0,005 %, base cz/sx/x/rz, qubits en línea). Buscar `ruido.py` en el repo de workshops; si no aparece, implementarlo con `NoiseModel` y despolarizante. Devolver un `AerSampler` con ruido que se pueda pasar al estimador.
- **P2** Extra **menos qubits**: mapeo de paridad + reducción de 2 qubits (y tapering con simetrías Z₂ si da el tiempo) → 6 qubits. Medir cuántos grupos y CZ se ahorran. Ojo: el ansatz se tiene que mapear igual (generadores con el mismo mapeo y la misma reducción).
- **P2** Extra **ruido de hardware**: un punto estirado (2,4 Å) con `FakeFez`.
- **Extra** otra geometría (H₄ cuadrado o zigzag), cambiando solo `h4_geometry`.

### Gonzalo: ansatz con menos CZ
- **P1** `h4vqe/adapt.py`: **ADAPT-VQE** con pool fermiónico (las 26 excitaciones de UCCSD, que pueden repetirse). En cada paso se elige el operador de gradiente máximo con `FastAnsatz` y se para con error < 0,5 mHa o gradiente < 1e-4. Devuelve la lista de excitaciones → `build_ansatz(m, excitations=lista)`. Meta: precisión química en toda la curva con muchas menos CZ que las ~584 de UCCSD×2.
- **P1** Una sola lista de operadores para toda la curva (o por zona), justificada: el notebook tiene que reportar las CZ por punto.
- **P2** Excitaciones de qubits (QEB, sin cadenas Z: simple = 2 CNOT, doble = 13 CNOT) como pool alternativo. `FastAnsatz` sigue sirviendo, porque G³ = G también vale.
- **P2** Extra **comparar ansätze**: tabla con UCCSD×1, UCCSD×2, ADAPT fermiónico, ADAPT-QEB y hardware-efficient (`efficient_su2`): error máximo en la zona estirada, parámetros, CZ todos con todos y en línea, profundidad.
- Avisar en el sync de las 2:00 qué ansatz queda como default (`config.ANSATZ_REPS` / lista).

### Santino: mediciones, σ y mitigación
- **P1** Calibración de σ: ≥200 repeticiones a 2,0 Å y a 3,0 Å, histograma de z = (E − E_exacta)/σ. Explicar o corregir la dispersión de 1,16. Agregar σ por **bootstrap** como chequeo cruzado.
- **P1** Bajar shots para el mismo σ: (a) agrupamiento por conmutación general (no QWC) con diagonalización Clifford; (b) comparar las asignaciones Σ|c|, Neyman y Neyman con varianza teórica del estado HF; (c) curva σ contra shots y shots necesarios para σ = 0,5 mHa. Si conviene, cambiar `SHOTS_PER_POINT` (avisar).
- **P2** `h4vqe/mitigation.py`: con el ruido de Octavio, (a) corrección de **readout**; (b) **post-selección por simetría**: en los grupos que miden solo Z, descartar los bitstrings con N_α ≠ 2 o N_β ≠ 2; (c) **ZNE** con plegado de CZ (factores 1, 3, 5) y extrapolación lineal o de Richardson, con σ propagada.
- **P2** Medir cuánto mejora cada técnica en 1 o 2 puntos estirados.

### Joaquín: integración, análisis, notebook y Colab
- **P1** Integrar PRs y mantener `main` verde (test + notebook en modo QUICK).
- **P1** R_e y D_e robustos: comparar el ajuste cuadrático contra el cúbico o spline y los puntos usados; reportar R_e(VQE) contra R_e(FCI) y D_e(VQE) contra D_e(FCI) con σ.
- **P1** Sección de confiabilidad (`06_reliability.md`): varianza ⟨H²⟩ − ⟨H⟩² **medida con shots** (no con el estado), simetrías medidas, convergencia con la profundidad o el número de operadores de ADAPT, y la calibración de Santino.
- **P2** Extra **VQE contra SQD** en la zona estirada, si queda tiempo.
- **Cierre** Notebook completo en Colab limpio (Entorno de ejecución → Reiniciar y ejecutar todo), tabla final, `00_intro.md` con los números, subida a Drive compartida con el jurado.

## Definición de terminado
- [ ] El notebook corre en Colab limpio, de punta a punta, en menos de ~20 min (hoy: ~25 min, ver arriba)
- [ ] 18 distancias con E ± σ, shots y CZ por punto
- [ ] |E_VQE − E_FCI| + σ < |E_CCSD − E_FCI| para d ≥ 1,8 Å
- [ ] |E_VQE − E_FCI| < 1,6 mHa en la zona estirada
- [ ] R_e y D_e con incertidumbre, comparados con FCI
- [ ] Argumento de confiabilidad sin usar FCI
- [ ] Decisiones, qué no funcionó y limitaciones
- [ ] Al menos 2 extras (objetivo: ruido + mitigación, y menos CZ o comparar ansätze)
