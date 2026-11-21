# Extra: ¿cuánto costaría la molécula real? (escalado de H₄ a FeMoco)

> **Estado:** idea para implementar **al final**. Es un extra (§05 de la consigna) y solo se evalúa si el desafío central está completo. No arrancar hasta tener terminados la curva de H₄, R_e, D_e y la comparación con CCSD.

## 1. La idea en una línea

Medir cómo crecen los recursos de VQE al alargar la cadena de hidrógenos (H₄, H₆, …, H₁₆) y cuánto cuesta simular eso en una computadora clásica. Después, extrapolar ambas curvas al tamaño del **FeMoco** (el núcleo de la nitrogenasa, que es la motivación de la consigna) para mostrar con números hasta dónde llega el método y por qué.

## 2. Por qué suma

- Cierra la historia de la consigna: nitrogenasa → H₄ → de vuelta a la nitrogenasa.
- Responde con números el criterio de **análisis de limitaciones**: hasta dónde funciona el método, dónde falla y por qué.
- Suma en **creatividad** y es poco probable que otro equipo lo haga.
- Es barato: la mayor parte consiste en **contar** recursos, no en simular.

## 3. Contexto: la molécula real

No se simula la proteína entera (decenas de miles de átomos), sino el cofactor **FeMoco** (núcleo MoFe₇S₉C, ~18 átomos en el cúmulo central). Y de él, solo un *espacio activo* de orbitales:

| Modelo | Electrones / orbitales | Qubits (espín-orbitales) |
|---|---|---|
| Reiher et al. 2017 | 54e / 54o | 108 |
| Li & Chan 2019 (sostienen que el de 54 no es representativo) | 113e / 76o | 152 |
| **Nuestra H₄ (STO-3G)** | 4e / 4o | 8 |

- **Vía clásica:** una función de onda para el modelo (54, 54) requeriría ~10²⁹ variables (Montgomery & Mazziotti).
- **Vía cuántica con corrección de errores** (Lee et al. 2021, QPE + tensor hypercontraction, **no VQE**): ~2.200 qubits lógicos, ~3,2×10¹⁰ Toffolis, **~4 millones de qubits físicos y menos de 4 días** (con ciclos de 1 µs).
- **VQE en general:** Gonthier et al. 2022 muestran que la cantidad de mediciones es el cuello de botella para llegar a precisión química en moléculas de tamaño real, aun con buenas técnicas de agrupamiento.

## 4. Diseño

### 4.1 La perilla: cadenas H_N

- Cadena lineal equiespaciada en STO-3G, igual que en el desafío central.
- **N par:** 4, 6, 8, 10, 12, 14, 16. Con N impar la capa queda abierta y la comparación se complica.
- N átomos → N orbitales → **2N qubits** (Jordan-Wigner) y N electrones. Hay un electrón por orbital, igual que en el modelo (54e, 54o) de Reiher, así que **las mismas fórmulas valen para N = 54**.
- **Dos distancias:**
  - **d_eq:** el R_e de H₄ que salga del desafío central.
  - **d_est = 2,0 Å:** zona estirada.

  La cantidad de qubits y de compuertas no depende de d. Los shots y el error de CCSD sí.

### 4.2 Qué se mide para cada N

| Magnitud | Cómo | Hasta N | Dónde |
|---|---|---|---|
| Qubits | 2N (fórmula) | todos | — |
| Dimensión FCI | C(N, N/2)² (fórmula) | todos | — |
| Memoria del statevector | 16 · 2^(2N) bytes (fórmula) | todos | — |
| Términos de Pauli de Ĥ | construir Ĥ y contar | 16 | Colab (≤12), VM (14–16) |
| Grupos de medición (QWC) | `H.group_commuting(qubit_wise=True)` | 12 (14–16 si entra) | Colab / VM |
| Σ\|c_k\| (sin la identidad) | de Ĥ | 16 | Colab / VM |
| Parámetros UCCSD | fórmula (§4.3); chequear con `saltos_posibles` para N ≤ 8 | todos | — |
| CZ del ansatz | transpilar a `cz, sx, x, rz` con conectividad en línea | 10 completo; 12–16 por excitación | Colab / VM |
| P(circuito sin errores) con `ruido.py` | 0,999^(n_CZ) · 0,99995^(n_1q) | todos | — |
| Shots para σ = 1,6 mHa | cota (§4.3) y valor exacto para N chico | todos | Colab / VM |
| E_HF, E_CCSD, E_FCI | PySCF | FCI hasta 16 | Colab (≤12), VM (14–16) |
| Error de CCSD | \|E_CCSD − E_FCI\| | 16 | ídem |
| VQE real (E, σ, shots, CZ) | el mismo pipeline del desafío central | 6, 8 (10 opcional) | Colab |
| Tiempo y memoria de simulación clásica | benchmark (§6, fase 3b) | 20–32 qubits | VM |

### 4.3 Fórmulas

**Parámetros de UCCSD.** Sea o_σ / v_σ la cantidad de orbitales ocupados / vacíos con espín σ:

- simples = o_α·v_α + o_β·v_β
- dobles αβ = o_α·v_α·o_β·v_β
- dobles αα + ββ = C(o_α, 2)·C(v_α, 2) + C(o_β, 2)·C(v_β, 2)

Para H_N vale o = v = N/2 en ambos espines. Con N = 4 da 26, que tiene que coincidir con lo que devuelva `saltos_posibles` en el desafío central.

**Shots por evaluación de energía.** Si el grupo g se mide con n_g shots, la varianza del estimador es σ² = Σ_g Var_g / n_g. Con el reparto óptimo (n_g ∝ √Var_g), los shots totales para llegar a σ = ε son

  M = (Σ_g √Var_g)² / ε².

Como √Var_g ≤ Σ_{k∈g} |c_k|, se obtiene la cota del peor caso

  M ≤ (Σ_{k≠I} |c_k|)² / ε².

- Para N ≤ 12 se calcula también M **exacto**, usando el estado del VQE o el de FCI como statevector, y se compara con la cota.
- Para N grande se usa la cota, o la cota corregida por el cociente exacto/cota ajustado en N chico.
- Ojo: esto es **por evaluación de energía**. Un VQE completo lo multiplica por la cantidad de evaluaciones del optimizador.

**Probabilidad de que el circuito no tenga errores con `ruido.py`:** P = (1 − 0,001)^(n_CZ) · (1 − 0,00005)^(n_1q).

### 4.4 Tabla de referencia (ya calculada)

| N | Qubits | Dimensión FCI | Parámetros UCCSD | Statevector (complex128) |
|---|---|---|---|---|
| 4 | 8 | 36 | 26 | 4 KiB |
| 6 | 12 | 400 | 117 | 64 KiB |
| 8 | 16 | 4.900 | 360 | 1 MiB |
| 10 | 20 | 63.504 | 875 | 16 MiB |
| 12 | 24 | 853.776 | 1.818 | 256 MiB |
| 14 | 28 | 11.778.624 | 3.381 | 4 GiB |
| 16 | 32 | 165.636.900 | 5.792 | 64 GiB |
| 54 (FeMoco, Reiher) | 108 | ≈ 3,8×10³⁰ | 779.301 | ≈ 5×10³³ bytes |

Colab tiene ~12 GB de RAM, así que el statevector entra hasta H₁₄ (en ajustado) y H₁₆ solo en una VM.

## 5. Por qué no hacemos VQE completo en H₁₆ (y cómo lo convertimos en resultado)

- 32 qubits entran en memoria en una VM de 128 GB. El problema es el **tiempo**:
  - Cada compuerta recorre todo el statevector (64 GiB), lo que tarda del orden de 1 s en CPU.
  - UCCSD para H₁₆ tiene 5.792 excitaciones, del orden de 10⁵ o más compuertas de dos qubits, así que **una sola evaluación de energía tarda horas**. A eso se suma medir decenas de miles de términos de Pauli.
  - Solo un gradiente por *parameter shift* pide ~2 × 5.792 evaluaciones. Un VQE completo llevaría meses o años.
- **No lo estimamos a ojo, lo medimos:** el benchmark de la fase 3b da el tiempo de una evaluación con 32 qubits, y lo multiplicamos por la cantidad de evaluaciones de un VQE. Ese número es un resultado del extra.
- Este es el argumento de fondo de la computación cuántica: simular clásicamente un sistema cuántico cuesta ~2^n.

## 6. Plan de implementación

### Fase 0: prerrequisitos
- El desafío central está terminado.
- Se reutilizan sus funciones (molécula, hamiltoniano, ansatz, conteo de CZ) **con las mismas convenciones**: orden de espín-orbitales α en un bloque y β en otro, little-endian y base `cz, sx, x, rz` con conectividad en línea.
- Sanity check obligatorio: para N = 4, todos los conteos del extra tienen que coincidir exactamente con los del desafío central.

### Fase 1: conteo de recursos (Colab, sin GCP) · ~3–4 h
1. Función `hidrogeno(N, d)` que arme la geometría `"H 0 0 0; H 0 0 d; H 0 0 2d; …"`.
2. Para N = 4…12: construir Ĥ y contar términos, grupos QWC y Σ|c_k|.
3. Parámetros UCCSD con la fórmula, verificados contra `saltos_posibles` para N ≤ 8.
4. CZ:
   - Para N ≤ 10, transpilar el circuito completo.
   - Para N ≥ 12, transpilar cada excitación (`PauliEvolutionGate`) por separado con la misma conectividad y sumar. Es una cota superior porque no hay cancelaciones entre excitaciones, y hay que aclararlo.
5. Calcular P(sin errores) y la cota de shots.
6. Salida: `resultados/recursos.csv`.

**Cuidado con la construcción de Ĥ:** `hamiltoniano_qubits` (de la skill) hace productos de `SparsePauliOp` en bucles de orden N⁴ y se vuelve lenta desde ~24 qubits. Si no da, usar OpenFermion (`jordan_wigner` sobre un `InteractionOperator`) y convertir a `SparsePauliOp`. OpenFermion **intercala** α y β (α₀, β₀, α₁, …), así que hay que reordenar a bloques para mantener la convención. Verificar con N = 4 que ambas construcciones dan el mismo espectro.

### Fase 2: VQE real en H₆ y H₈ (Colab) · ~3–4 h
- El mismo pipeline del desafío central, en d_eq y d_est. Las energías salen de medir con EstimatorV2/SamplerV2 de `qiskit_aer`, igual que en el central.
- Registrar E, σ, shots, CZ, error contra FCI y error de CCSD.
- Calcular M exacto con las varianzas reales para compararlo con la cota.
- Estos son los puntos **medidos** del ajuste.
- Salida: `resultados/vqe_hn.csv`.

### Fase 3: VM en GCP · ~medio día, menos de US$30
- **3a. FCI y CCSD de H₁₀ a H₁₆** (d_eq y d_est) con PySCF:
  - Fijar `max_memory`.
  - Guardar el resultado de cada N apenas termina (checkpoint).
  - Si RHF o CCSD no convergen en la zona estirada, **se reporta como resultado**: es justamente la falla del método clásico.
  - Salida: `resultados/fci_ccsd.csv`.
- **3b. Benchmark de la simulación clásica:**
  - Circuitos de 20, 22, …, 32 qubits, por ejemplo k capas de rotaciones + CZ en línea y/o las primeras k excitaciones UCCSD de H_N.
  - Aer con `method="statevector"`, en precisión doble y simple.
  - Medir tiempo por capa (3 repeticiones) y memoria pico (`psutil`).
  - Medir además el tiempo de calcular el valor esperado de un lote de términos de Pauli (por ejemplo, 1.000) para extrapolar el costo de ⟨Ĥ⟩ completo.
  - Salida: `resultados/bench.csv`.
- **3c (opcional).** Términos y grupos de H₁₄ y H₁₆ si no entraron en Colab. El agrupamiento arma un grafo de conmutación con decenas de miles de nodos (~T² comparaciones). Si es muy pesado, muestrear términos o extrapolar.

### Fase 4: ajustes y extrapolación · ~2 h
- Ajuste **log-log** (ley de potencias) de términos, CZ, parámetros y shots vs. N. Reportar el exponente con su error y compararlo con lo esperado en teoría (~N⁴ términos y parámetros, ~N⁵ CZ con Jordan-Wigner).
- Ajuste **semilog** (exponencial) del tiempo y la memoria de simulación vs. qubits.
- Extrapolar a **N = 54** (Reiher) y **N = 76** (Li & Chan) **con bandas de incertidumbre**.
- Calcular:
  - P(sin errores) para FeMoco con `ruido.py`;
  - el tiempo que tardaría un VQE completo de H₁₆ en la VM;
  - los shots totales;
  - la comparación con la estimación con corrección de errores de Lee et al. 2021.

### Fase 5: sección del notebook · ~2 h
- Título sugerido: **«Extra: ¿y la nitrogenasa?»**.
- Una variable `RECALCULAR_PESADO = False`. Por defecto, la sección carga los CSV precalculados para que el notebook siga corriendo de punta a punta en Colab **sin la VM**.
- Recalcular en vivo los casos chicos (N ≤ 8) para mostrar que es reproducible.
- Para lo pesado, documentar el script, la máquina usada, el tiempo y el costo.
- Gráficos:
  1. Recursos cuánticos vs. N en log-log (términos, CZ, parámetros, shots) con los ajustes y líneas verticales en N = 54 y 76.
  2. Costo de la simulación clásica vs. qubits en semilog (memoria, tiempo por capa, dimensión FCI) extrapolado a 108 qubits.
  3. Error de CCSD vs. N en la zona estirada.
  4. P(sin errores) con `ruido.py` vs. N.
- Una tabla final que compare H₄, H₁₆ y FeMoco (extrapolado), más la estimación con corrección de errores.

**Hipótesis a confirmar o refutar con los datos** (no escribirla como conclusión antes de tenerlos): VQE escala de forma polinomial en qubits y compuertas, pero la cantidad de mediciones y la acumulación de ruido lo vuelven impracticable para FeMoco. Simularlo clásicamente escala de forma exponencial. Por eso las estimaciones serias para FeMoco apuntan a QPE con corrección de errores.

## 7. Plan de GCP

- **Créditos:** ~US$350. Primero hay que verificar de qué tipo son. Con el *Free Trial* no se pueden usar GPUs ni pedir aumentos de cuota. Para este plan alcanza con CPU y mucha memoria.
- **Máquina:** `n2-highmem-16` (16 vCPU, 128 GB). Cuesta ~US$1,05/h on-demand o ~US$0,37–0,45/h como Spot en us-central1, según calculadoras de terceros; verificar en la consola.
  - Plan B: `n2-highmem-8` (64 GB). Ahí 32 qubits solo entran en precisión simple (32 GiB).
- **Uso estimado:** 10–20 h de VM, menos de US$30.
- **Antes de crear la VM:** configurar alertas de presupuesto en US$50 y US$100.
- **Spot puede interrumpirse:** los scripts tienen que guardar el resultado de cada N al terminarlo y poder relanzarse salteando lo ya hecho. Correr dentro de `tmux`.
- **Borrar la VM al terminar.** Una `n2-highmem-16` on-demand prendida dos semanas ≈ US$350, o sea, todos los créditos.

```bash
# Crear la VM (Spot)
gcloud compute instances create escalado-h16 \
  --zone=us-central1-a --machine-type=n2-highmem-16 \
  --provisioning-model=SPOT --instance-termination-action=STOP \
  --image-family=debian-12 --image-project=debian-cloud \
  --boot-disk-size=100GB

# En la VM
sudo apt-get update && sudo apt-get install -y python3-venv git tmux
python3 -m venv ~/venv && source ~/venv/bin/activate
pip install qiskit qiskit-aer pyscf numpy scipy pandas psutil

# Bajar resultados a la máquina local
gcloud compute scp --recurse escalado-h16:~/quantum-toffolis/resultados ./resultados --zone=us-central1-a

# Al terminar: borrar
gcloud compute instances delete escalado-h16 --zone=us-central1-a
```

## 8. Riesgos y limitaciones (para escribir en el notebook)

- **H_N no es FeMoco.** FeMoco tiene metales de transición, orbitales d y muchos estados de espín casi degenerados, y su estado fundamental no es un singlete. La extrapolación vale para **contar recursos** (qubits, términos, compuertas), **no** para predecir precisión ni qué tan difícil es optimizar (*barren plateaus*, cantidad de iteraciones).
- **Pocos puntos:** el exponente del ajuste tiene incertidumbre. Reportar bandas, no un número único.
- **Shots:** para N grande son una cota del peor caso. Con mejores técnicas de medición (por ejemplo, medición basada en factorizaciones) podrían ser menos, y hay que decirlo.
- **CZ:** dependen del ansatz. UCCSD es el caso pesimista. Si en el desafío central usamos ADAPT o rotaciones de Givens, mostrar ambos.
- **Costo por evaluación:** un VQE completo lo multiplica por la cantidad de evaluaciones del optimizador.
- **Zona estirada:** RHF y CCSD pueden no converger para N grande. Eso se reporta, no se esconde.
- **Agrupamiento de H₁₆:** puede ser muy costoso de calcular. Si se muestrea o se extrapola, aclararlo.

## 9. Definición de terminado

- [ ] `resultados/recursos.csv` con N = 4…16
- [ ] `resultados/vqe_hn.csv` con H₆ y H₈ (H₁₀ opcional)
- [ ] `resultados/fci_ccsd.csv` con N = 4…16 (o hasta donde haya entrado, aclarado)
- [ ] `resultados/bench.csv` con 20…32 qubits
- [ ] Los 4 gráficos y la tabla de extrapolación a N = 54 y 76, comparada con Lee et al. 2021
- [ ] La sección del notebook corre en Colab sin la VM (carga los CSV) e incluye el párrafo de limitaciones
- [ ] VM borrada y costo gastado anotado en el notebook

## 10. Referencias

- M. Reiher et al., «Elucidating reaction mechanisms on quantum computers», *PNAS* 114, 7555 (2017). [arXiv:1605.03590](https://arxiv.org/abs/1605.03590)
- Z. Li y G. K.-L. Chan et al., «The electronic complexity of the ground-state of the FeMo cofactor of nitrogenase as relevant to quantum simulations», *J. Chem. Phys.* 150, 024302 (2019). [arXiv:1809.10307](https://arxiv.org/abs/1809.10307)
- J. M. Montgomery y D. A. Mazziotti, «Strong electron correlation in nitrogenase cofactor, FeMoco». [arXiv:1805.08746](https://arxiv.org/abs/1805.08746)
- J. Lee et al., «Even more efficient quantum computations of chemistry through tensor hypercontraction», *PRX Quantum* 2, 030305 (2021). [doi:10.1103/PRXQuantum.2.030305](https://doi.org/10.1103/PRXQuantum.2.030305)
- J. F. Gonthier et al., «Measurements as a roadblock to near-term practical quantum advantage in chemistry: resource analysis», *Phys. Rev. Research* 4, 033154 (2022). [arXiv:2012.04001](https://arxiv.org/abs/2012.04001)
- Google Cloud, [Free Trial](https://docs.cloud.google.com/free) y [precios de Spot VMs](https://cloud.google.com/spot-vms/pricing).
