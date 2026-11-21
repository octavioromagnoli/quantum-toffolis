# El diferencial: la misma exactitud con 100× menos CZ, y sobrevive al ruido

> **Estado:** implementado y probado en la rama `octavio/hamiltoniano`: `h4vqe/reduction.py`, `h4vqe/noise.py`, `tests/test_reduction.py` y `scripts/diferencial_demo.py`. Todos los números de este documento salen de `python scripts/diferencial_demo.py` (guarda `docs/ideas/diferencial_resultados.json`). Las distancias se miden en Å, las energías en Ha y los errores en mHa contra FCI.

## La idea en tres líneas

1. Cualquier equipo con IA llega a FCI sin ruido con UCCSD en 8 qubits. Nuestra base ya lo hace, pero con **1648 CZ en línea**. Con nuestro modelo de ruido (sección 2) ese circuito da **561 mHa de error, peor que Hartree-Fock**.
2. Usando las simetrías de la cadena (paridad de N_α, paridad de N_β y la **inversión**), el problema baja a **5 qubits**. Ahí un ansatz en ladrillo nativo de la línea llega a FCI (< 0,001 mHa en las 18 distancias) con **16 CZ y sin SWAP**.
3. Con 16 CZ el ruido deja un sesgo de solo ~14 mHa. La mitigación bien hecha (**plegado global**, no solo de CZ) lo baja a **0,01–0,02 mHa**. Medido con shots (42M por punto) da **−1,1 ± 0,9, +0,4 ± 0,9 y −0,8 ± 0,8 mHa** a 1,8, 2,4 y 3,0 Å: **precisión química con ruido**, y lejos de CCSD (−9 a −31 mHa).

La consigna lo dice en la motivación: *"no alcanza con llegar al número correcto: importa cómo se llega, con cuántas compuertas y con cuántas mediciones"*. Ahí está el diferencial.

## 1. Compuertas y precisión (a 2,0 Å)

| ansatz | qubits | parámetros | CZ (todos con todos) | CZ (línea) | error sin ruido | error con ruido |
|---|---|---|---|---|---|---|
| UCCSD ×1 | 8 | 26 | 292 | 806 | 4,46 mHa | 340 mHa |
| UCCSD ×2 (base actual) | 8 | 52 | 582 | **1648** | 0,000 | **561 mHa** |
| compacto L=6 | 5 | 35 | 12 | 12 | 3,29 | 12,6 |
| compacto L=7 | 5 | 40 | 14 | 14 | 0,000 | 12,8 |
| **compacto L=8 (propuesto)** | **5** | **45** | **16** | **16** | **0,000** | **13,5** |

Como referencia, a 2,0 Å HF tiene 322 mHa de error y CCSD −18,3 mHa. La columna "con ruido" es el valor de expectación exacto con nuestro modelo de ruido y *sin* mitigación, así que mide solo el sesgo.

- La curva completa con el compacto L=8 da error < 0,001 mHa en las 18 distancias, con ⟨S²⟩ = 0 y N = 4. Los barridos de ida (0,6 → 3,0) y de vuelta coinciden. Con L=10 también da < 0,0001 mHa. **Con L=6 falla** (hasta 3,8 mHa, con mínimos locales): eso es un buen análisis de límites para el notebook.
- **Trampa que conviene contar en "qué no funcionó":** si en cada capa se ponen *todas* las CZ a la vez (0-1, 1-2, 2-3, 3-4), el ansatz tiene un **techo de fidelidad de 0,961** que no mejora con más capas (probado con L = 2…8). La fidelidad 1 se alcanza con el patrón en ladrillo (pares alternados) o con una cadena de CX.
- UCCSD ×1 *no* llega a precisión química en la zona estirada (3,3 mHa a 1,8 Å y 5,2 mHa a 2,4 Å). Por eso la base pasó a ×2 y duplicó las CZ.

### Cómo se llega a 5 qubits (`h4vqe/reduction.py`)

- Simetrías Z₂ de tipo Z: P_α = Z₀Z₁Z₂Z₃, P_β = Z₄Z₅Z₆Z₇ y P_u = ΠZ sobre los orbitales impares ante la inversión (σu, detectados automáticamente con los coeficientes de los orbitales moleculares; dan [+, −, +, −] en las 18 distancias).
- Tapering de Bravyi *et al.* (2017): los generadores se llevan a forma escalonada sobre GF(2) para que cada uno tenga un qubit exclusivo; U = (X_q + g)/√2; X_q se reemplaza por el autovalor del sector HF. Se eliminan los qubits 0, 4 y 7.
- Validado: E₀ (5 qubits) = FCI con error de 10⁻¹⁵ en las 18 distancias; ⟨HF|H₅|HF⟩ = E_HF. El estado HF sigue siendo un bitstring.
- Esto va un paso más allá del extra "6 qubits": el mapeo de paridad elimina 2 qubits y la inversión, uno más.

## 2. Ruido y mitigación (`h4vqe/noise.py`)

La consigna ya no incluye `ruido.py`: los organizadores aclararon que quedó de una versión anterior y que no entra. Igual usamos un modelo de ruido **propio**, simple y controlado, con los parámetros que describía esa versión (despolarizante: 0,1 % por cz y 0,005 % por compuerta de 1 qubit; base cz/sx/x/rz; qubits en línea). Sirve para medir cuánto daño hace cada compuerta y para probar la mitigación. Para ruido de hardware real queda el extra con `FakeFez`.

Sesgo de cada técnica, como valor de expectación exacto con ruido. Esto sirve para *elegir* la técnica; no son energías reportadas.

| técnica | 1,8 Å | 3,0 Å | comentario |
|---|---|---|---|
| sin mitigar | +15,5 | +13,4 | 16 CZ ≈ 1,6 % + ~150 compuertas de 1 qubit ≈ 0,7 % |
| ZNE **plegando solo las CZ** (λ = 1, 3, 5) | +4,0 | +3,5 | ❌ las compuertas de 1 qubit no se escalan |
| ZNE global lineal (λ = 1, 3) | +0,42 | +0,36 | |
| **ZNE global Richardson (λ = 1, 3, 5)** | **+0,01** | **+0,01** | ✅ pero amplifica σ ×3,5 (12× shots) |
| eco C·C† solo (reescalado depolarizante) | −1,53 | −1,52 | ❌ sobrecorrige: queda **por debajo de FCI** |
| **eco + ZNE lineal (λ = 1, 3)** | **+0,02** | **+0,02** | ✅ amplifica σ ×2 (4× shots): **el recomendado** |

- **Plegado global:** C → C (C† C)^k. Así se escalan todas las compuertas. Es clave verificar que C† tenga las mismas compuertas que C. Transpilado con `optimization_level=0`, C† salía con **368 rz contra 72**; el factor de ruido dejaba de ser λ y la extrapolación se pasaba unos 4 mHa. `noise._inverse()` ahora lo chequea con un `assert`.
- **Eco:** C·C†|0⟩ debería volver a |00000⟩. P(0) da la fidelidad F, y se reescala E′ = (E − (1 − F^λ)·Tr(H)/2ⁿ)/F^λ. No usa ningún dato clásico. Solo no alcanza, porque el ruido es local y no global, pero combinado con ZNE lineal sí.
- Para repartir los shots entre niveles: N_λ ∝ |γ_λ|, o sea 3:1 con λ = (1, 3) (ver `noise.zne_shot_split`).
- **Resultado con shots** (eco + ZNE lineal: 30M de shots a λ=1, 10M a λ=3 y 2M en el eco; 42M por punto):

  | d | sin mitigar (10M shots) | eco + ZNE (1, 3) | CCSD − FCI |
  |---|---|---|---|
  | 1,8 | +14,8 ± 0,9 | **−1,1 ± 0,9** | −9,1 |
  | 2,4 | +11,5 ± 0,9 | **+0,4 ± 0,9** | −30,6 |
  | 3,0 | +13,0 ± 0,8 | **−0,8 ± 0,8** | −23,0 |

  La fidelidad del eco da F = 0,981 por pasada. Ya sin mitigar, el circuito compacto le gana a CCSD desde 2,4 Å.

## 3. Curva medida sin ruido y shots

**Curva** (`SamplerV2` + QWC + Neyman, `estimator.py` del equipo sin cambios): 8M de shots por punto con d < 1,8 Å (σ ≈ 0,7–1,0 mHa) y 24M en la zona estirada (σ ≈ 0,51–0,57 mHa).
- Zona estirada: **le gana a CCSD contando la σ en 7/7 puntos y tiene |E − FCI| < 1,6 mHa en 7/7**. El peor caso de |E − FCI| + σ es 2,1 mHa, a 1,8 Å.
- ⟨S²⟩ medido = 0,00 ± 0,01 y N = 4,000 en los 18 puntos.
- Calibración: el circuito es exacto (|E_circuito − FCI| < 10⁻⁸ mHa), así que el desvío de cada punto es estadístico. A 0,6 Å el valor medido salió **3,7σ por debajo de FCI**. Repetí la medición y el estimador está bien calibrado: std(z) = 0,93 con 40 repeticiones de 500k shots y 0,98 con 16 repeticiones de 8M shots. Fue una fluctuación poco probable. **Se reporta como salió**, con la calibración al lado, sin volver a correr hasta que dé bien.

### El costo real y el intercambio

- La reducción a 5 qubits baja las CZ 100×, pero **sube unas 2× los shots** para la misma σ con QWC + Neyman. El costo (Σ_g √Var_g)² con el estado exacto da 7,7 Ha² en 5 qubits contra 3,9 en 8 qubits a 2,0 Å. Ninguna de las elecciones de qubits a eliminar lo arregla: la mejor ahorra solo 7 %. Lo declaramos: **sin ruido, los shots solo cuestan tiempo de simulación; con ruido, el circuito de 8 qubits no sirve con ningún número de shots.**
- Cifras útiles: σ = 0,5 mHa necesita ~3×10⁷ shots por punto en 5 qubits. La demo usa 8M de shots (σ ≈ 1 mHa) cerca del equilibrio y 24M (σ ≈ 0,55 mHa) en la zona estirada, que es donde hace falta margen para la precisión química.
- **Punto para el pitch:** `EstimatorV2` de qiskit_aer no mide. Calcula el valor exacto y le suma ruido gaussiano N(0, *precision*), y devuelve `stds = precision` como constante (ver `qiskit_aer/primitives/estimator_v2.py`, líneas 140–151). Los equipos que lo usen van a reportar "shots" y σ que son un parámetro elegido a mano. Nosotros medimos de verdad con `SamplerV2`, con σ sacada de las muestras.
- Mi calibración independiente (30 repeticiones a 2,0 Å): la σ reportada coincide con la dispersión real, con 70 % de cobertura en 1σ (esperado: 68 %) y 100 % en 2σ.

## 4. Confiabilidad sin mirar FCI (para `06_reliability.md`)

1. **A 3,0 Å la energía sola no prueba nada.** En el sector de simetría, el triplete está a 0,95 mHa del fundamental, el segundo singlete a 1,20 y el quintete a 1,65: los 16 estados de espín de 4 átomos de H caen dentro de 1,65 mHa. Hasta 4 átomos aislados, sin entrelazar, quedan a 0,96 mHa de FCI. Por eso **medimos ⟨S²⟩ y N con shots** (S² = 0,00 ± 0,01 y N = 4,000 en toda la curva), y por eso la búsqueda de parámetros penaliza μ⟨S²⟩ (μ = 0,5). **Lo vimos pasar:** sin la penalización, el barrido de vuelta a 3,0 Å convergió al **triplete (⟨S²⟩ = 1,99) con error de 0,95 mHa**. Eso es "precisión química" con el estado equivocado, y solo se detecta midiendo S².
2. **Límite de átomos separados:** E_H(STO-3G) = −0,4665818 Ha se calcula exacto (un solo electrón), así que 4E_H = −1,866327 Ha. A 3,0 Å, **CCSD queda 24 mHa por debajo de 4 átomos aislados**. Eso es no físico y se ve sin calcular FCI. El VQE medido queda 1,2 ± 0,5 mHa por debajo (FCI: 0,96), una interacción residual chica y física.
3. Variacional: todas nuestras energías están por encima de las de otros métodos variacionales (HF, UCCSD ×1), y entre variacionales gana el más bajo. CCSD no es variacional.
4. Robustez de la búsqueda: el barrido de ida y el de vuelta coinciden y los reinicios aleatorios llegan a la misma energía. La convergencia con la profundidad da L=6 ≠ L=7 = L=8 = L=10.
5. Por qué falla CCSD: las ocupaciones naturales a 3,0 Å son [1,09, 1,05, 0,95, 0,91]. Son cuatro electrones prácticamente desapareados: correlación fuerte y multirreferencia, el régimen donde una referencia única con excitaciones dobles no alcanza.

## 5. Qué hace cada uno (compatible con `PLAN.md`)

- **Octavio:** este PR (`reduction.py`, `noise.py`, test y demo). Su sección de `02_hamiltonian.md` lleva las simetrías y la reducción a 5 qubits; los extras, el ruido y "menos qubits". Siguiente paso: FakeFez a 2,4 Å con el compacto (16 CZ entran en una línea de 5 qubits de Fez).
- **Gonzalo:** sumar el compacto L = 6, 7, 8 a la tabla de "comparar ansätze". ADAPT sigue siendo una buena comparación físicamente motivada. **Propuesta para el sync de las 2:00: el compacto L=8 como ansatz principal y UCCSD ×2 como referencia de 8 qubits.**
- **Santino:** ZNE con **plegado global** (`noise.fold_global`), no solo de CZ; eco + ZNE lineal con reparto 3:1 (`noise.zne`). Probar agrupamiento por conmutación general en 5 qubits para recuperar el 2× de shots. `noise.measure_energy_noisy` reutiliza las funciones de `estimator.py` sin tocarlo.
- **Joaquín:** integrar el compacto al pipeline y agregar `reduction` y `noise` a `MODULES` en `build_notebook.py`. Este fragmento en `run_point` alcanza:

  ```python
  from .reduction import tapering, build_compact_ansatz, optimize_compact, reduced_number_operators
  if config.ANSATZ == "compact":               # nuevo en config.py, junto con COMPACT_LAYERS = 8
      H = tapering(m, H).apply(H)              # 5 qubits
      ans = build_compact_ansatz(m, layers=config.COMPACT_LAYERS)
      opt = optimize_compact(ans, H, reduced_number_operators(m)["S2"], x0=x0)
  ```

  `measure_energy`, `cz_count` y la tabla funcionan igual: el `Ansatz` compacto tiene la misma interfaz.

## 6. Limitaciones (para declararlas, no para esconderlas)

- El ansatz en ladrillo es *hardware-efficient*: no conserva N ni S² por construcción (lo verificamos midiéndolos) y con L ≤ 6 aparecen mínimos locales. En 5 qubits no hay *barren plateaus*, pero esto no escala a moléculas grandes. El tapering, en cambio, sí escala.
- El modelo de ruido es propio (despolarizante, no el de un hardware real). Si las rz no tienen ruido (en IBM son virtuales), el sesgo sin mitigar baja de ~15 a ~13 mHa a 1,8 Å. Las conclusiones no cambian.
- Unas 2× más shots que en 8 qubits para la misma σ (sección 3).
