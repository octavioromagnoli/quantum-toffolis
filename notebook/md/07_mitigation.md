### Extra: mitigación de errores  *(Santino)*

**Planteo.** Usamos el ansatz compacto de 5 qubits con 16 CZ (sección 1) en el modelo de ruido propio: despolarizante, 0,1 % por CZ y 0,005 % por compuerta de 1 qubit. Le agregamos un **error de lectura de 2 % por qubit**, porque el modelo base no trae y el hardware real sí. Probamos tres ingredientes y sus combinaciones (`h4vqe/mitigation.py`):

1. **Corrección de lectura (RO).** Calibramos la matriz de confusión de cada qubit midiendo todo en |0⟩ y todo en |1⟩, con 20 M de shots por estado. Después corregimos los *valores* de cada grupo con (A⁻¹)ᵀ, sin tocar los histogramas: así la σ sigue saliendo de las muestras.
2. **Post-selección por simetría (PS).** N_α y N_β siguen siendo diagonales (solo Z) después de la reducción a 5 qubits. En un grupo que mide en Z todos los qubits del soporte de N_α, cada shot da su N_α, y si N_α ≠ 2 ese shot tiene un error seguro y se descarta. Pasa lo mismo con N_β. Se puede chequear en 5 de los 35 grupos, que son los de más peso (~45 % de Σ|c|).
   Para combinarla con RO usamos un **estimador de cociente** por grupo: E_g = media(u)/media(w), con u = keep·v y w = keep, las dos corregidas por (A⁻¹)ᵀ. La σ sale por el método delta. Sin mitigación, u = v y w = 1, y queda el estimador de siempre.
3. **Extrapolación a ruido cero (ZNE)** con plegado global C → C(C†C)^k, que escala el ruido por λ = 2k+1 (`noise.fold_global`, de Octavio). También probamos el **reescalado por la fidelidad del eco** C·C†|0⟩ (idea de Octavio, ahora con la lectura corregida en el eco).

**Hallazgo: hay que plegar también la rotación de medición.** Cada grupo agrega 2–3 compuertas de 1 qubit por qubit para rotar a su base (R). Si se pliega solo el ansatz C y se agrega R después, el ruido de R no escala con λ y queda un sesgo que la extrapolación no saca. Plegando R·C entero, ese sesgo desaparece:

| plegado (sesgo exacto, 1,8 Å) | ZNE lineal (1,3) | eco + ZNE lineal | Richardson (1,3,5) |
|---|---|---|---|
| solo C, R después | +0,53 mHa | +0,19 mHa | +0,16 mHa |
| **R·C entero** | **+0,38 mHa** | **+0,030 mHa** | **+0,002 mHa** |

**Qué aporta cada técnica (sesgo exacto, sin shots).** Calculamos las probabilidades exactas con matriz densidad, para elegir la técnica antes de gastar shots (`scripts/mitigation_study.py`). "Shots" es lo necesario para σ = 0,5 mHa, con el reparto óptimo entre niveles de λ.

| técnica (con 2 % de error de lectura) | sesgo 1,8 Å | sesgo 2,4 Å | sesgo 3,0 Å | shots (2,4 Å) |
|---|---|---|---|---|
| sin mitigar | +110,5 | +109,1 | +110,8 | 10 M |
| RO | +14,7 | +13,8 | +13,3 | 14 M |
| PS + RO | +10,9 | +11,0 | +11,4 | 13 M |
| ZNE lineal (1,3) + PS + RO | +0,27 | +0,29 | +0,30 | 54 M |
| **eco + ZNE lineal (1,3) + RO** | **+0,030** | **+0,008** | **+0,006** | **61 M** |
| ZNE Richardson (1,3,5) + PS + RO | −0,001 | −0,003 | −0,004 | 164 M |

(Valores en mHa respecto de la energía sin ruido del mismo circuito, que coincide con FCI dentro de 10⁻¹¹ mHa. Los θ* quedan guardados en `results/compact_params/` junto con el H₅ con el que se optimizaron, así que el estudio exacto y la corrida con shots usan exactamente el mismo circuito y el mismo H. Con otro θ* de la misma energía, el sesgo sin mitigar cambia ~0,5 mHa.)

**Con shots** (`scripts/mitigation_shots.py`; todas las energías medidas con `SamplerV2`, error = E − FCI):

| d (Å) | sin mitigar (10 M) | RO (14 M) | PS + RO (14 M) | ZNE lineal + PS + RO (54 M) | **eco + ZNE lineal + RO (60 M)** | CCSD − FCI |
|---|---|---|---|---|---|---|
| 1,8 | +110,4 ± 0,6 | +14,2 ± 0,6 | +10,3 ± 0,5 | −0,5 ± 0,5 | **+0,3 ± 0,6** | −9,1 |
| 2,4 | +109,0 ± 0,5 | +14,4 ± 0,5 | +11,3 ± 0,5 | −0,1 ± 0,5 | **−1,0 ± 0,5** | −30,6 |
| 3,0 | +110,5 ± 0,5 | +14,2 ± 0,5 | +10,6 ± 0,5 | +0,3 ± 0,5 | **+0,1 ± 0,5** | −23,0 |

(E − FCI en mHa ± σ; entre paréntesis, los shots por punto, incluidos los del eco. La matriz de lectura se calibra con shots aparte: 20 M por estado.)

- **Con mitigación, las 6 energías quedan dentro de la precisión química** (|E − FCI| < 1,6 mHa) y dentro de 2σ de FCI. Las dos técnicas recomendadas le ganan a CCSD por 9 a 31 mHa.
- **Las 15 mediciones coinciden con el sesgo exacto** que predecía el estudio sin shots: χ² = 15,4 con 15 grados de libertad (p = 0,42) y ningún desvío pasa de 1,9σ. La σ con ruido y mitigación también está bien calibrada.
- **La calibración de la lectura deja un error sistemático de −0,15, −0,05 y +0,00 mHa** a 1,8, 2,4 y 3,0 Å, comparando la matriz calibrada contra la exacta. No entra en σ, pero es menor que σ/3.

**Conclusiones**

- **La corrección de lectura es imprescindible y barata.** Con 2 % de error por qubit, la lectura agrega ~95 mHa, mucho más que las 16 CZ. Plegar no amplifica el error de lectura, así que la ZNE sola no lo puede sacar. Cuesta ×1,4 en shots.
- **La post-selección es gratis pero rinde poco:** baja el sesgo un 15–26 % sin costo en shots, porque se descarta menos del 1 % de los shots. Rinde poco porque solo se puede chequear en 5 de 35 grupos, y porque los errores de fase (Z) no cambian N y no se detectan.
- **Recomendado: eco + ZNE lineal + RO.** Deja ≤ 0,03 mHa de sesgo y cuesta ~4× los shots de una medición con solo RO (~60 M de shots por punto para σ = 0,5 mHa). Richardson (1,3,5) llega a 0,004 mHa, pero cuesta 12×.

**Qué no funcionó**

- Plegar solo el ansatz (lo que hace `noise.zne`): deja 0,16–0,23 mHa de sesgo con Richardson o con eco + lineal, escondido debajo de σ ≈ 1 mHa.
- Combinar el eco con la post-selección: el reescalado supone que el ruido lleva al estado totalmente mezclado, y post-seleccionar cambia esa referencia.
- Guardar solo θ* para reusarlo en otro proceso: el signo de los orbitales de PySCF cambia según cómo se reparta el cálculo entre hilos, y con otro signo cambian términos de H. Un θ* reusado así dio +668 mHa en vez de 0. Por eso guardamos θ* junto con su H.

**Limitaciones**

- El modelo de ruido es propio, y la corrección de lectura supone errores independientes por qubit. El error de calibración de la lectura es sistemático y no entra en σ: lo medimos comparando contra la matriz exacta y da hasta 0,15 mHa con 20 M de shots por estado.
- La ZNE supone que el ruido escala con el plegado, lo que es exacto en nuestro modelo y aproximado en hardware real.
- La post-selección solo sirve en los grupos que miden en Z el soporte de N_α o N_β.
