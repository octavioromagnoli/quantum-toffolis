**Sesgo exacto** (tabla de arriba, 2,4 Å; los valores cambian unas décimas de mHa con el $\boldsymbol\theta^*$ que encuentre el optimizador): sin mitigar, ~+109 mHa, que es casi todo error de lectura; con RO, ~+14 mHa, el ruido de las 16 CZ y de las compuertas de un qubit; con PS + RO, ~+11 mHa; con ZNE lineal + PS + RO, ~+0,3 mHa; **con eco + ZNE lineal + RO, ~0,01 mHa**; y con Richardson, ~0,003 mHa, pero a 12 veces el costo en shots.

**Con shots** (`scripts/mitigation_shots.py`, ~20 min por punto, por eso no se recalcula acá; todas las energías medidas con `SamplerV2`, error = $E - E_{\rm FCI}$ en mHa ± σ; entre paréntesis, los shots por punto, incluidos los del eco):

| d (Å) | sin mitigar (10 M) | RO (14 M) | PS + RO (14 M) | ZNE lineal + PS + RO (54 M) | **eco + ZNE lineal + RO (60 M)** | CCSD − FCI |
|---|---|---|---|---|---|---|
| 1,8 | +110,4 ± 0,6 | +14,2 ± 0,6 | +10,3 ± 0,5 | −0,5 ± 0,5 | **+0,3 ± 0,6** | −9,1 |
| 2,4 | +109,0 ± 0,5 | +14,4 ± 0,5 | +11,3 ± 0,5 | −0,1 ± 0,5 | **−1,0 ± 0,5** | −30,6 |
| 3,0 | +110,5 ± 0,5 | +14,2 ± 0,5 | +10,6 ± 0,5 | +0,3 ± 0,5 | **+0,1 ± 0,5** | −23,0 |

- **Con mitigación, las 6 energías quedan dentro de la precisión química** ($\lvert E - E_{\rm FCI}\rvert$ < 1,6 mHa) y a menos de 2σ de FCI, **con ruido**. Las dos técnicas recomendadas le ganan a CCSD por 9 a 31 mHa. Esto cumple el extra "precisión química en toda la zona estirada" también con ruido, con un circuito de 16 CZ.
- **Las 15 mediciones coinciden con el sesgo exacto** que predecía el estudio sin shots: $\chi^2$ = 15,4 con 15 grados de libertad ($p$ = 0,42), y ningún desvío pasa de 1,9σ. La $\sigma$ con ruido y mitigación también está bien calibrada.
- **La calibración de la lectura deja un error sistemático** de −0,15, −0,05 y +0,00 mHa a 1,8, 2,4 y 3,0 Å (comparando la matriz calibrada con la exacta). No entra en la $\sigma$, pero es menor que $\sigma/3$.

**Por qué funcionó**

- **La corrección de lectura es imprescindible y barata.** Con 2 % de error por qubit, la lectura agrega ~95 mHa, mucho más que las 16 CZ. Plegar el circuito no repite la medición, así que la ZNE sola no puede sacar ese error. Cuesta 1,4× en shots.
- **La post-selección es gratis, pero rinde poco:** baja el sesgo un 15–26 % sin costo en shots (se descarta menos del 1 % de los shots). Solo se puede chequear en 5 de 35 grupos, y los errores de fase ($Z$) no cambian $N$, así que no se detectan.
- **La ZNE funciona porque el ruido escala con el plegado** (en nuestro modelo, exactamente), y el eco corrige la curvatura que la extrapolación lineal no ve. **Recomendado: eco + ZNE lineal + RO**: deja ≤ 0,03 mHa de sesgo y cuesta ~4 veces los shots de una medición con solo RO (~60 M por punto para $\sigma$ = 0,5 mHa).

**Qué no funcionó**

- **Plegar solo las CZ:** deja 4 mHa de sesgo, porque las ~150 compuertas de un qubit también suman ruido y no se escalan. El plegado tiene que ser global.
- **Plegar solo el ansatz y no la rotación de medición** (lo que hacía `noise.zne` al principio): deja 0,16–0,23 mHa de sesgo, escondido debajo de $\sigma \approx$ 1 mHa.
- **Un $C^\dagger$ con otras compuertas que $C$:** transpilado con `optimization_level=0`, $C^\dagger$ salía con 368 `rz` contra 72, y el factor de ruido dejaba de ser $\lambda$: la extrapolación se pasaba ~4 mHa. Ahora `noise._inverse` lo chequea.
- **El eco solo, sin ZNE:** sobrecorrige (queda ~1,5 mHa **por debajo** de FCI), porque el ruido es local y no global.
- **Combinar el eco con la post-selección:** el reescalado supone que el ruido lleva al estado totalmente mezclado, y post-seleccionar cambia esa referencia.
- **Guardar solo $\boldsymbol\theta^*$ para usarlo en otro proceso:** el signo de los orbitales de PySCF puede cambiar entre procesos, y con otro signo cambian términos de $\hat H$. Un $\boldsymbol\theta^*$ reusado así dio +668 mHa en vez de 0. Por eso acá $\boldsymbol\theta^*$ y $\hat H_5$ salen siempre del mismo proceso.

**Limitaciones:** el modelo de ruido es propio, y la corrección de lectura supone errores independientes por qubit. La ZNE supone que el ruido escala con el plegado, algo exacto en nuestro modelo y aproximado en hardware real. La post-selección solo sirve en los grupos que miden en $Z$ el soporte de $N_\alpha$ o $N_\beta$.
