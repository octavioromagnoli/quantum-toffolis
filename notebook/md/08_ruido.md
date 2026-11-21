## 7 · Problema 5, el ruido: las CZ deciden
Desde acá empiezan los **extras** de la consigna. El desafío central ya está cumplido (§5–§6). La pregunta ahora es si el método sobrevive al ruido de una computadora real.

### 7.1 El modelo de ruido

La consigna describe un modelo con compuertas `cz`, `sx`, `x` y `rz`, qubits conectados **en línea** y error de 0,1 % por CZ y 0,005 % por compuerta de un qubit. El archivo `ruido.py` que menciona no se distribuye (los organizadores aclararon que quedó de una versión anterior), así que implementamos el **mismo modelo como propio** (`h4vqe/noise.py`): canales despolarizantes después de cada compuerta,

$$\mathcal E_p(\rho) = (1 - p)\,\rho + p\,\frac{I}{2^n}, \qquad p_{\rm CZ} = 10^{-3}\ (n = 2), \qquad p_{1q} = 5\cdot 10^{-5}\ (n = 1).$$

Con $n_{\rm CZ}$ compuertas de dos qubits y $n_{1q}$ de uno, la probabilidad de que el circuito corra sin ningún error es

$$P_{\rm limpio} = (1 - p_{\rm CZ})^{n_{\rm CZ}}\,(1 - p_{1q})^{n_{1q}} \approx e^{-n_{\rm CZ}/1000}.$$

**Observación.** Con este modelo, el VQE "de libro" no sirve: UCCSD ×2 en línea tiene ~1550 CZ, $P_{\rm limpio} \approx$ 0,2, y el estado se mezcla con ruido. ADAPT baja a 313 CZ con el ruteo fSWAP, pero sigue siendo demasiado.

**Hipótesis.** El sesgo escala con la cantidad de compuertas. Para bajarlo a algo que la mitigación pueda corregir (unos pocos mHa) hacen falta **unas 100 veces menos CZ**. Un ansatz químico en 8 qubits no llega ahí; uno nativo de la línea, en menos qubits, sí.

La celda de abajo lo verifica a 2,4 Å con tres ansätze que sin ruido llegan a ≤ 0,4 mHa de FCI. Calcula el valor esperado **exacto** con ruido (probabilidades con matriz densidad, sin shots), o sea, el sesgo que el ruido deja y que más shots no sacan. Los circuitos se transpilan a la línea y $\hat H$ se reordena con el layout final.
