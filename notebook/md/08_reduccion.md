**Qué dio.** Sin mitigar, UCCSD ×2 (lo mínimo) queda **peor que no hacer nada**: su error con ruido supera al de HF, que no necesita ningún circuito. ADAPT mejora tres veces, pero sigue a cientos de mHa. Solo el compacto de 16 CZ deja un sesgo de pocos mHa, que la mitigación (§8) puede sacar. Las dos subsecciones que siguen explican cómo se llega a ese circuito.

### 7.2 Menos qubits: tapering con simetrías ($8 \to 6 \to 5$)

$\hat H$ conmuta con tres operadores de tipo $Z$ (`h4vqe/reduction.py`):

$$P_\alpha = Z_0Z_1Z_2Z_3 = (-1)^{N_\alpha}, \qquad P_\beta = Z_4Z_5Z_6Z_7 = (-1)^{N_\beta}, \qquad P_u = \prod_{p\,\in\,u} Z_{p}\,Z_{p+4},$$

donde $u$ son los orbitales moleculares **impares** ante la inversión de la cadena (átomo $i \leftrightarrow$ átomo $3-i$). Cada orbital es par o impar, y el producto de las paridades de los ocupados se conserva. En el estado que buscamos, las tres tienen un valor fijo (el del determinante HF).

El *tapering* (Bravyi, Gambetta, Mezzacapo y Temme, 2017) usa cada simetría $g_i$ para eliminar un qubit $q_i$. Primero se llevan los generadores a forma escalonada en $\mathrm{GF}(2)$, para que cada uno tenga un qubit "exclusivo". Con el Clifford

$$U_i = \frac{X_{q_i} + g_i}{\sqrt 2}, \qquad U_i\, g_i\, U_i = X_{q_i},$$

en $U\hat HU$ el qubit $q_i$ aparece solo como $I$ o $X_{q_i}$, y $X_{q_i}$ se reemplaza por el autovalor $\pm 1$ de $g_i$ en el sector HF. El hamiltoniano queda con un qubit menos por simetría, y su espectro en ese sector es el mismo.

- **6 qubits** (solo $P_\alpha$ y $P_\beta$): equivale al mapeo de paridad con la reducción de 2 qubits, que es lo que pide el extra. Toda excitación conserva $N_\alpha$ y $N_\beta$, así que **los 26 operadores del pool de UCCSD/ADAPT se reducen** y el ansatz químico funciona igual. No depende de los orbitales.
- **5 qubits** (además $P_u$): elimina un qubit más, pero exige orbitales canónicos (con paridad definida), y solo se reducen las excitaciones simétricas ante la inversión (14 de 26).

### 7.3 Un ansatz nativo de la línea: 16 CZ

En 5 qubits usamos un ansatz **en ladrillo** (*brickwork*): $X$ sobre el HF reducido, una capa de $R_Y$ y $L$ capas de [CZ entre vecinos en pares alternados + $R_Y$ en cada qubit]:

$$\lvert\psi(\boldsymbol\theta)\rangle = \prod_{l=1}^{L}\Big[\textstyle\bigotimes_q R_Y(\theta_{l,q})\ \prod_{(q,q+1)\,\in\,\mathcal B_l}\mathrm{CZ}_{q,q+1}\Big]\ \bigotimes_q R_Y(\theta_{0,q})\,\lvert \mathrm{HF}_5\rangle ,$$

donde $\mathcal B_l$ alterna entre los pares $(0,1),(2,3)$ y $(1,2),(3,4)$. Con $L = 8$ tiene **16 CZ, ninguna SWAP con qubits en línea** y 45 parámetros. Como no conserva $N$ ni $S^2$ por construcción, la búsqueda de $\boldsymbol\theta^*$ minimiza $E(\boldsymbol\theta) + \mu\,\langle S^2\rangle$ con $\mu$ = 0,5 (la energía que se reporta es solo $E$). **Sin esa penalización, una corrida convergió al triplete con 0,95 mHa de error**: "precisión química" con el estado equivocado. Por eso también hay que medir $S^2$.

La celda de abajo valida el tapering en las 18 distancias, muestra cuánto cambia la medición en 8, 6 y 5 qubits y compara el compacto con UCCSD ×2 sin ruido a 2,0 Å.
