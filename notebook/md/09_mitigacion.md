## 8 · Mitigación de errores: precisión química con ruido
**Observación.** Con 16 CZ, el compacto queda a unos +14 mHa de FCI (§7.1). En hardware real se suma el **error de lectura**: cada qubit medido puede devolver el bit equivocado.

**Hipótesis.** Los tres sesgos tienen orígenes distintos y se pueden corregir por separado: la lectura con una calibración, los errores que cambian $N_\alpha$ o $N_\beta$ descartando esos shots, y el ruido de compuertas extrapolando a ruido cero. Combinados, deberían bajar el sesgo por debajo de la $\sigma$ de los shots.

Agregamos al modelo de §7.1 un **error de lectura de 2 % por qubit** (el modelo base no trae y el hardware real sí) y probamos (`h4vqe/mitigation.py`):

**1. Corrección de lectura (RO).** Para cada qubit, la matriz de confusión $A^{(q)}$ con $A^{(q)}_{ij} = P(\text{leer } i \mid \text{estado } j)$ se calibra preparando todo en $\lvert 0\rangle$ y todo en $\lvert 1\rangle$ (20 M de shots por estado). Si los errores son independientes, $\vec p_{\rm leído} = \big(A^{(n-1)}\otimes\cdots\otimes A^{(0)}\big)\,\vec p_{\rm real} \equiv A\,\vec p_{\rm real}$. En lugar de invertir histogramas, corregimos la **tabla de valores** de cada grupo:

$$\langle\hat O_g\rangle = \vec v_g\cdot\vec p_{\rm real} = \big((A^{-1})^{T}\,\vec v_g\big)\cdot\vec p_{\rm leído},$$

así cada shot sigue dando una muestra (con el valor corregido) y la $\sigma$ sigue saliendo de las muestras.

**2. Post-selección por simetría (PS).** $\hat N_\alpha$ y $\hat N_\beta$ siguen siendo diagonales (solo $Z$) después de la reducción a 5 qubits. En un grupo que mide en $Z$ todos los qubits de su soporte, cada shot trae su $N_\alpha$ y $N_\beta$, y si alguno es distinto de 2 el shot tuvo un error seguro y se descarta. Con lectura corregida usamos un **estimador de cociente** por grupo:

$$E_g = \frac{\overline{u}}{\overline{w}}, \qquad u(x) = \mathrm{keep}(x)\,v_g(x),\quad w(x) = \mathrm{keep}(x)\quad(\text{los dos corregidos con } (A^{-1})^T),$$

con la $\sigma$ por el método delta. Se puede chequear en 5 de los 35 grupos, que son los de más peso (~45 % de $\sum\lvert c\rvert$).

**3. Extrapolación a ruido cero (ZNE).** Se **pliega** el circuito, $C \to C\,(C^\dagger C)^k$, que hace la misma operación con $\lambda = 2k+1$ veces más compuertas, y por lo tanto más ruido. Después se extrapola a $\lambda = 0$:

$$E(0) \approx \sum_i \gamma_i\, E(\lambda_i), \qquad \gamma_i = \prod_{j\neq i}\frac{\lambda_j}{\lambda_j - \lambda_i}, \qquad \sigma_0^2 = \sum_i \gamma_i^2\,\sigma_i^2 .$$

Lineal con $\lambda = (1, 3)$: $E(0) = \tfrac32 E(1) - \tfrac12 E(3)$. Richardson con $(1, 3, 5)$: $\gamma = \big(\tfrac{15}{8}, -\tfrac{5}{4}, \tfrac38\big)$. El precio es que los $\gamma_i$ amplifican la $\sigma$.

**4. Eco + ZNE.** El circuito $C\,C^\dagger\lvert 0\rangle$ debería volver a $\lvert 0\cdots 0\rangle$. Su probabilidad de volver mide la fidelidad por pasada $F$ ($P_0 = F^2 + (1 - F^2)/2^n$). Con un ruido despolarizante global, $E_\lambda = F^\lambda E + (1 - F^\lambda)\,E_{\rm mix}$, con $E_{\rm mix} = \mathrm{Tr}(\hat H)/2^n$, así que se puede reescalar antes de extrapolar:

$$E'_\lambda = \frac{E_\lambda - (1 - F^{\lambda})\,E_{\rm mix}}{F^{\lambda}} .$$

**Hallazgo: hay que plegar también la rotación de medición.** Cada grupo agrega 2–3 compuertas de un qubit por qubit para rotar a su base ($R$). Si se pliega solo el ansatz $C$ y se agrega $R$ después, el ruido de $R$ no escala con $\lambda$ y deja un sesgo que la extrapolación no saca. Plegando $R\cdot C$ entero, ese sesgo desaparece (sesgo exacto a 1,8 Å):

| plegado | ZNE lineal (1, 3) | eco + ZNE lineal | Richardson (1, 3, 5) |
|---|---|---|---|
| solo $C$, $R$ después | +0,53 mHa | +0,19 mHa | +0,16 mHa |
| **$R\cdot C$ entero** | **+0,38 mHa** | **+0,030 mHa** | **+0,002 mHa** |

La celda de abajo calcula el **sesgo exacto** de cada combinación a 2,4 Å (probabilidades con matriz densidad, sin shots), con el mismo circuito, $\hat H_5$ y $\boldsymbol\theta^*$ de §7.1, y los shots que necesitaría cada una para $\sigma$ = 0,5 mHa. Sirve para **elegir la técnica antes de gastar shots**.
