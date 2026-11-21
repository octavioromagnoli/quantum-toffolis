### 5.1 $R_e$ y $D_e$, y el problema 3: el ajuste tenía sesgo
El mínimo sale de ajustar un polinomio a las energías medidas cerca de 0,9 Å, pesando cada punto con $1/\sigma_i$:

$$E(d) \approx \sum_{j=0}^{n} a_j\, d^{\,j}, \qquad R_e:\ E'(R_e) = 0,\ E''(R_e) > 0, \qquad D_e = E(3{,}0\ \text{Å}) - E(R_e), \qquad \sigma_{D_e} = \sqrt{\sigma_{3,0}^2 + \sigma_{E(R_e)}^2}.$$

La incertidumbre sale por **Monte Carlo**: re-muestreamos cada energía con su $\sigma$ ($E_i \to E_i + \mathcal N(0, \sigma_i^2)$), volvemos a ajustar (incluida la elección de la ventana de puntos) y tomamos la dispersión de $R_e$ y de $E(R_e)$ en 2000 réplicas.

**Observación.** En lo mínimo usábamos un cuadrático en 5 puntos y comparábamos con FCI *ajustado con el mismo modelo*. Así el error del ajuste se cancela en la comparación y no se ve.

**Hipótesis.** La curva es asimétrica (sube más rápido hacia la izquierda), y una parábola corre el mínimo hacia afuera. Para medir ese sesgo sin mirar el VQE, aplicamos cada modelo a las energías **FCI de la misma grilla** y lo comparamos con el **mínimo FCI verdadero**, que sacamos minimizando PySCF en $d$ ($R_e$ = 0,8886 Å). Esa diferencia es el sesgo del *método de ajuste*. La celda de abajo prueba 7 modelos (columnas `sesgo R_e` y `sesgo E_min`):

- el cuadrático en 5 puntos (el de lo mínimo) tiene **+7,8 mÅ** de sesgo en $R_e$, unas 3σ: ese era el problema;
- el cúbico en 7 puntos casi no tiene sesgo en $R_e$, pero sí en la energía del mínimo (−575 µHa), que entra en $D_e$;
- el spline no tiene sesgo, pero pasa por cada punto con su ruido y queda 3 veces más ruidoso;
- **el cuártico en 7 puntos (0,7–1,1 Å)** es el único con sesgo chico en $R_e$ **y** en $E_{\min}$ a la vez (−0,9 mÅ y −16 µHa), con una $\sigma$ de las más chicas. **Lo elegimos.**

La elección tampoco depende de FCI: el cuártico en 7 puntos recupera el mínimo de una curva de Morse conocida en esta misma grilla, con error < 2 mÅ en $R_e$ y < 0,1 mHa en $E_{\min}$ (`tests/test_analysis.py`). La comparación con FCI solo elige un método numérico; no corrige ningún resultado.
