**Resultado** (ADAPT con 15 operadores, 194 CZ, 8 M de shots por punto; la celda de arriba lo recalcula):

| | VQE (medido) | FCI (referencia) | diferencia |
|---|---|---|---|
| $R_e$ | 0,8873 ± 0,0018 Å | 0,8886 Å | −1,2 mÅ (−0,7σ) |
| $D_e = E(3{,}0\ \text{Å}) - E(R_e)$ | 313,00 ± 0,66 mHa | 313,21 mHa | −0,21 mHa (−0,3σ) |

- **Zona estirada ($d \ge 1{,}8$ Å):** le ganamos a CCSD en **7 de 7** puntos contando la $\sigma$, y tenemos **precisión química en 7 de 7** (error máximo: 0,56 mHa, a 2,4 Å, contra 9–31 mHa de CCSD).
- **Toda la curva:** precisión química en las 18 distancias, con error máximo de 1,15 mHa a 0,7 Å, dentro de 2σ.
- Las dos diferencias con FCI están dentro de 1σ. El sesgo del ajuste (−0,9 mÅ) y el error de ADAPT sin ruido (≤ 0,38 mHa) quedan por debajo de la $\sigma$ de los shots. Para bajar la $\sigma$ de $D_e$ hacen falta más shots en el mínimo y a 3,0 Å, no un ajuste mejor.
- **Comparado con lo mínimo:** con el ajuste cuadrático de antes, los mismos datos darían $R_e$ con +7 mÅ de error (fila `cuadrático 5`). El cuártico lo baja a ~1 mÅ.
