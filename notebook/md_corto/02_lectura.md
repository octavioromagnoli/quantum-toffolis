### 1.4 Cómo leer el hamiltoniano

Los 185 términos se agrupan en tres tipos:

| términos | cuántos | qué significan |
|---|---|---|
| $I$, $Z_p$, $Z_pZ_q$ | 37 | energía de cada configuración (diagonal) |
| 2 X/Y | 56 | **un** electrón salta de orbital |
| 4 X/Y | 92 | **dos** electrones saltan a la vez |

**Al estirar**, los términos diagonales se achican (los $Z_p$ pasan de 1,64 a 0,31 Ha) pero los saltos dobles no (2,13 → 2,22 Ha). Muchas configuraciones quedan casi con la misma energía y bien conectadas: el estado exacto se vuelve una mezcla. Eso es la **correlación estática**, lo que rompe a CCSD.

Todos los términos conservan $N_\alpha$ y $N_\beta$; lo aprovechamos en §7 para sacar qubits.

### 1.5 Validación

En las 18 distancias (primera línea de la salida de abajo):
- $\langle \mathrm{HF}\rvert\hat H\lvert \mathrm{HF}\rangle = E_{\rm HF}$ de PySCF, a menos de $10^{-14}$ Ha.
- La energía más baja de $\hat H$ con 2 electrones $\alpha$ y 2 $\beta$ coincide con $E_{\rm FCI}$, y es un singulete.

Así, cualquier diferencia entre VQE y FCI viene del ansatz o de las mediciones, no del hamiltoniano.
