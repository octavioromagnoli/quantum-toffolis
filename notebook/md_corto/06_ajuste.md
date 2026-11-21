### 5.1 $R_e$ y $D_e$, y el problema 3: el ajuste tenía sesgo

$R_e$ sale de ajustar un polinomio a las energías medidas cerca del mínimo; $D_e = E(3{,}0\ \text{Å}) - E(R_e)$. Las incertidumbres salen por Monte Carlo (2000 réplicas re-muestreando cada energía con su σ).

**El problema.** Antes usábamos una parábola en 5 puntos y la comparábamos con FCI ajustado igual, así que el error del ajuste se cancelaba y no se veía. La curva es asimétrica, y una parábola corre el mínimo.

**Cómo elegimos el ajuste sin mirar el VQE:** aplicamos cada modelo a FCI en la misma grilla y lo comparamos con el mínimo FCI verdadero. La diferencia es el sesgo del ajuste (tabla de abajo):

- **parábola en 5 puntos** (la de antes): **+7,8 mÅ** de sesgo, unas 3σ;
- **cuártico en 7 puntos** (0,7–1,1 Å): −0,9 mÅ y −16 µHa, el único con sesgo chico en $R_e$ y en la energía. **Es el que usamos.**

El cuártico también recupera el mínimo de una curva de Morse conocida (`tests/test_analysis.py`), así que la elección no depende de FCI.
