## 4 · Resultados: curva, R_e y D_e  *(Joaquín)*

Para cada distancia: hamiltoniano → ansatz → θ* (sin ruido, con warm start desde el punto anterior) → **energía medida** con shots. La tabla tiene, por punto, la energía, su σ, los shots y las CZ: todos con todos (`cz_all`), en línea con el transpiler de Qiskit (`cz_line`) y en línea con nuestro ruteo fSWAP (`cz_line_fswap`, sección 2). El circuito ruteado prepara el mismo estado con los modos reordenados; medimos el circuito sin rutear con Ĥ original porque Ĥ reordenado necesita ~1,9× más grupos QWC.

**Cómo sacamos R_e y D_e.** El mínimo sale de ajustar las energías medidas cerca de 0,9 Å. La σ se obtiene por Monte Carlo: re-muestreamos cada energía con su σ y re-ajustamos, incluida la elección de los puntos.

**Cómo elegimos el ajuste.** Probamos 7 modelos: cuadrático en 3 y 5 puntos, cúbico en 6 y 7, cuártico en 7 y 8, y spline cúbico entre 0,6 y 1,6 Å. Para elegir **sin mirar el VQE**, aplicamos cada modelo a las energías FCI de la *misma grilla* y lo comparamos con el mínimo FCI verdadero, que sacamos minimizando PySCF en d (R_e = 0,8886 Å). Esa diferencia es el sesgo del *método de ajuste*:

| modelo | sesgo en R_e | sesgo en E_min | σ(R_e) con nuestras barras |
|---|---|---|---|
| cuadrático, 5 puntos (el que usábamos) | **+7,8 mÅ** | +126 µHa | 2,2 mÅ |
| cúbico, 7 puntos | +0,4 mÅ | **−575 µHa** | 2,5 mÅ |
| **cuártico, 7 puntos (0,7–1,1 Å)** | **−0,9 mÅ** | **−16 µHa** | 2,3 mÅ |
| spline 0,6–1,6 Å | −0,5 mÅ | −5 µHa | 7,6 mÅ |

El cuadrático tiene un sesgo de 3σ: la curva es asimétrica (sube más rápido hacia la izquierda) y una parábola corre el mínimo hacia afuera. El spline no tiene sesgo, pero pasa por cada punto con su ruido y queda 3 veces más ruidoso. La elección tampoco depende de FCI: el cuártico en 7 puntos recupera el mínimo de una curva de Morse conocida en esta misma grilla, con error < 2 mÅ en R_e y < 0,1 mHa en E_min (`tests/test_analysis.py`). La comparación con FCI solo elige un método numérico; no corrige ningún resultado.

**Elegimos el cuártico en 7 puntos**: es el único con sesgo chico en R_e *y* en E_min a la vez, y su σ es de las más chicas. Las diferencias entre modelos de bajo sesgo (≈ 2 mÅ en R_e, ≲ 0,1 mHa en D_e) quedan dentro de la σ estadística.

**Resultado** (corrida completa, ADAPT 15 operadores, 194 CZ, 8M shots por punto):

| | VQE (medido) | FCI | diferencia |
|---|---|---|---|
| R_e | 0,8873 ± 0,0018 Å | 0,8886 Å | −1,2 mÅ (−0,7σ) |
| D_e = E(3,0 Å) − E(R_e) | 313,00 ± 0,66 mHa | 313,21 mHa | −0,21 mHa (−0,3σ) |

Las dos diferencias están dentro de 1σ. El sesgo del ajuste (−0,9 mÅ) y el error de ADAPT sin ruido (≤ 0,38 mHa) quedan debajo de la σ de shots.

En la zona estirada (d ≥ 1,8 Å) le ganamos a CCSD en 7/7 puntos contando la σ, y tenemos precisión química en 7/7 (error máximo: 0,56 mHa, a 2,4 Å). En toda la curva el error máximo es 1,15 mHa (0,7 Å). La σ por punto es 0,50–0,60 mHa. La σ de D_e está dominada por las dos energías que entran en la resta. Para bajarla hacen falta más shots en el mínimo y a 3,0 Å, no un ajuste mejor.

> Las celdas de abajo recalculan todo esto en cada ejecución. Si cambian el ansatz o los shots, los números de esta sección se actualizan ahí.
