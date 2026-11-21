**Qué dio.** *Sorted insertion* necesita **1,6–1,9 veces menos shots** que el coloreo para la misma $\sigma$, sin cambiar el circuito. Con los mismos 8 M de shots, la $\sigma$ predicha baja de ~0,70 a ~0,52 mHa en la zona estirada, que es lo que después se mide en la curva (§5). GC necesita otras 2–3 veces menos shots, pero agrega de 4 a 23 CZ por grupo (de 4 a 77 con qubits en línea).

**Decisión.** La curva se mide con **QWC + *sorted insertion* + dos etapas**, con 8 M de shots por punto ($\sigma$ = 0,50–0,60 mHa). Con ruido, cada CZ extra cuesta precisión, así que QWC además permite post-seleccionar por $N_\alpha$ y $N_\beta$ en los grupos que miden solo $Z$ (§8). GC queda como opción para una simulación sin ruido: 4 M de shots darían la misma $\sigma$.

**Por qué funcionó.** $K = (\sum_g \sqrt{V_g})^2$ crece con la cantidad de grupos y con la varianza de cada uno. *Sorted insertion* arma menos grupos y pone juntos los términos de coeficiente grande. Como dentro de un grupo se miden en el mismo shot, sus covarianzas entran en $V_g$ y en parte se compensan. Con 20 % menos grupos, $K$ baja casi a la mitad.

### 4.3 ¿La barra de error es honesta? Calibración de σ

Repetimos la medición del **mismo circuito** con semillas independientes y miramos

$$z = \frac{E - E_{\rm exacta}}{\sigma}, \qquad E_{\rm exacta} = \langle\psi\rvert\hat H\lvert\psi\rangle \text{ del circuito medido}.$$

Si $\sigma$ está bien estimada, $z \sim \mathcal N(0,1)$. Hicimos 500 repeticiones por configuración (300 para $\sum\lvert c\rvert$), a 2,0 y 3,0 Å: **17,8 mil millones de shots** en total (`scripts/calib_sigma.py`, ~35 min en 12 núcleos, por eso no se recalcula acá). También calculamos $\sigma$ por **bootstrap**: re-muestrear los shots de cada grupo y recalcular $E$, sin suponer normalidad.

| configuración | d (Å) | σ media (mHa) | σ predicha $\sqrt{K/N}$ | media de z | std de z [IC 95 %] | σ_boot / σ | \|z\| > 2 |
|---|---|---|---|---|---|---|---|
| QWC, 8 M, dos etapas | 2,0 | 0,519 | 0,519 | −0,02 ± 0,05 | 1,01 [0,95; 1,07] | 0,997 | 5,0 % |
| QWC, 8 M, dos etapas | 3,0 | 0,500 | 0,500 | +0,04 ± 0,05 | 1,00 [0,94; 1,06] | 0,998 | 4,4 % |
| QWC, 1 M, dos etapas | 2,0 | 1,467 | 1,467 | +0,04 ± 0,05 | 1,03 [0,97; 1,09] | 0,999 | 4,8 % |
| QWC, 1 M, dos etapas | 3,0 | 1,413 | 1,413 | −0,02 ± 0,05 | 1,04 [0,98; 1,11] | 0,999 | 6,2 % |
| QWC, 8 M, Σ\|c\| | 2,0 | 0,679 | 0,679 | −0,16 ± 0,06 | 0,98 [0,91; 1,06] | 0,997 | 5,3 % |
| QWC, 8 M, Σ\|c\| | 3,0 | 0,681 | 0,681 | −0,02 ± 0,06 | 1,03 [0,95; 1,12] | 1,001 | 5,0 % |
| GC, 4 M, dos etapas | 2,0 | 0,507 | 0,507 | +0,01 ± 0,04 | 0,96 [0,91; 1,03] | 0,996 | 3,4 % |
| GC, 4 M, dos etapas | 3,0 | 0,424 | 0,424 | +0,01 ± 0,05 | 1,04 [0,98; 1,11] | 0,994 | 5,8 % |

- **$\sigma$ está calibrada.** La dispersión de $z$ es compatible con 1 en las 8 configuraciones. Con las 3600 repeticiones juntas: media de $z$ = −0,006 ± 0,017, dispersión = 1,011 [0,988; 1,035], $\lvert z\rvert > 2$ en el 5,0 % de los casos (lo normal es 4,6 %), 11 casos con $\lvert z\rvert > 3$ (se esperan ~10) y la prueba KS contra $\mathcal N(0,1)$ da $p$ = 0,78. No hay correlación entre repeticiones seguidas ($\lvert r\rvert$ < 0,12).
- **$\sigma = \sqrt{K/N}$ predice la $\sigma$ medida hasta el tercer dígito**, y de 1 M a 8 M de shots la $\sigma$ baja en $\sqrt 8$ = 2,83 (medido: 2,83). Por eso la usamos para planificar.
- **El 1,16 de antes era ruido estadístico** (hipótesis b confirmada): con 60 repeticiones, aunque $\sigma$ sea exacta, la dispersión de $z$ cae en [0,82; 1,18] el 95 % de las veces.
- **El bootstrap coincide** con la $\sigma$ analítica dentro de 0,7 %.
- La media de $z$ = −0,16 ± 0,06 ($\sum\lvert c\rvert$ a 2,0 Å) es un desvío de 2,7σ entre unas 18 pruebas, algo que por azar pasa ~10 % de las veces. El mismo circuito medido en dos etapas da −0,02 ± 0,05.

### 4.4 Semillas: el error que sí encontramos

Al principio cada grupo usaba semillas consecutivas de Aer, y vimos **errores correlacionados entre grupos**, con desvíos de 4σ: la $\sigma$ quedaba subestimada. Ahora cada grupo y cada etapa tienen su propia semilla, derivada con `numpy.random.SeedSequence`. Más tarde vimos ~3σ seguidos en 0,7–0,9 Å con 200 k shots: las 18 distancias usaban **la misma** semilla, y el ruido de shots quedaba correlacionado entre puntos vecinos. Eso sesga el ajuste de $R_e$ y achica su $\sigma$ de Monte Carlo, que supone errores independientes. Ahora cada distancia y cada diagnóstico tienen su propia semilla.

### 4.5 Qué no funcionó

- **Neyman con las varianzas del estado HF**: 11 veces peor que el piloto. HF es un estado de la base computacional, así que da varianza cero a todos los grupos que miden solo $Z$ y les asigna muy pocos shots.
- **Reusar shots entre grupos** (cada base QWC también estima los Pauli compatibles de otros grupos): sobre el coloreo ganaba 1,4–2,3×, pero sobre *sorted insertion* casi no suma, y en la zona estirada empeora un poco porque el promedio ignora las covarianzas. No lo implementamos.
- **Calibrar con 60 repeticiones**: no alcanza para distinguir 1,16 de 1.

**Limitación.** La calibración es **sin ruido**: verifica la estadística de los shots, no el sesgo del ruido de compuertas, que se trata en §7–§8.
