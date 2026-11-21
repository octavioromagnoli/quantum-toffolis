**Qué dio.** *Sorted insertion* necesita **1,6–1,9× menos shots** que el coloreo para la misma σ: con 8 M de shots, σ baja de ~0,70 a ~0,52 mHa. GC ahorra otras 2–3×, pero agrega CZ.

**Decisión:** QWC + *sorted insertion* + dos etapas, 8 M de shots por punto (σ = 0,50–0,60 mHa). Sin CZ extra, y permite post-seleccionar por simetría con ruido (§8).

### 4.3 ¿La barra de error es honesta?

Medimos el **mismo circuito** 3600 veces con semillas independientes (17,8 mil millones de shots, `scripts/calib_sigma.py`; no se recalcula acá). Si σ está bien, $z = (E - E_{\rm exacta})/\sigma$ tiene que tener dispersión 1.

| configuración | σ (mHa), 2,0 / 3,0 Å | dispersión de $z$, 2,0 / 3,0 Å |
|---|---|---|
| QWC, 8 M, dos etapas (la curva) | 0,519 / 0,500 | 1,01 / 1,00 |
| QWC, 1 M, dos etapas | 1,467 / 1,413 | 1,03 / 1,04 |
| GC, 4 M, dos etapas | 0,507 / 0,424 | 0,96 / 1,04 |

- **σ está calibrada:** con todo junto, dispersión de $z$ = 1,011 [0,988; 1,035] y $\lvert z\rvert > 2$ en el 5,0 % de los casos.
- **$\sqrt{K/N}$ predice la σ medida** hasta el tercer dígito.
- **El 1,16 era ruido estadístico:** con 60 repeticiones, cae en [0,82; 1,18] el 95 % de las veces.
- **El bootstrap coincide** con la σ analítica (diferencia < 0,7 %).

### 4.4 Lo que sí estaba mal: las semillas

- Semillas consecutivas entre grupos daban errores correlacionados (desvíos de 4σ).
- La misma semilla en las 18 distancias correlacionaba puntos vecinos y sesgaba el ajuste de $R_e$.

Ahora cada grupo, etapa y distancia tiene su propia semilla independiente.

**Qué no funcionó:** repartir los shots con las varianzas de Hartree-Fock (11× peor) y reusar shots entre grupos (casi no suma con *sorted insertion*).
