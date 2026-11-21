### 13.1 La línea de tiempo

| hora | qué hicimos | quién | qué cambió |
|---|---|---|---|
| 11:09 | **Lo mínimo:** UCCSD ×2, 8 M de shots, ajuste cuadrático | Joaquín | Cumple la consigna (7/7 contra CCSD), pero con 1648 CZ en línea |
| 11:42 | Ajuste elegido por su sesgo contra el mínimo FCI verdadero; diagnósticos medidos con shots | Joaquín | $R_e$: sesgo de +7,8 a −0,9 mÅ |
| 11:46 | **ADAPT-VQE** con lista común (13 operadores) | Gonzalo | 582 → 185 CZ |
| 11:49 | ***Sorted insertion*** + estudio analítico de shots | Santino | 1,6–1,9× menos shots para la misma σ |
| 12:02 | Semilla propia por punto; medición por conmutación general | Joaquín, Santino | Ruido de shots independiente entre distancias |
| 12:07 | **Tapering a 5 qubits**, ansatz compacto de 16 CZ, modelo de ruido; otra geometría | Octavio | Primer circuito que sobrevive al ruido |
| 12:20 | **Ruteo fSWAP** para qubits en línea | Gonzalo | 515 → 370 CZ en línea |
| 12:31 | **FakeFez** con lectura + ZNE | Octavio | Compatible con FCI con ruido de hardware real |
| 12:40 | **Calibración de σ** (3600 repeticiones, 17,8 G shots) | Santino | La barra de error es honesta (std de z = 1,01) |
| 13:03 | Mitigación: post-selección, lectura combinable, **plegado de R·C** | Santino | Sesgo de la ZNE de 0,16 a 0,002 mHa |
| 13:04 | ADAPT de **15 operadores** (pares de espín) y doble arranque | Gonzalo | Error de 0,59 a 0,38 mHa; $S^2$ a la mitad |
| 13:23 | Reducción a **6 qubits** (lo que pide el extra) | Octavio | Tabla 8 → 6 → 5 |
| 13:41 | ADAPT 15 integrado como ansatz de la curva | Joaquín | 194 CZ (378 en línea) en las 18 distancias |
| 13:44 | **Mitigación con shots** | Santino | 6/6 con precisión química con ruido |
| 14:27 | Comparación de ansätze sin y con ruido | Gonzalo | El criterio "cuál usar y cuándo" |
| 14:37 | **Escalado** H₄ → H₁₄ → FeMoco | Santino | Hasta dónde llega el método |

### 13.2 Las decisiones y por qué mejoraron el resultado

| decisión | por qué funcionó | alternativa descartada |
|---|---|---|
| **ADAPT-VQE** con una lista común de 15 operadores | En H₄ la correlación la llevan unas pocas dobles $\alpha\beta$, y ADAPT las encuentra por gradiente. El mismo circuito en toda la curva, con 3 veces menos CZ que UCCSD ×2 | UCCSD ×2: exacto, pero con 582 CZ (1648 en línea) y 52 parámetros |
| Cortar ADAPT con la **varianza** y completar los **pares de espín** | La varianza no necesita FCI y cae un orden de magnitud justo cuando el error lo hace; los pares de espín arreglan $S^2$ y bajan el error a 0,38 mHa | Cortar por gradiente < 10⁻⁴: se frena en mesetas |
| **Ruteo fSWAP** propio | Un fSWAP mueve modos fermiónicos con su signo, en 2 CZ, y no hay que deshacerlo | SWAP del transpiler: ~2,8× más CZ en línea |
| **QWC con *sorted insertion*** + Neyman en dos etapas | Menos grupos y covarianzas que se compensan: 1,6–1,9× menos shots, sin CZ extra | Coloreo QWC (46 grupos); GC: menos shots pero agrega CZ |
| σ **desde las muestras**, semilla propia por grupo y por punto, calibrada con 3600 repeticiones | Barras honestas (std de z = 1,01) y errores independientes entre puntos | EstimatorV2 con `precision` fija; semillas consecutivas (desvíos de 4σ) |
| $\boldsymbol\theta^*$ **sin ruido** con gradiente analítico y doble arranque | Lo permite la consigna; segundos por punto y sin mínimos malos | SPSA con shots: lento y ruidoso; solo *warm start*: mínimos malos |
| **Cuártico en 7 puntos** elegido por su sesgo contra el mínimo FCI exacto | Mide el error del método de ajuste sin mirar el VQE | Cuadrático en 5 puntos: +7,8 mÅ de sesgo; comparar contra FCI ajustado igual |
| Confiabilidad **medida**: eco de Loschmidt, $S^2$, profundidad, límite físico | Cuenta eventos raros (resuelve $10^{-6}$ Ha²) y detecta el estado de espín equivocado | Varianza por grupos sola: $\sigma \approx 2\cdot10^{-3}$ Ha², no separa ansätze |
| Con ruido: **tapering a 5 qubits + compacto de 16 CZ** | Las simetrías son información conocida; sacarlas deja un espacio de 32 estados que un circuito nativo de la línea recorre sin SWAP | ADAPT o UCCSD con ruido: 181 y 548 mHa de error |
| Mitigación: **lectura + eco + ZNE con plegado de R·C** | Cada técnica ataca un sesgo distinto; plegar también la rotación de medición escala todo el ruido | ZNE sola (no toca la lectura); plegar solo las CZ (4 mHa de sesgo) |

### 13.3 Hasta dónde llega el método (limitaciones)

- **La curva central es sin ruido de compuertas.** La $\sigma$ es solo estadística. Con 0,1 % de error por CZ, el ansatz de la curva (194 CZ, 378 en línea) tendría un sesgo de ~180 mHa. Con ruido, el resultado que vale es el del compacto de 5 qubits con mitigación, que **no escala** a moléculas grandes.
- **Los parámetros se buscan sin ruido**, como permite la consigna. La curva mide qué tan expresivo es el ansatz más el ruido de los shots, no la dificultad de optimizar con mediciones ruidosas.
- **Shots:** 8 M por punto (144 M en la curva) para $\sigma \approx$ 0,5 mHa; la mitigación multiplica eso por ~4.
- **El eco de Loschmidt** es un diagnóstico de simulador (más de 4000 CZ). En hardware, lo equivalente es la varianza por grupos con más shots.
- **FakeFez** no modela los tiempos muertos ni la diafonía, y su $\sigma \approx$ 1,2 mHa no alcanza para afirmar precisión química.
- **Escalado (§12):** los conteos crecen como potencias de $N$, pero los shots y la probabilidad de un circuito sin errores vuelven inviable un VQE así para FeMoco. Las estimaciones serias para FeMoco usan estimación de fase con corrección de errores.

### 13.4 Lo que no hicimos y haríamos después

- **VQE contra SQD** (extra de la consigna): no llegamos. Con 8 qubits, SQD es trivial (hasta muestras al azar dan la respuesta), así que el experimento interesante es el de §12: un tope al subespacio en cadenas más largas.
- **Un pool adaptado al espín** (pares de dobles $\alpha\beta$ con un mismo parámetro) para que $S^2$ sea exactamente 0, y operadores CEO de 9 CNOT para bajar más las CZ.
- **GC en 5 qubits** para el compacto con ruido: menos shots a cambio de unas pocas CZ.
- **Hardware real**: la receta de §9 (línea elegida por calibración, lectura + ZNE) está lista para un backend de IBM con `QiskitRuntimeService`.
