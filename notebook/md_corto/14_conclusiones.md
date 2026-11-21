### 13.1 Qué hicimos y qué cambió

| qué hicimos | qué cambió |
|---|---|
| **Lo mínimo:** UCCSD ×2 | cumple la consigna, pero con 1648 CZ en línea |
| **ADAPT-VQE** (15 operadores) + **ruteo fSWAP** | 582 → 194 CZ; 1648 → 378 en línea; error 0,38 mHa |
| ***Sorted insertion*** + Neyman | 1,6–1,9× menos shots para la misma σ |
| **Calibración de σ** (3600 repeticiones) | la barra de error es honesta |
| Ajuste elegido por su sesgo | sesgo en $R_e$: +7,8 → −0,9 mÅ |
| Confiabilidad **medida** con shots | eco, simetrías y profundidad, sin FCI |
| **5 qubits + 16 CZ** + mitigación | precisión química **con ruido** |
| FakeFez, geometrías, ansätze, FeMoco | la receta funciona y sabemos hasta dónde llega |

### 13.2 Decisiones clave

| decisión | por qué | descartado |
|---|---|---|
| ADAPT 15, una lista para toda la curva | pocas dobles llevan la correlación; 3× menos CZ | UCCSD ×2 (582 CZ) |
| ruteo fSWAP propio | mueve modos con su signo, en 2 CZ | SWAP del transpiler (~2,8× más CZ) |
| QWC + *sorted insertion* | menos shots sin CZ extra | coloreo; GC (agrega CZ) |
| σ desde las muestras, semillas independientes | barras honestas | σ fija; semillas correlacionadas |
| $\boldsymbol\theta^*$ sin ruido, dos arranques | rápido y sin mínimos malos | optimizar con shots (lento) |
| cuártico en 7 puntos | sesgo mínimo, elegido sin mirar el VQE | parábola (+7,8 mÅ) |
| con ruido: 5 qubits + 16 CZ + lectura + eco + ZNE | cada técnica ataca un error distinto | ADAPT/UCCSD con ruido (181 y 548 mHa) |

### 13.3 Limitaciones

- La curva central es **sin ruido de compuertas**: con ruido, el ansatz de la curva tendría ~180 mHa de sesgo. Con ruido vale el compacto de 5 qubits, que no escala.
- Los parámetros se buscan sin ruido, como permite la consigna.
- 8 M de shots por punto; la mitigación multiplica eso por ~4.
- El eco de Loschmidt es un diagnóstico de simulador (más de 4000 CZ).
- FakeFez no modela tiempos muertos ni diafonía.

### 13.4 Qué haríamos después

- **VQE contra SQD** (extra de la consigna) en cadenas más largas: en 8 qubits SQD es trivial.
- Un pool adaptado al espín, para que $S^2$ sea exactamente 0, y operadores más baratos en CZ.
- Correr la receta de §9 en hardware real de IBM.
