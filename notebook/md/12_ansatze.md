## 11 · Comparación de ansätze: ¿cuál conviene y cuándo?  *(14:27, Gonzalo)*

Para cerrar la parte de circuitos, juntamos en una tabla todos los ansätze que probamos. **Sin ruido, en las 18 distancias** (`scripts/compare_ansatze.py`, varios minutos, por eso no se recalcula acá; los $\boldsymbol\theta$ se optimizan como en la curva y FCI se usa solo para reportar el error; CZ y grupos a 2,4 Å):

| ansatz | qubits | θ | CZ todos con todos | CZ en línea | grupos | error máx. | error máx. ($d \ge 1{,}8$ Å) | $\langle S^2\rangle$ máx. |
|---|---|---|---|---|---|---|---|---|
| UCCSD ×1 | 8 | 26 | 292 | 804 | 37 | 4,46 mHa | 4,46 mHa | 0,85 |
| UCCSD ×2 (lo mínimo) | 8 | 52 | 582 | 1648 | 37 | 0,001 mHa | < 0,001 mHa | ≈ 0 |
| ADAPT 13 | 8 | 13 | 185 | 370 (fSWAP) | 37 | 0,59 mHa | 0,59 mHa | 0,013 |
| **ADAPT 15 (la curva)** | 8 | 15 | 194 | **378 (fSWAP)** | 37 | **0,38 mHa** | 0,38 mHa | ≤ 7·10⁻³ |
| ADAPT 15 reducido | 6 | 15 | 194 | ~510 | 37 | 0,38 mHa | 0,38 mHa | ≤ 7·10⁻³ |
| ADAPT 15 reducido | 5 | 15 | 172 | 427 | 35 | 0,38 mHa | 0,38 mHa | ≤ 7·10⁻³ |
| ADAPT-QEB 15 | 8 | 15 | 159 | 453 | 37 | 0,81 mHa | 0,81 mHa | 0,07–0,14 |
| compacto, $L$ = 5 | 5 | 30 | 10 | 10 | 35 | 56 mHa | 50 mHa | 0,025 |
| compacto, $L$ = 6 | 5 | 35 | 12 | 12 | 35 | 3,8 mHa | 3,8 mHa | 3·10⁻³ |
| **compacto, $L$ = 7** | 5 | 40 | **14** | **14** | 35 | **< 0,001 mHa** | < 0,001 mHa | ≈ 0 |
| compacto, $L$ = 8 (el que usamos con ruido) | 5 | 45 | 16 | 16 | 35 | < 0,001 mHa | < 0,001 mHa | ≈ 0 |

**Con ruido, a 2,4 Å** (`scripts/noise_ansatze.py`; modelo de §7.1 con qubits en línea; 500 k shots sin mitigar y 2 M con ZNE de plegado global, $\lambda$ = 1, 3; CCSD − FCI = −30,6 mHa). El sesgo exacto de §7.1 coincide con estas mediciones dentro de 2σ:

| ansatz | CZ (línea) | compuertas de 1 qubit | error sin ruido | con ruido, sin mitigar | con ZNE |
|---|---|---|---|---|---|
| compacto 5q, $L$ = 8 | 16 | 141 | 0,0 mHa | +12,9 ± 2,2 mHa | **+0,2 ± 2,2 mHa** |
| ADAPT 15 + fSWAP | 313 | 1377 | 0,2 mHa | +178 ± 3 mHa | +70 ± 3 mHa |
| UCCSD ×2 | 1546 | 3801 | 0,0 mHa | +548 ± 2 mHa | — |

**Qué muestran las dos tablas**

- **Sin ruido, el ansatz químico gana en significado y escala:** ADAPT 15 llega a 0,38 mHa con 3,5 veces menos parámetros y 3 veces menos CZ que UCCSD ×2. Conserva $N$ exactamente y $S^2$ casi exactamente, sus parámetros tienen sentido físico (cada uno es una excitación) y se puede justificar sin FCI con la varianza. QEB ahorra CZ con conectividad total, pero en línea queda peor que el ruteo fSWAP y rompe más el espín.
- **Reducir qubits no le sirve a ADAPT:** la energía es idéntica en 8, 6 y 5 qubits (la reducción está bien hecha), pero las CZ casi no bajan, y en línea quedan peor que con el ruteo fSWAP de 8 qubits. La reducción mezcla los modos y las excitaciones dejan de ser locales.
- **La reducción sí rinde con el compacto:** en 5 qubits, un ansatz en ladrillo de $R_Y$ + CZ entre vecinos (sin SWAP) llega a FCI con **14 CZ** ($L$ = 7). Con 12 ya falla (3,8 mHa). Es *hardware-efficient*, no químico: no conserva $N$ ni $S^2$ por construcción, solo por la penalización.
- **Con ruido, la cantidad de compuertas decide todo:** los tres ansätze son igual de buenos sin ruido, pero con 0,1 % por CZ solo el de 16 CZ le gana a CCSD. La ZNE no alcanza a salvar 313 CZ.

**Cuál usar.** Sin ruido (la curva central), **ADAPT 15**: químico, pocos parámetros, justificable sin FCI y con un método que escala a moléculas más grandes. Con ruido (los extras), **el compacto de 5 qubits**.
