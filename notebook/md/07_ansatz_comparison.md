## Comparación de ansätze  *(Gonzalo)*

**Sin ruido, en las 18 distancias** (`scripts/compare_ansatze.py`). Los θ se optimizan como en la curva del notebook principal. FCI se usa solo para reportar el error. S² y las CZ en línea varían un poco entre corridas (el optimizador puede llegar a θ distintos con la misma energía, y el transpiler usa semillas). CZ y grupos a 2,4 Å (el circuito es el mismo en toda la curva).

| ansatz | qubits | θ | CZ todos con todos | CZ en línea | grupos | error máx. | error máx. (d ≥ 1,8 Å) | S² máx. |
|---|---|---|---|---|---|---|---|---|
| UCCSD ×1 | 8 | 26 | 292 | 804 | 37 | 4,46 mHa | 4,46 mHa | 0,85 |
| UCCSD ×2 | 8 | 52 | 582 | 1648 | 37 | 0,001 mHa | < 0,001 mHa | ≈ 0 |
| ADAPT 13 | 8 | 13 | 185 | 370 (fSWAP) | 37 | 0,59 mHa | 0,59 mHa | 0,013 |
| **ADAPT 15 (curva central)** | 8 | 15 | 194 | **378 (fSWAP)** | 37 | **0,38 mHa** | 0,38 mHa | ≤ 7·10⁻³ |
| ADAPT 15 reducido | 6 | 15 | 194 | ~510 | 37 | 0,38 mHa | 0,38 mHa | ≤ 7·10⁻³ |
| ADAPT 15 reducido | 5 | 15 | 172 | 427 | 35 | 0,38 mHa | 0,38 mHa | ≤ 7·10⁻³ |
| ADAPT-QEB 15 | 8 | 15 | 159 | 453 | 37 | 0,81 mHa | 0,81 mHa | 0,07-0,14 |
| compacto (Octavio), L = 5 | 5 | 30 | 10 | 10 | 35 | 56 mHa | 50 mHa | 0,025 |
| compacto, L = 6 | 5 | 35 | 12 | 12 | 35 | 3,8 mHa | 3,8 mHa | 3·10⁻³ |
| **compacto, L = 7** | 5 | 40 | **14** | **14** | 35 | **< 0,001 mHa** | < 0,001 mHa | ≈ 0 |
| compacto, L = 8 (default) | 5 | 45 | 16 | 16 | 35 | < 0,001 mHa | < 0,001 mHa | ≈ 0 |

Qué muestra la tabla:

- **Ansatz químico:** ADAPT 15 llega a 0,38 mHa con 3,5 veces menos parámetros y 3 veces menos CZ que UCCSD ×2. QEB ahorra CZ con conectividad total, pero en línea queda peor que el ruteo fSWAP y rompe más el espín.
- **Reducir qubits no le sirve a ADAPT:** la energía es idéntica en 8, 6 y 5 qubits (la reducción está bien hecha), pero las CZ casi no bajan. En línea quedan peor que el ruteo fSWAP de 8 qubits, porque la reducción mezcla los modos y las excitaciones dejan de ser locales.
- **La reducción rinde con el compacto:** en 5 qubits, una brickwork de RY + CZ entre vecinos (sin SWAP) llega a FCI con **14 CZ** (L = 7). Con 12 CZ ya falla (3,8 mHa). No conserva N ni S² por construcción (solo por la penalización de S² en la optimización): es "hardware-efficient", no químico.

**Con ruido, a 2,4 Å** (`scripts/noise_ansatze.py`; modelo de `noise.py`, qubits en línea; 500 k shots sin mitigar, 2 M con ZNE de plegado global con λ = 1, 3). CCSD − FCI = −30,6 mHa.

| ansatz | CZ (línea) | compuertas de 1 qubit | error sin ruido | con ruido, sin mitigar | con ZNE |
|---|---|---|---|---|---|
| compacto 5q L = 8 | 16 | 141 | 0,0 mHa | +12,9 ± 2,2 mHa | **+0,2 ± 2,2 mHa** |
| ADAPT 15 + fSWAP | 313 | 1377 | 0,2 mHa | +178 ± 3 mHa | +70 ± 3 mHa |
| UCCSD ×2 | 1546 | 3801 | 0,0 mHa | +548 ± 2 mHa | — |

**Con ruido, la cantidad de compuertas decide todo.** Los tres ansätze son igual de buenos sin ruido, pero con 0,1 % por CZ solo el de 16 CZ le gana a CCSD; ADAPT y UCCSD ×2 quedan a cientos de mHa. La ZNE no alcanza a salvar 313 CZ.

Con muchos shots, el compacto con mitigación llega a **precisión química con ruido, contando la σ** (sección de mitigación de Santino, `scripts/mitigation_shots.py`: ZNE lineal + post-selección por simetría + corrección de lectura, 54 M de shots, con 2 % de error de lectura): −0,48 ± 0,54 mHa a 1,8 Å, −0,09 ± 0,50 a 2,4 Å y +0,31 ± 0,48 a 3,0 Å.

**Cuál usar:** sin ruido (la curva central), ADAPT 15. Es químico (conserva N exacto y casi S²), tiene pocos parámetros con significado físico, y se puede justificar sin FCI con la varianza. Con ruido (extras), el compacto de 5 qubits.

