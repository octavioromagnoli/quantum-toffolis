## 11 · Comparación de ansätze: ¿cuál conviene y cuándo?

Todos los ansätze que probamos (`scripts/compare_ansatze.py` y `scripts/noise_ansatze.py`; no se recalculan acá).

**Sin ruido, 18 distancias** (CZ a 2,4 Å):

| ansatz | qubits | CZ | CZ en línea | error máx. | $\langle S^2\rangle$ máx. |
|---|---|---|---|---|---|
| UCCSD ×1 | 8 | 292 | 804 | 4,46 mHa | 0,85 |
| UCCSD ×2 (lo mínimo) | 8 | 582 | 1648 | < 0,001 mHa | ≈ 0 |
| **ADAPT 15 (la curva)** | 8 | 194 | **378** | **0,38 mHa** | ≤ 0,007 |
| ADAPT 15 reducido | 5 | 172 | 427 | 0,38 mHa | ≤ 0,007 |
| ADAPT-QEB 15 | 8 | 159 | 453 | 0,81 mHa | 0,14 |
| **compacto, 7 capas** | 5 | **14** | **14** | < 0,001 mHa | ≈ 0 |

**Con ruido, a 2,4 Å** (CCSD: −30,6 mHa):

| ansatz | CZ en línea | sin mitigar | con ZNE |
|---|---|---|---|
| compacto 5 qubits | 16 | +12,9 ± 2,2 mHa | **+0,2 ± 2,2 mHa** |
| ADAPT 15 + fSWAP | 313 | +178 ± 3 mHa | +70 ± 3 mHa |
| UCCSD ×2 | 1546 | +548 ± 2 mHa | — |

**Cuál usar:**
- **Sin ruido → ADAPT 15:** químico (cada parámetro es una excitación), conserva las simetrías, se justifica sin FCI y escala a moléculas más grandes.
- **Con ruido → el compacto de 5 qubits:** la cantidad de CZ decide todo. Pero es una solución para H₄: no conserva las simetrías por sí solo y no escala.
- Reducir qubits **no** le sirve a ADAPT (las CZ casi no bajan); sí al compacto.
