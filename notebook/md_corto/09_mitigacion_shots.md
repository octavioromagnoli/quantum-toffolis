**Sesgo exacto** (2,4 Å): sin mitigar ~+109 mHa (casi todo lectura) → RO ~+14 → ZNE ~+0,3 → **eco + ZNE + RO ~0,01 mHa**.

**Con shots** (`scripts/mitigation_shots.py`, ~20 min por punto; no se recalcula acá). Error contra FCI en mHa:

| d (Å) | sin mitigar | RO | **eco + ZNE lineal + RO** | CCSD |
|---|---|---|---|---|
| 1,8 | +110,4 ± 0,6 | +14,2 ± 0,6 | **+0,3 ± 0,6** | −9,1 |
| 2,4 | +109,0 ± 0,5 | +14,4 ± 0,5 | **−1,0 ± 0,5** | −30,6 |
| 3,0 | +110,5 ± 0,5 | +14,2 ± 0,5 | **+0,1 ± 0,5** | −23,0 |

- **Con ruido y mitigación, precisión química en la zona estirada** con un circuito de 16 CZ, y le gana a CCSD por 9–31 mHa.
- Las mediciones coinciden con el sesgo exacto predicho ($\chi^2$ = 15,4 con 15 grados de libertad): la σ con ruido también está calibrada.
- **Recomendado: eco + ZNE lineal + RO.** Cuesta ~60 M de shots por punto (4× más que solo RO).

**Lecciones:** la corrección de lectura es imprescindible (la lectura pesa más que las 16 CZ) y barata; la post-selección es gratis pero rinde poco (15–26 %); hay que plegar todo el circuito, no solo las CZ ni solo el ansatz.

**Limitación:** el modelo de ruido es propio; en hardware real la ZNE es aproximada.
