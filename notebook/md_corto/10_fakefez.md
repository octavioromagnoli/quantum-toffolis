## 9 · Ruido de hardware real: FakeFez

**Pregunta.** ¿La receta (5 qubits, 16 CZ, lectura + ZNE) sobrevive al ruido de un chip real, con errores distintos en cada qubit?

`FakeFez` simula el procesador IBM Fez (156 qubits) con su calibración real. Corrimos 2,4 Å con el compacto de 16 CZ:

- **Ubicación:** el circuito entra en una línea del chip **sin SWAP**. Elegimos, entre 856 líneas, la de menos error esperado (la mitad que la mediana). UCCSD ×2 necesitaría ~1940 CZ.
- **Sin mitigar:** ~+76 mHa, peor que CCSD.
- **Con lectura corregida:** ~+49 mHa (el error de las compuertas).
- **ZNE + lectura:** **compatible con FCI** dentro de 1σ (σ ≈ 1,2 mHa), lejos de CCSD (−30,6 mHa).

Igual que con nuestro modelo, hay que corregir **las dos cosas**: la ZNE sola deja +32 mHa.

**Limitación:** sin ruido de espera ni diafonía. Para afirmar precisión química con confianza harían falta ~6× más shots.
