### Extra: ruido de hardware real con FakeFez  *(Octavio)*

`FakeFez` es una copia simulada del procesador IBM Fez (Heron r2, 156 qubits) con su calibración real: el error de cada compuerta, T1/T2 y el error de lectura de cada qubit. Corrimos un punto estirado (2,4 Å) con el ansatz compacto de 5 qubits y 16 CZ (ver la reducción de qubits en la sección 1).

**Cómo lo ubicamos en el chip**

- El ansatz compacto solo usa CZ entre vecinos de una línea, así que entra en el chip sin SWAP. Probamos las 856 líneas de 5 qubits de Fez y elegimos la de menos errores esperados (la suma del error de cada compuerta más el error de lectura, con la calibración del backend): 0,087 contra 0,172 de la mediana, o sea la mitad.
- UCCSD ×2 ubicado en Fez necesita 1939 CZ (por los SWAP). Con ~1 % de error por CZ, ese circuito no sirve.
- Simulamos con el modelo de ruido de Fez restringido a esos 5 qubits, con matriz densidad. Da exactamente el mismo valor esperado que el chip completo de 156 qubits (lo verifica `tests/test_hardware.py`) y es unas 20 veces más rápido.

**Mitigación** (energías medidas con shots, tabla de abajo)

- **Sin mitigar:** +76,5 ± 1,9 mHa, peor que CCSD.
- **Corrigiendo la lectura:** calibramos cada qubit midiendo todo en |0⟩ y todo en |1⟩ (errores de 0,4 % a 0,8 %) y corregimos los valores de cada grupo con la inversa de la matriz de confusión. Queda +47,1 ± 2,0 mHa, que es el sesgo de las compuertas: el valor esperado exacto con ruido y sin errores de lectura da 49 mHa.
- **ZNE global (λ = 1, 3, 5) + lectura:** −2,9 ± 2,1 mHa, compatible con FCI (a 1,4σ) y muy lejos de CCSD (−30,6 mHa). Usa 12 veces más shots porque la extrapolación amplifica σ.
- **Qué no funcionó:** la ZNE sola no alcanza. En una prueba aparte, ZNE (1, 3) sin corregir la lectura dejó +32 mHa, porque plegar el circuito no amplifica el error de lectura. Hay que corregir las dos cosas.

**Limitaciones:** la simulación no incluye el ruido de los qubits que esperan sin compuertas (no agregamos *delays*) ni la diafonía entre qubits. La corrección de lectura supone errores independientes por qubit. Con σ ≈ 2 mHa no alcanza para afirmar precisión química en este punto.
