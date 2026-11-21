## 9 · Ruido de hardware real: FakeFez
**Pregunta.** Nuestro modelo de ruido es simple y controlado. ¿La receta (5 qubits, 16 CZ, lectura + ZNE) sobrevive al ruido de un chip real, con errores distintos en cada qubit, relajación ($T_1$, $T_2$) y lectura asimétrica?

`FakeFez` es una copia simulada del procesador IBM Fez (Heron r2, 156 qubits) con su calibración real: el error de cada compuerta, $T_1$/$T_2$ y el error de lectura de cada qubit. Corrimos un punto estirado (2,4 Å) con el ansatz compacto de 5 qubits y 16 CZ (§7.3).

**Cómo lo ubicamos en el chip.**

- El compacto solo usa CZ entre vecinos de una línea, así que entra en el chip **sin SWAP**. Probamos las 856 líneas de 5 qubits de Fez y elegimos la de menos errores esperados, $\sum_{\rm compuertas}\epsilon_g + \sum_{\rm qubits}\epsilon^{\rm lectura}_q$ con la calibración del backend: la mitad que la mediana.
- UCCSD ×2 ubicado en Fez necesita ~1940 CZ por los SWAP. Con los errores de CZ de Fez, a ese circuito no le queda señal.
- Simulamos con el modelo de ruido de Fez restringido a esos 5 qubits, con matriz densidad. Da el mismo valor esperado que el chip completo de 156 qubits (`tests/test_hardware.py`) y es unas 20 veces más rápido.

**Mitigación**, con las energías medidas con shots (tabla de abajo):

- **Sin mitigar:** unos +76 mHa, peor que CCSD.
- **Corrigiendo la lectura** (errores de 0,4 % a 0,8 % por qubit en esa línea, calibrados con 500 k shots por estado): unos +49 mHa, que es el sesgo de las compuertas. El valor esperado exacto con el ruido de Fez y sin error de lectura da 49 mHa: la corrección hace exactamente lo que tiene que hacer.
- **ZNE global ($\lambda$ = 1, 3, 5, Richardson) + lectura:** **compatible con FCI** dentro de 1σ ($\sigma \approx$ 1,2 mHa) y muy lejos de CCSD (−30,6 mHa). Usa 12 veces más shots, porque la extrapolación amplifica la $\sigma$.

**Qué no funcionó:** la ZNE sola no alcanza. En una prueba aparte, ZNE (1, 3) sin corregir la lectura dejó +32 mHa, porque plegar el circuito no amplifica el error de lectura. Hay que corregir las dos cosas, igual que con nuestro modelo (§8).

**Limitaciones:** la simulación no incluye el ruido de los qubits que esperan sin compuertas (no agregamos *delays*) ni la diafonía entre qubits, y la corrección de lectura supone errores independientes por qubit. Con $\sigma \approx$ 1,2 mHa el valor central queda dentro de la precisión química, pero la barra no alcanza para afirmarlo con confianza: harían falta ~6 veces más shots.
