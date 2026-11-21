## 7 · Problema 5, el ruido: las CZ deciden

**Desde acá, extras.** El desafío central ya está cumplido; ahora vemos si el método sobrevive al ruido.

**Modelo de ruido** (el de la consigna, implementado en `noise.py`): qubits en línea, error de 0,1 % por CZ y 0,005 % por compuerta de un qubit. La probabilidad de que el circuito corra sin errores es $\approx e^{-n_{\rm CZ}/1000}$.

**Qué vimos.** UCCSD ×2 en línea tiene ~1550 CZ: solo ~20 % de los shots salen limpios. ADAPT con ruteo baja a 313, todavía demasiado.

**Idea.** Hacen falta **~100 veces menos CZ**. Un ansatz químico en 8 qubits no llega; uno nativo de la línea, en menos qubits, sí.

La celda de abajo calcula, a 2,4 Å, el sesgo **exacto** que deja el ruido en tres ansätze (sin shots: ese sesgo no baja con más shots).
