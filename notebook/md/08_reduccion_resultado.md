**Qué dio la reducción.**

- **Validación** (primera línea de la salida): en las 18 distancias, el fundamental de $\hat H_6$ y $\hat H_5$ coincide con $E_{\rm FCI}$, y $\langle \mathrm{HF}\rvert\hat H\lvert \mathrm{HF}\rangle$ con $E_{\rm HF}$, con error < $10^{-13}$ Ha. `tests/test_reduction.py` chequea además que un ansatz con los 26 operadores reducidos a 6 qubits da la misma energía que en 8 (diferencia de $10^{-16}$ Ha).
- **Medición:** de 8 a 5 qubits los términos de Pauli bajan de 185 a 132, pero los grupos QWC casi no cambian (37 → 35 a 2,4 Å), y los términos que quedan tienen más varianza. Para la misma $\sigma$ hacen falta ~1,1 veces **más** shots en 5 qubits. La reducción no ahorra shots: **ahorra compuertas**.
- **Ansatz:** el compacto de 5 qubits llega a FCI sin ruido (error < $10^{-9}$ mHa a 2,0 Å) con **16 CZ, 100 veces menos** que las 1648 de UCCSD ×2 en línea. El mismo ansatz en 6 qubits necesita $L$ = 10 (25 CZ). La inversión ahorra 9 CZ y 21 parámetros.

**Por qué funcionó.** Las simetrías son información que ya conocemos (cuántos electrones $\alpha$ y $\beta$ hay, y la paridad ante la inversión), y guardarla en qubits es un desperdicio. Al sacarla, el espacio que el ansatz tiene que recorrer baja de 256 estados a 32, y un circuito genérico en ladrillo alcanza para llegar al fundamental. La geometría en línea de ese circuito coincide con la del chip, así que no hace falta ninguna SWAP.

**Qué no funcionó y qué cuesta**

- Con $L$ = 6 capas aparecen mínimos locales (hasta 3,8 mHa de error). Con $L$ = 7 (14 CZ) ya llega a FCI; $L$ = 8 da margen.
- Si cada capa pone **todas** las CZ a la vez (0-1, 1-2, 2-3, 3-4), el ansatz se traba en una fidelidad de 0,961, que no mejora con más capas. Con el patrón en ladrillo llega a 1.
- Los orbitales localizados (enlaces 1-2 y 3-4) no tienen paridad definida ante la inversión: con ellos solo se puede bajar a 6 qubits.
- El ansatz en ladrillo **no conserva $N$ ni $S^2$** por construcción: por eso la penalización de $S^2$ y la medición de $N$ y $S^2$. A diferencia del tapering, **no escala** a moléculas grandes: es una solución para H₄, no un método general.
