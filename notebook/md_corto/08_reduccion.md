**Qué dio.** Con ruido, UCCSD ×2 queda **peor que Hartree-Fock**, que no usa ningún circuito. ADAPT mejora, pero sigue a cientos de mHa. Solo el compacto de 16 CZ deja un sesgo de pocos mHa, que la mitigación (§8) puede sacar.

### 7.2 Menos qubits con simetrías (8 → 6 → 5)

$\hat H$ conserva la paridad de $N_\alpha$, la de $N_\beta$ y la paridad ante la inversión de la cadena. Cada simetría con valor conocido permite **eliminar un qubit** (*tapering*, Bravyi *et al.*, 2017) sin cambiar la energía.

- **6 qubits** (las dos paridades): lo que pide el extra. El ansatz químico funciona igual.
- **5 qubits** (además la inversión): un qubit menos.

### 7.3 Un ansatz nativo de la línea: 16 CZ

En 5 qubits usamos un ansatz **en ladrillo**: capas de $R_Y$ en cada qubit y CZ entre vecinos, alternando pares. Con 8 capas: **16 CZ, sin SWAP** y 45 parámetros.

No conserva $N$ ni $S^2$ por sí solo, así que la búsqueda penaliza $\langle S^2\rangle$. **Sin esa penalización, una corrida terminó en el triplete**, con un error de apenas 0,95 mHa pero el estado equivocado.

La celda de abajo valida la reducción y compara el compacto con UCCSD ×2.
