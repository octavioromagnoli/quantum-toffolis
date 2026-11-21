## 2 · El ansatz  *(Gonzalo)*

UCCSD con Trotter: 8 simples + 18 dobles = 26 excitaciones por paso. Cada excitación es exp(−iθG) con G = i(T − T†); sus términos de Pauli conmutan, así que el circuito es exactamente Π exp(−iθ_k G_k). Eso nos deja optimizar con un simulador numpy (gradiente analítico) y medir **el mismo circuito**.

| ansatz | parámetros | CZ (todos con todos) | CZ (línea) | error máx. vs FCI |
|---|---|---|---|---|
| UCCSD ×1 | 26 | 292 | 806 | 4,5 mHa (2,0 Å) |
| UCCSD ×2 | 52 | ~584 | ~1600 | < 0,001 mHa |

> TODO (Gonzalo): bajar CZ manteniendo < 1,6 mHa: ADAPT-VQE (pool fermiónico o de excitaciones de qubits), compuertas de Givens / double-excitation eficientes, poda de excitaciones con θ≈0.
