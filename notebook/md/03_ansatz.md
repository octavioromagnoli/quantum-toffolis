## 2 · El ansatz  *(Gonzalo)*

Partimos de UCCSD con Trotter: 8 simples + 18 dobles = 26 excitaciones por paso. Cada excitación es exp(−iθG) con G = i(T − T†). Sus términos de Pauli conmutan, así que el circuito es exactamente Π exp(−iθ_k G_k). Eso nos deja optimizar con un simulador numpy (gradiente analítico) y medir **el mismo circuito**.

El problema: UCCSD×2 llega a FCI, pero cuesta ~584 CZ (~1600 en línea). Con 0,1 % de error por CZ eso no sirve.

### ADAPT-VQE: solo las excitaciones que hacen falta

ADAPT (Grimsley et al., 2019) arranca de |HF⟩ y agrega de a una la excitación del pool con gradiente máximo, g_k = ⟨ψ|i[G_k, H]|ψ⟩. Después de cada agregado reoptimiza todos los parámetros. Pool: las 26 excitaciones de UCCSD, que se pueden repetir (`h4vqe/adapt.py`).

**Una lista común para toda la curva.** En cada paso elegimos el operador con mayor gradiente en *cualquiera* de las 18 distancias y reoptimizamos cada una por separado. Así el circuito (y la cantidad de CZ) es el mismo en todos los puntos, y solo cambian los θ.

**Dónde cortar sin usar FCI.** Medimos la varianza de la energía ⟨H²⟩ − ⟨H⟩², que vale 0 solo en un autoestado y se puede medir en una computadora cuántica. Al pasar de 12 a 13 operadores, la varianza máxima sobre la curva cae de 3,2·10⁻³ a 2,6·10⁻⁴ Ha², y el error máximo cae de 3,7 a 0,6 mHa. Cortamos en **13 operadores** (`ADAPT_FERMIONIC`).

| ansatz | parámetros | CZ (todos con todos) | CZ (línea) | error máx. vs FCI (sin ruido) |
|---|---|---|---|---|
| UCCSD ×1 | 26 | 292 | 806 | 4,5 mHa (2,0 Å) |
| UCCSD ×2 | 52 | ~584 | ~1600 | < 0,001 mHa |
| **ADAPT fermiónico, 13 ops (default nuevo)** | **13** | **185** | **~515** | **0,84 mHa (1,4 Å); ≤ 0,6 mHa para d ≥ 1,8 Å** |
| ADAPT-QEB, 15 ops (extra) | 15 | 159 | ~451 | 0,81 mHa con la lista; ver abajo |

El pool QEB usa excitaciones de *qubit* (sin cadenas Z de Jordan-Wigner). Conservan N_α y N_β, tocan solo 2 o 4 qubits y cuestan menos CZ.

### Qué no funcionó

- **Umbral de gradiente 1e-4:** ADAPT se frena en mesetas a 1,6, 2,8 y 3,0 Å con 0,6-2 mHa de error, y después el gradiente vuelve a crecer. Con 1e-5 llega a FCI en toda la curva.
- **QEB con warm start entre distancias:** a 2,0 Å cae en un mínimo local malo (9,8 mHa, S² = 0,78). Es más barato pero menos robusto; queda como extra.
- **Probar layouts y semillas del transpiler** para conectividad en línea: el mejor resultado fue ~515 CZ, contra ~526 por defecto. El sobrecosto (~2,8× respecto de todos con todos) viene de los SWAP que necesitan las excitaciones αβ.

### Limitaciones y próximos pasos

- S² ≤ 0,013 (no es exactamente 0), porque el pool no está adaptado al espín.
- Para bajar las CZ en línea: elegir operadores por gradiente / CZ en línea (Co-ADAPT, Ramôa et al. 2026), podar operadores de θ chico, operadores CEO de 9 CNOT y orbitales localizados.
