### Menos qubits: de 8 a 6 y a 5 con simetrías

El hamiltoniano conmuta con tres simetrías de tipo Z (`h4vqe/reduction.py`):

- la paridad de N_α (Z₀Z₁Z₂Z₃) y la de N_β (Z₄Z₅Z₆Z₇);
- la **inversión** de la cadena: cada orbital molecular es par o impar al dar vuelta la cadena, y se conserva el producto de Z sobre los orbitales impares (α y β).

Con el *tapering* de Bravyi *et al.* (2017), cada simetría se rota a una X sobre un qubit, y ese qubit se reemplaza por su valor en el sector del estado HF. Hay dos niveles (`tapering(m, n_qubits=6 o 5)`):

- **6 qubits:** solo las dos paridades. Es lo mismo que el mapeo de paridad con la reducción de 2 qubits, y es lo que pide el extra. Toda excitación conserva N_α y N_β, así que **los 26 operadores del pool de UCCSD/ADAPT se reducen** y el ansatz de excitaciones sigue funcionando igual. Tampoco depende de los orbitales: anda con orbitales rotados o localizados.
- **5 qubits:** además la inversión. Elimina un qubit más, pero tiene dos condiciones: orbitales canónicos (con paridad definida) y solo se reducen las excitaciones simétricas ante la inversión (14 de 26). Es la versión del ansatz compacto.

**Cuánto mejora** (2,4 Å; la celda de abajo lo recalcula):

| qubits | simetrías | términos de Pauli | grupos QWC | shots para la misma σ | excitaciones del pool | orbitales rotados | ansatz y CZ en línea |
|---|---|---|---|---|---|---|---|
| 8 | — | 185 | 37 | 1 | 26/26 | sí | UCCSD ×2: 1648 · ADAPT (sección 3): 378 |
| 6 | N_α, N_β | 165 | 37 | ×1,02 | 26/26 | sí | compacto L = 10: 25 · ADAPT reducido: lo mide la sección 3 |
| 5 | + inversión | 132 | 35 | ×1,13 | 14/26 | **no** | **compacto L = 8: 16** |

Los shots salen de (Σ_g √Var_g)² en el estado exacto, con los grupos y el reparto de Neyman de `estimator.py` (sorted insertion).

**Validación** (primera línea de la salida, en las 18 distancias): el fundamental de Ĥ₆ y de Ĥ₅ coincide con E_FCI, y ⟨HF|Ĥ|HF⟩ con E_HF, con error < 10⁻¹⁴ Ha. `tests/test_reduction.py` chequea además que un ansatz con los 26 operadores reducidos a 6 qubits da la misma energía que en 8 (diferencia de 10⁻¹⁶ Ha).

**Ansatz compacto.** En 5 qubits usamos un ansatz en ladrillo nativo de la línea: X sobre HF, una capa de RY y L capas de [CZ entre pares vecinos alternados + RY]. Con L = 8 tiene 16 CZ, ninguna SWAP con qubits en línea, y llega a FCI sin ruido (tabla de abajo, a 2,0 Å). UCCSD ×2 necesita 1648 CZ en línea: son **100 veces menos CZ** con la misma exactitud. Es el circuito que usamos en el extra de FakeFez. El mismo ansatz en 6 qubits necesita L = 10 (25 CZ); con L = 8 (20 CZ) ya no llega exacto (0,1 mHa a 0,7 Å). La inversión ahorra 9 CZ y 21 parámetros.

**Qué no funcionó y qué cuesta**

- Con L = 6 aparecen mínimos locales (hasta 3,8 mHa de error). Con L = 7 ya llega a FCI.
- Si cada capa pone todas las CZ a la vez (0-1, 1-2, 2-3, 3-4), el ansatz se traba en una fidelidad de 0,961 que no mejora con más capas. Con el patrón en ladrillo llega a 1.
- La reducción no ahorra shots. Hay menos términos, pero casi los mismos grupos, y los términos que quedan tienen más varianza. En 5 qubits hacen falta 1,13 veces más shots para la misma σ a 2,4 Å, y 1,20 veces a 2,0 Å. Con el agrupamiento por coloreo que usábamos antes eran 2 veces: el *sorted insertion* de la sección 4 achicó esa diferencia. Sin ruido, eso solo cuesta tiempo de simulación; con ruido, el circuito de 8 qubits no sirve con ningún número de shots.
- Los orbitales localizados (enlaces 1-2 y 3-4) no tienen paridad ante la inversión. Con ellos, la reducción a 5 qubits no existe y hay que quedarse en 6.
- El ansatz en ladrillo no conserva N ni S² por construcción. Por eso la búsqueda de θ* penaliza ⟨S²⟩, y hay que medir N y S². Además, a diferencia del *tapering*, no escala a moléculas grandes.
