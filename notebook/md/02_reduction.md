### Menos qubits: de 8 a 5 con simetrías

El hamiltoniano conmuta con tres simetrías de tipo Z (`h4vqe/reduction.py`):

- la paridad de N_α (Z₀Z₁Z₂Z₃) y la de N_β (Z₄Z₅Z₆Z₇);
- la **inversión** de la cadena: cada orbital molecular es par o impar al dar vuelta la cadena, y se conserva el producto de Z sobre los orbitales impares (α y β).

Con el *tapering* de Bravyi *et al.* (2017), cada simetría se rota a una X sobre un qubit, y ese qubit se reemplaza por su valor en el sector del estado HF. Se eliminan 3 qubits: **de 8 a 5**. Esto va un paso más allá del extra pedido (6 qubits): el mapeo de paridad solo elimina 2, y la inversión elimina uno más.

**Validación** (primera línea de la salida, en las 18 distancias): el fundamental de Ĥ₅ coincide con E_FCI y ⟨HF|Ĥ₅|HF⟩ con E_HF, con error < 10⁻¹⁴ Ha. Ĥ₅ tiene 132 términos de Pauli (contra 185) y 34 grupos de medición QWC (contra 46).

**Ansatz compacto.** En 5 qubits usamos un ansatz en ladrillo nativo de la línea: X sobre HF, una capa de RY y L capas de [CZ entre pares vecinos alternados + RY]. Con L = 8 tiene 16 CZ, ninguna SWAP con qubits en línea, y llega a FCI sin ruido (tabla de abajo, a 2,0 Å). UCCSD ×2 necesita 1648 CZ en línea: son **100 veces menos CZ** con la misma exactitud. Es el circuito que usamos en el extra de FakeFez.

**Qué no funcionó y qué cuesta**

- Con L = 6 aparecen mínimos locales (hasta 3,8 mHa de error). Con L = 7 ya llega a FCI.
- Si cada capa pone todas las CZ a la vez (0-1, 1-2, 2-3, 3-4), el ansatz se traba en una fidelidad de 0,961 que no mejora con más capas. Con el patrón en ladrillo llega a 1.
- Menos grupos no significa menos shots: para la misma σ, en 5 qubits hacen falta unas 2 veces más shots, porque los términos que quedan tienen más varianza. Sin ruido eso solo cuesta tiempo de simulación; con ruido, el circuito de 8 qubits no sirve con ningún número de shots.
- El ansatz en ladrillo no conserva N ni S² por construcción. Por eso la búsqueda de θ* penaliza ⟨S²⟩, y hay que medir N y S². Además, a diferencia del *tapering*, no escala a moléculas grandes.
