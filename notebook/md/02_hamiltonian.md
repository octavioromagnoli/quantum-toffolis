## 1 · De la molécula a los qubits  *(Octavio)*

- Geometría: 4 H en línea sobre z, separados por d.
- RHF → integrales en la base de orbitales moleculares → Ĥ con operadores fermiónicos → Jordan-Wigner (8 qubits, 185 términos de Pauli).
- Convención: qubits 0–3 = espín-orbitales α, 4–7 = β; HF = `|00110011⟩`.

### Validación del hamiltoniano

En las 18 distancias de la curva comparamos Ĥ, ya escrito en qubits, contra PySCF:

- ⟨HF|Ĥ|HF⟩ coincide con E_HF, con una diferencia máxima de 7·10⁻¹⁵ Ha.
- El autovalor más bajo de Ĥ **restringido al sector N_α = N_β = 2** (36 de los 256 estados de la base) coincide con E_FCI (CASCI completo), con la misma diferencia máxima, y su estado es un singulete (S² = 0). Restringimos al sector porque Ĥ en qubits también describe estados con otra cantidad de electrones (H₄⁺, H₄⁻, …), y el que nos interesa es el de 4 electrones con S_z = 0.

Con esto, cualquier diferencia entre VQE y FCI viene del ansatz, de la optimización o de las mediciones, no del hamiltoniano.

### Por qué CCSD falla al estirar

El gráfico de abajo muestra el error de CCSD contra FCI (izquierda) y el peso del determinante HF en el estado exacto (derecha). La zona gris es la zona estirada (d ≥ 1,8 Å).

- **Cerca del equilibrio** (d ≈ 0,9 Å) el estado exacto es casi el determinante HF (peso 0,95) y CCSD es prácticamente exacto (+0,007 mHa).
- **Al estirar**, cada electrón se queda en su átomo y el estado exacto pasa a ser una mezcla de muchos determinantes con pesos parecidos: el peso de HF baja a 0,58 a 1,8 Å y a 0,27 a 3,0 Å. Esto es **correlación estática**. CCSD parte de un solo determinante y solo incluye las excitaciones simples y dobles conectadas. Las triples y cuádruples aparecen solo como productos de las dobles (para 4 electrones, FCI equivale a CCSDTQ), y eso deja de alcanzar cuando HF ya no domina.
- **CCSD no es variacional.** Su energía sale de proyectar la ecuación de Schrödinger sobre HF y sus excitaciones simples y dobles, no del valor esperado ⟨ψ|Ĥ|ψ⟩ de un estado. Por eso no es una cota superior y puede quedar por debajo de FCI, que es lo que pasa acá: a partir de 1,1 Å queda por debajo, con −9,1 mHa a 1,8 Å y un máximo de −30,6 mHa entre 2,4 y 2,6 Å, unas 19 veces la precisión química. Quedar por debajo no es "mejor": el error en valor absoluto es grande y el signo muestra que el método salió de su rango de validez. CCSD convergió en todas las distancias, así que no es un problema numérico.
- **VQE sí es variacional:** E(θ) = ⟨ψ(θ)|Ĥ|ψ(θ)⟩ ≥ E₀ para cualquier θ, así que no puede quedar por debajo de FCI salvo por el ruido estadístico de las mediciones. Eso es lo que le permite ganarle a CCSD en la zona estirada.
