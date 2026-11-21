## 1 · De la molécula a los qubits  *(Octavio)*

- Geometría: 4 H en línea sobre z, separados por d.
- RHF → integrales en la base de orbitales moleculares → Ĥ con operadores fermiónicos → Jordan-Wigner.
- Convención: qubits 0–3 = espín-orbitales α, 4–7 = β; HF = `|00110011⟩`.
- Chequeos: ⟨HF|Ĥ|HF⟩ = E_HF de PySCF y referencia FCI (CASCI completo).

> TODO (Octavio): explicar por qué CCSD falla al estirar (es no variacional: queda *por debajo* de FCI) y, como extra, reducción de qubits (paridad + simetrías → 6 qubits).
