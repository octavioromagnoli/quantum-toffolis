# H₄ con VQE — extras

**Quantum Jam 2026 · equipo quantum-toffolis** (Octavio, Gonzalo, Santino, Joaquín)

Este notebook complementa a `H4_VQE.ipynb`, que tiene lo que pide la consigna: la curva con energías medidas, R_e, D_e y la confiabilidad. Acá van los extras. Cada uno corre solo, sin depender de la curva principal.

1. **Menos qubits** *(Octavio)*: tapering con simetrías Z₂ de 8 a 6 y a 5 qubits, y un ansatz compacto de 5 qubits con 16 CZ que llega a FCI sin ruido.
2. **Mitigación de errores** *(Santino)*: corrección de lectura, post-selección por simetría y ZNE, sobre el modelo de ruido propio, con el ansatz compacto.
3. **Ruido de hardware real** *(Octavio)*: un punto estirado en `FakeFez` (IBM Heron r2), con corrección de lectura y ZNE.
4. **Otra geometría** *(Octavio)*: H₄ cuadrado y en zigzag, con el mismo flujo que la cadena recta.

La comparación de ansätze (UCCSD, ADAPT fermiónico y QEB, ruteo fSWAP) está en la sección 2 del notebook principal.
