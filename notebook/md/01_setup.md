## 0 · Instalación y código

La primera celda instala lo que falte (Qiskit 2.6, Qiskit Aer 0.18, PySCF 2.15 y `qiskit-ibm-runtime` para FakeFez). La segunda escribe el paquete `h4vqe/`, así el notebook corre de punta a punta en un Colab nuevo. Las dos tienen el código oculto: con doble clic (o *Mostrar código* en Colab) se ve. Los módulos:

| módulo | qué hace | capítulo |
|---|---|---|
| `config` | distancias, shots por punto, ansatz de la curva | — |
| `hamiltonian` | PySCF → Ĥ en qubits (Jordan-Wigner); HF, FCI y CCSD de referencia | §1 |
| `ansatz`, `optimize` | excitaciones (UCCSD) como circuito, simulador numpy exacto y búsqueda de θ con gradiente analítico | §2, §3 |
| `adapt`, `routing` | ADAPT-VQE con lista común y ruteo fSWAP para qubits en línea | §3 |
| `estimator` | energía **medida** con `SamplerV2` de Aer: grupos, reparto de shots y σ | §2, §4 |
| `analysis`, `pipeline` | curva con warm start, tabla, $R_e$, $D_e$ y diagnósticos medidos | §5, §6 |
| `reduction`, `noise` | tapering 8 → 6 → 5 qubits, ansatz compacto y modelo de ruido | §7 |
| `mitigation` | corrección de lectura, post-selección por simetría y ZNE | §8 |
| `hardware`, `geometry` | FakeFez y otras geometrías | §9, §10 |
| `scaling` | escalado H₄ → H₁₄ → FeMoco | §12 |
| `historia` | diagnósticos de este relato (lo mínimo, ADAPT, shots, ruido por ansatz) | §1–§7, §13 |
