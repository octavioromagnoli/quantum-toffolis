## 0 · Instalación y código

Las dos celdas siguientes están colapsadas (doble clic para verlas): instalan las dependencias y escriben el paquete `h4vqe/`, así el notebook corre solo en un Colab nuevo.

| módulo | qué hace | capítulo |
|---|---|---|
| `hamiltonian` | PySCF → Ĥ en qubits; HF, FCI y CCSD | §1 |
| `ansatz`, `optimize` | circuitos de excitaciones y búsqueda de θ | §2, §3 |
| `adapt`, `routing` | ADAPT-VQE y ruteo fSWAP en línea | §3 |
| `estimator` | energía **medida** con `SamplerV2`, con su σ | §2, §4 |
| `analysis`, `pipeline` | curva, $R_e$, $D_e$ y diagnósticos | §5, §6 |
| `reduction`, `noise`, `mitigation` | 5 qubits, modelo de ruido y mitigación | §7, §8 |
| `hardware`, `geometry`, `scaling` | FakeFez, otras geometrías, FeMoco | §9–§12 |
| `historia` | diagnósticos de este relato | §1–§7, §13 |
