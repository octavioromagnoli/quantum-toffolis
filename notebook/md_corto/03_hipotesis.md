### 2.5 Qué dio lo mínimo

- **UCCSD ×1** (292 CZ) le gana a CCSD, pero **no llega a precisión química** (4,5 mHa a 2,0 Å) y **rompe el espín** ($\langle S^2\rangle$ = 0,85 a 3,0 Å).
- **UCCSD ×2** (582 CZ) llega a FCI. Con 8 M de shots por punto ya cumplía la consigna: 7/7 contra CCSD y precisión química en los 7 puntos estirados, $R_e$ = 0,8932 ± 0,0018 Å, $D_e$ = 312,95 ± 0,83 mHa.

**La consigna mínima ya estaba cumplida.** Pero quedaron cinco problemas, uno por capítulo:

| # | problema | qué hicimos | § |
|---|---|---|---|
| 1 | **582 CZ** (1648 en línea): con ruido, casi ningún shot sale limpio | ADAPT-VQE y ruteo fSWAP | §3 |
| 2 | **8 M de shots** por punto y una σ dudosa (dispersión de $z$ = 1,16 en vez de 1) | mejor agrupamiento y 3600 repeticiones para calibrar σ | §4 |
| 3 | $R_e$ con un ajuste que tenía **sesgo** | elegir el ajuste por su sesgo contra FCI | §5 |
| 4 | la confiabilidad se argumentaba con el **estado del simulador** | medir varianza y simetrías con shots | §6 |
| 5 | con ruido, UCCSD ×2 es **peor que Hartree-Fock** | 5 qubits, 16 CZ y mitigación | §7–§9 |
