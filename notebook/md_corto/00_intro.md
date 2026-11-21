# Romper una molécula, un qubit a la vez: H₄ con VQE

**Quantum Jam 2026 · equipo quantum-toffolis** (Octavio, Gonzalo, Santino, Joaquín)

### Resumen

Calculamos con VQE la curva de energía de la cadena lineal H₄ (STO-3G, 8 qubits) en 18 distancias entre 0,6 y 3,0 Å. Todas las energías se **miden** con `SamplerV2` de `qiskit_aer`.

- **Ansatz:** ADAPT-VQE con 15 excitaciones, **194 CZ** (378 en línea). El UCCSD "de libro" usa 582 (1648).
- **Resultado:** donde CCSD se equivoca entre 9 y 31 mHa, el VQE queda a ≤ 0,56 mHa de la energía exacta. Le gana a CCSD en 7/7 puntos estirados y tiene **precisión química en las 18 distancias**.
- **$R_e$** = 0,8873 ± 0,0018 Å (FCI: 0,8886) y **$D_e$** = 313,00 ± 0,66 mHa (FCI: 313,21).
- **Confiabilidad sin mirar FCI:** varianza de la energía medida, simetrías medidas y convergencia con la profundidad.
- **Extras:** con ruido de 0,1 % por CZ, bajando a 5 qubits y mitigando errores, también llegamos a precisión química (16 CZ). Además: un chip de IBM simulado (FakeFez), otras geometrías, comparación de ansätze y escalado hasta FeMoco.

El notebook cuenta **cómo llegamos**: arranca del VQE mínimo que pide la consigna y sigue los problemas en el orden en que aparecieron. Cada capítulo: **qué vimos → qué hicimos → qué dio**.

### Dónde está cada pedido de la consigna

| pedido | dónde | resultado |
|---|---|---|
| Curva $E(d)$, ≥ 10 distancias, con $\pm\sigma$, contra CCSD y FCI | §5 | 18 distancias, σ = 0,50–0,60 mHa |
| Energías medidas con `SamplerV2` | §2, §4, §5 | todas |
| $R_e$ y $D_e$ | §5 | ver arriba |
| Ganarle a CCSD en $d \ge 1{,}8$ Å / precisión química | §5 | 7/7 y 7/7 |
| Confiabilidad antes de comparar con FCI | §6 | varianza, simetrías, profundidad |
| Decisiones y qué no funcionó | cada capítulo y §13 | |
| Extras | §7–§12 | ruido, mitigación, FakeFez, geometrías, ansätze, FeMoco |

> **Para correrlo:** *Entorno de ejecución → Ejecutar todo*. Tarda varios minutos (más en Colab gratis). Las energías reportadas siempre salen de shots; el estado exacto solo se usa para buscar parámetros y comparar alternativas.
