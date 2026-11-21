# Romper una molécula, un qubit a la vez: H₄ con VQE

**Quantum Jam 2026 · equipo quantum-toffolis** (Octavio, Gonzalo, Santino, Joaquín)

> **En una frase.** Donde CCSD se equivoca hasta 31 mHa, nuestro VQE mide la energía de la cadena H₄ estirada con un error de 0,56 mHa como máximo (precisión química en las 18 distancias) y con 3 veces menos compuertas CZ que el VQE "de libro". Con ruido, bajando de 8 a 5 qubits, llega a precisión química con 16 CZ y mitigación de errores.

Este notebook cuenta **cómo llegamos**: arranca de lo mínimo que pide la consigna (un VQE con UCCSD que ya le ganaba a CCSD, pero con 1648 CZ en línea) y sigue el orden en que fuimos encontrando problemas. Cada capítulo tiene la misma estructura: **qué observamos → qué hipótesis planteamos → qué hicimos → qué dio y por qué funcionó (o no)**. Las horas son las de los commits del 21/11/2026.

### Qué pide la consigna y dónde está

| pedido de la consigna | dónde | resultado |
|---|---|---|
| Curva $E(d)$ con ≥ 10 distancias entre 0,6 y 3,0 Å, con $\pm\sigma$, contra CCSD y FCI | §5 | 18 distancias, $\sigma$ = 0,50–0,60 mHa |
| Energías medidas con `SamplerV2` de `qiskit_aer` (no leyendo el estado) | §2.3, §4, §5 | todas las energías reportadas |
| $R_e$ y $D_e = E(3{,}0\ \text{Å}) - E(R_e)$ | §5 | $R_e$ = 0,8873 ± 0,0018 Å (FCI: 0,8886); $D_e$ = 313,00 ± 0,66 mHa (FCI: 313,21) |
| Ganarle a CCSD para $d \ge 1{,}8$ Å contando la barra de error | §5 | 7 de 7 puntos |
| Meta máxima: precisión química (< 1,6 mHa) en esa zona | §5 | 7 de 7 puntos (error máximo: 0,56 mHa) |
| Argumentar la confiabilidad **antes** de comparar con FCI | §6 | varianza medida (eco de Loschmidt), simetrías, profundidad y límite físico |
| Tabla con energía, error, shots y CZ por punto | §5 | columnas `E_vqe`, `sigma_mHa`, `shots`, `cz_all`, `cz_line_fswap` |
| Qué decidimos, por qué y qué no funcionó | cada capítulo y §13 | |
| Extras (opcionales) | §7–§12 | ruido, mitigación, FakeFez, 5 qubits, precisión química con ruido, comparar ansätze, otra geometría y escalado a FeMoco |

### El recorrido

1. **El problema** (§1): el hamiltoniano de H₄ escrito en qubits, término por término, y por qué CCSD falla al estirar.
2. **Lo mínimo** (§2, 11:09): VQE con UCCSD. Cumple la consigna, pero deja **cinco problemas** a la vista, y cada uno se vuelve una hipótesis.
3. **Problema 1, demasiadas CZ** (§3): ADAPT-VQE elige solo las excitaciones que hacen falta, y un ruteo propio con fSWAP baja las CZ con qubits en línea: de 1648 a 378.
4. **Problema 2, demasiados shots y una σ dudosa** (§4): mejor agrupamiento y reparto de Neyman (hasta 1,9× menos shots) y 3600 repeticiones para calibrar la σ.
5. **Resultado central** (§5): la curva, la tabla, $R_e$ y $D_e$. **Problema 3**: el ajuste del mínimo tenía un sesgo de +8 mÅ.
6. **Problema 4, confiabilidad sin FCI** (§6): todo medido con shots, incluida la varianza por eco de Loschmidt.
7. **Problema 5, el ruido** (§7–§9): con 0,1 % de error por CZ, el VQE "de libro" es peor que no hacer nada. Simetrías (8 → 5 qubits), un ansatz de 16 CZ, mitigación y una copia de un chip de IBM.
8. **Más allá de la consigna** (§10–§12): otras geometrías, comparación de ansätze y cuánto costaría la molécula real (FeMoco).
9. **Resumen final** (§13): lo mínimo contra lo final, número por número, y las decisiones que explican la mejora.

> **Para correrlo:** *Entorno de ejecución → Ejecutar todo*. Instala sus dependencias, escribe el paquete `h4vqe/` y recalcula todo. Es largo: el resumen final (§13) dice cuánto tardó (en Colab gratis tarda bastante más que en una PC). Las energías reportadas salen siempre de medir circuitos con shots. Los cálculos con el estado exacto (sin shots) solo se usan para buscar parámetros y elegir entre alternativas, como permite la consigna, y están marcados como tales.
