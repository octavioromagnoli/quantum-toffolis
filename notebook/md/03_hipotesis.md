### 2.5 Qué dio lo mínimo y qué nos dejó pensando

**UCCSD ×1** (26 parámetros, 292 CZ) ya le gana a CCSD en toda la zona estirada, pero **no llega a precisión química** (hasta 4,5 mHa a 2,0 Å) y **rompe el espín**: a 3,0 Å $\langle S^2\rangle$ llega a 0,85 (0 es singulete y 2 es triplete), porque converge a una mezcla de los dos. Duplicar los pasos de Trotter (**UCCSD ×2**, 52 parámetros) llega a FCI con error < 0,001 mHa en toda la curva. Esa fue nuestra primera entrega completa, y en esa corrida (11:09; 18 distancias, 8 M de shots por punto):

- le ganaba a CCSD en 7 de 7 puntos estirados, con precisión química en los 7 (error máximo medido: 0,74 mHa);
- $\sigma$ = 0,65–0,78 mHa por punto;
- $R_e$ = 0,8932 ± 0,0018 Å con un ajuste cuadrático de 5 puntos, y $D_e$ = 312,95 ± 0,83 mHa;
- la varianza, las simetrías y la comparación con UCCSD ×1 se calculaban con el estado del simulador.

O sea: **la consigna mínima ya estaba cumplida**. Pero al mirar los números aparecieron cinco problemas, y cada uno se convirtió en una hipótesis y en un capítulo de este notebook:

| # | qué observamos | hipótesis | qué probamos | capítulo |
|---|---|---|---|---|
| 1 | UCCSD ×2 usa **582 CZ** (1648 con qubits en línea). Con 0,1 % de error por CZ, menos de 1 de cada 5 shots sale sin error. | No todas las excitaciones aportan lo mismo: alcanza con las que más bajan la energía. | ADAPT-VQE (elige por gradiente) y un ruteo propio con fSWAP para la línea | §3 |
| 2 | Hacen falta **8 M de shots** por punto para $\sigma \approx$ 0,7 mHa, y en 60 repeticiones la dispersión de $z = (E - E_{\rm exacta})/\sigma$ dio **1,16** en lugar de 1. | La varianza depende de cómo se agrupan y reparten los shots. Y el 1,16 puede ser fluctuación estadística de solo 60 repeticiones. | *Sorted insertion*, Neyman en dos etapas, conmutación general y 3600 repeticiones con bootstrap | §4 |
| 3 | $R_e$ salía de un ajuste cuadrático (0,8932 ± 0,0018 Å) y lo comparábamos con FCI **ajustado con el mismo modelo** (0,8964 Å). Esa comparación esconde el error del ajuste: el mínimo FCI verdadero está en 0,8886 Å. | El modelo de ajuste tiene un sesgo propio, porque la curva es asimétrica alrededor del mínimo. | Aplicar cada modelo a FCI en la misma grilla y compararlo con el mínimo FCI exacto | §5 |
| 4 | La confiabilidad se argumentaba con el **estado del simulador**, que en hardware no existe. | La varianza y las simetrías se pueden **medir**, y un eco de Loschmidt resuelve varianzas de $10^{-6}$ Ha². | Varianza por grupos, eco, $N_\alpha$, $N_\beta$, $S^2$ medidos y convergencia con la profundidad | §6 |
| 5 | Con el modelo de ruido, UCCSD ×2 queda **peor que Hartree-Fock**. | Hacen falta unas 100 veces menos CZ. En 8 qubits un ansatz químico no llega; sacando qubits con simetrías, uno nativo de la línea sí. | Tapering 8 → 6 → 5 qubits, ansatz de 16 CZ, mitigación y FakeFez | §7–§9 |

Los cinco se atacaron en paralelo entre las 11:10 y las 14:40, cada persona en su rama. La versión final del ansatz de la curva (ADAPT con 15 operadores) quedó integrada a las 13:41.
