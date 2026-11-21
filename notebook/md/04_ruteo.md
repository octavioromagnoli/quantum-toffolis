### 3.2 Qué dio ADAPT y cómo lo cerramos

**El corte en 13 operadores.** Los primeros 8 son dobles $\alpha\beta$ (un electrón $\alpha$ y uno $\beta$ que saltan juntos), después vienen las dos dobles del mismo espín y recién ahí las simples. En la zona estirada, la varianza máxima cae **12 veces** al pasar de 12 a 13 operadores (de ~3·10⁻³ a ~3·10⁻⁴ Ha², tabla de arriba), y el error máximo de la curva, de 3,7 a 0,6 mHa. Cerca del equilibrio la varianza ya era chica y casi no cambia (con 13 operadores, el máximo de toda la curva queda a 1,2 Å, con 0,3 mHa de error). La primera lista cortó ahí: **13 operadores y 185 CZ**, contra 582 de UCCSD ×2.

**Completar los pares de espín.** Entre los 13 operadores, ADAPT eligió las simples $0_\alpha \to 2_\alpha$ y $1_\beta \to 3_\beta$, pero no sus pares $0_\beta \to 2_\beta$ y $1_\alpha \to 3_\alpha$. El estado rompía un poco el espín ($\langle S^2\rangle$ hasta 0,013). Agregamos esas dos simples **al lado de su pareja** (`ADAPT_FERMIONIC`, 15 operadores): el error máximo baja de 0,59 a **0,38 mHa** en toda la curva y $\langle S^2\rangle$ a ≤ 6,5·10⁻³, por 9 CZ más (194 en total). Corriendo ADAPT con dos operadores más, el algoritmo elige justamente esas dos simples: la corrección coincide con lo que pide el gradiente.

**Dos arranques en la optimización.** Solo con *warm start* (el $\boldsymbol\theta^*$ de la distancia anterior), el optimizador caía en mínimos malos en algunos puntos: 2 mHa a 2,8 Å con la lista de 15. Ahora `optimize_params` prueba siempre el *warm start* **y** $\boldsymbol\theta = 0$ ($\lvert \mathrm{HF}\rangle$), y se queda con la energía más baja. Como el método es variacional, la más baja es siempre la mejor.

### 3.3 Qubits en línea: ruteo propio con fSWAP
**Observación.** Con qubits en línea (como en el modelo de ruido y en muchos chips) el transpiler de Qiskit agrega SWAP y multiplica las CZ ~2,8×: ADAPT pasa de 194 a ~530.

**Hipótesis.** Una excitación fermiónica solo necesita que sus modos estén **contiguos**. Si movemos los modos con *fermionic SWAP* en lugar de SWAP comunes, el signo de Jordan-Wigner viaja con el modo y no hay que deshacer nada.

$$\mathrm{fSWAP}_{j,j+1} = \begin{pmatrix} 1&0&0&0\\ 0&0&1&0\\ 0&1&0&0\\ 0&0&0&-1 \end{pmatrix} = \mathrm{CZ}\cdot\mathrm{SWAP} .$$

Intercambia dos modos e incluye el signo $-1$ cuando los dos están ocupados, que es justo lo que pide el anticonmutador. Es un iSWAP salvo fases de un qubit, así que cuesta **2 CZ**, contra 3 de un SWAP común.

Nuestro ruteo (`h4vqe/routing.py`, inspirado en las redes fSWAP de Kivlichan *et al.*, 2018, y en Co-ADAPT de Ramôa *et al.*, 2026):

1. Antes de cada excitación junta sus modos en qubits **contiguos** con fSWAP.
2. **No deshace** la permutación: el orden de los modos se sigue de forma clásica.
3. Aplica cada doble con un circuito pensado para la línea (18 CZ, casi todas entre vecinos) y mide la energía con $\hat H$ escrito en el **orden final** de modos.

Los $\boldsymbol\theta$ son los mismos (las fSWAP solo reetiquetan modos), y la energía coincide hasta $10^{-15}$ (`tests/test_routing.py`). El orden inicial se eligió entre 72 candidatos: $[3_\alpha\, 3_\beta\, 1_\alpha\, 1_\beta\, 2_\alpha\, 2_\beta\, 0_\alpha\, 0_\beta]$. Resultado: **530 → 378 CZ en línea** (−29 %). El costo: $\hat H$ en el orden final se agrupa peor (69 grupos QWC contra 37), así que necesita más shots para la misma σ. Por eso la curva central se mide con el circuito sin rutear y $\hat H$ original (es el mismo estado), y el ruteado se usa para contar CZ en línea y para simular con ruido.

| ansatz | parámetros | CZ (todos con todos) | CZ (línea) | error máx. sin ruido | $\langle S^2\rangle$ máx. |
|---|---|---|---|---|---|
| UCCSD ×1 | 26 | 292 | 806 | 4,5 mHa (2,0 Å) | 0,85 |
| UCCSD ×2 (lo mínimo) | 52 | 582 | 1648 | < 0,001 mHa | ≈ 0 |
| ADAPT, 13 operadores | 13 | 185 | 526 (transpiler) / 370 (fSWAP) | 0,59 mHa (2,0 Å) | 0,013 |
| **ADAPT + pares de espín, 15 operadores (la curva)** | **15** | **194** | **~530 (transpiler) / 378 (fSWAP)** | **0,38 mHa** | **6,5·10⁻³** |

**Por qué funcionó.** En H₄ la correlación la llevan casi toda las dobles $\alpha\beta$ entre los dos orbitales ocupados y los dos vacíos. ADAPT las encuentra primero porque son las de mayor gradiente, y descarta las que UCCSD ×2 incluye "por las dudas". Bajar de 52 a 15 parámetros también hace la optimización más rápida y estable.

### 3.4 Qué no funcionó

- **Umbral de gradiente de 10⁻⁴** como criterio de parada: ADAPT se frena en mesetas a 1,6, 2,8 y 3,0 Å con 0,6–2 mHa de error, y después el gradiente vuelve a crecer. Con 10⁻⁵ llega a FCI en toda la curva. Por eso cortamos con la varianza.
- **Agregar las simples de espín al final de la lista** (y no al lado de su pareja), con el *warm start* entre distancias: cae en un mínimo malo (2,5 mHa, $\langle S^2\rangle$ = 0,52).
- **Reordenar los operadores** para necesitar menos fSWAP (314 CZ en línea): el error sube a 5,2 mHa. El orden de los operadores cambia el ansatz.
- **Co-ADAPT** (elegir cada operador por gradiente dividido por su costo en línea): necesita 16–19 operadores para precisión química y termina en ~460–470 CZ en línea. En H₄ casi todas las dobles hacen falta, así que penalizar por costo solo agrega operadores.
- **Circuitos explícitos de excitación sin ruteo propio** (escalera de CNOT + RY multicontrolada): exactos, pero el transpiler los rutea peor que `PauliEvolutionGate` (705 contra 526 CZ en línea).
- **Probar layouts y semillas del transpiler**: el mejor resultado fue ~515 CZ, contra ~526 por defecto.
- **Pool QEB** (excitaciones de qubit, sin cadenas Z): 159 CZ con conectividad total, pero ~451 en línea, peor que el ruteo fSWAP, y rompe más el espín ($\langle S^2\rangle$ hasta 0,074). Queda en la comparación de ansätze (§11).

La celda de abajo verifica el ansatz de la curva a 2,0 Å: $\hat H$ contra PySCF, que el circuito de Qiskit y el simulador numpy den el mismo estado, y las CZ en las tres variantes.
