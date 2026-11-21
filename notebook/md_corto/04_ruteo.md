### 3.2 Qué dio ADAPT

- **Corte en 13 operadores:** al pasar de 12 a 13, la varianza en la zona estirada cae 12 veces y el error máximo baja de 3,7 a 0,6 mHa. Son **185 CZ**.
- **Pares de espín:** ADAPT eligió dos excitaciones simples sin su pareja de espín opuesto, y el estado rompía un poco el espín. Agregando las dos parejas (**15 operadores, 194 CZ**), el error máximo baja a **0,38 mHa** y $\langle S^2\rangle$ a la mitad.
- **Dos arranques:** el optimizador prueba el *warm start* y también $\boldsymbol\theta = 0$, y se queda con la energía más baja (evita mínimos malos).

### 3.3 Qubits en línea: ruteo con fSWAP

Con qubits en línea, el transpiler de Qiskit agrega SWAP y ADAPT pasa de 194 a ~530 CZ. Nuestro ruteo (`routing.py`) mueve los modos con **fSWAP** (2 CZ, lleva consigo el signo fermiónico), los junta antes de cada excitación y no deshace la permutación. Mismo estado, **378 CZ en línea** (−29 %).

| ansatz | parámetros | CZ | CZ en línea | error máx. sin ruido |
|---|---|---|---|---|
| UCCSD ×1 | 26 | 292 | 806 | 4,5 mHa |
| UCCSD ×2 (lo mínimo) | 52 | 582 | 1648 | < 0,001 mHa |
| ADAPT, 13 ops | 13 | 185 | 370 | 0,59 mHa |
| **ADAPT, 15 ops (la curva)** | **15** | **194** | **378** | **0,38 mHa** |

**Por qué funcionó:** en H₄ casi toda la correlación la llevan unas pocas dobles; ADAPT las encuentra primero y descarta el resto.

**Qué no funcionó:** cortar ADAPT por gradiente < 10⁻⁴ (se frena en mesetas); reordenar operadores para ahorrar fSWAP (el error sube a 5,2 mHa); elegir operadores por gradiente/costo (termina con más CZ); y el pool de excitaciones de qubit (más CZ en línea y más ruptura de espín, ver §11).

La celda de abajo verifica el ansatz de la curva a 2,0 Å: $\hat H$ contra PySCF, circuito contra simulador y CZ.
