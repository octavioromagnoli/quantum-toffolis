## 8 · Mitigación de errores: precisión química con ruido

**Qué vimos.** Con 16 CZ, el compacto queda a ~+14 mHa de FCI. En hardware se suma el **error de lectura** (bits medidos al revés). Agregamos al modelo **2 % de error de lectura por qubit**.

**Cuatro técnicas** (`mitigation.py`):

| técnica | qué corrige | cómo |
|---|---|---|
| **RO**: corrección de lectura | bits leídos al revés | se calibra cada qubit midiendo todo en 0 y todo en 1, y se corrigen los valores de cada grupo |
| **PS**: post-selección | errores que cambian $N_\alpha$ o $N_\beta$ | se descartan los shots con $N \neq 2$ (solo en los grupos que miden en $Z$) |
| **ZNE**: extrapolación a ruido cero | ruido de compuertas | se repite el circuito con 3× y 5× más ruido ($C \to C\,C^\dagger C$) y se extrapola a ruido 0 |
| **eco + ZNE** | ídem, más fino | $C\,C^\dagger$ mide la fidelidad y corrige cada nivel antes de extrapolar |

**Hallazgo:** hay que plegar también la **rotación de medición** de cada grupo, no solo el ansatz. Si no, queda un sesgo que la extrapolación no saca (Richardson a 1,8 Å: de +0,16 a +0,002 mHa).

La celda de abajo calcula el sesgo **exacto** de cada combinación a 2,4 Å (sin shots) y cuántos shots necesita cada una para σ = 0,5 mHa: sirve para elegir antes de gastar shots.
