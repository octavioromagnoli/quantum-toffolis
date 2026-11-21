## 12 · ¿Y la nitrogenasa? Cuánto costaría la molécula real
La consigna arranca con la nitrogenasa: su cofactor **FeMoco** es el sistema que motiva todo esto. Como cierre, medimos **cómo crecen los recursos de VQE** al alargar la cadena de hidrógenos (H₄, H₆, …, H₁₄) y los extrapolamos al tamaño de FeMoco. Es la parte del *análisis de limitaciones*: hasta dónde llega el método y por qué.

**Hipótesis (a confirmar o refutar con los datos).** VQE escala de forma polinomial en qubits y compuertas, pero la cantidad de mediciones y la acumulación de ruido lo vuelven impracticable para FeMoco.

**La perilla: cadenas H$_N$** lineales y equiespaciadas en STO-3G, igual que el desafío: $N$ átomos → $N$ orbitales → $2N$ qubits y $N$ electrones. Hay un electrón por orbital, igual que en el modelo (54e, 54o) de FeMoco de Reiher *et al.* (2017), así que las fórmulas valen para $N$ = 54. Usamos dos distancias: $d$ = 0,894 Å (el $R_e$ de H₄) y 2,0 Å (estirada).

**Qué se cuenta para cada $N$** (`h4vqe/scaling.py`):

- **Términos de Pauli, grupos QWC (*sorted insertion*) y $\sum\lvert c_k\rvert$**: se construye $\hat H$ y se cuenta, hasta $N$ = 14 (28 qubits).
- **Parámetros de UCCSD** con la fórmula exacta. Con $o = v = N/2$ orbitales ocupados y vacíos por espín:

$$n_\theta = \underbrace{2\,o\,v}_{\text{simples}} + \underbrace{(o\,v)^2}_{\text{dobles }\alpha\beta} + \underbrace{2\binom{o}{2}\binom{v}{2}}_{\text{dobles }\alpha\alpha,\ \beta\beta}\qquad (N = 4:\ 8 + 16 + 2 = 26).$$

- **CZ de UCCSD ×1**, transpiladas como en §2 (hasta $N$ = 14), y una cota exacta por fórmula: 4 CZ por simple y 14 por doble, más 2 por cada qubit espectador de la cadena de Jordan-Wigner.
- **Shots por evaluación de energía** para $\sigma$ = 1,6 mHa, con Neyman (§4.1). La cota de peor caso, $\sqrt{V_g} \le \sum_{k\in g}\lvert c_k\rvert$, da

$$M \;=\; \frac{\big(\sum_g \sqrt{V_g}\big)^2}{\varepsilon^2} \;\le\; \frac{\big(\sum_{k}\lvert c_k\rvert\big)^2}{\varepsilon^2}.$$

  Calculamos el valor **exacto** con el estado FCI y los grupos reales hasta $N$ = 12 (24 qubits), y la cota para todo $N$.
- **Probabilidad de que ninguna CZ falle** con 0,1 % por CZ: $\log_{10} P = n_{\rm CZ}\,\log_{10}(1 - 10^{-3})$.
- **CCSD contra FCI** hasta $N$ = 14 (FCI de 11,8 millones de determinantes).

El cálculo pesado (N = 10 a 14, ~1 h en una PC) está precalculado (`scripts/scaling.py` → `results/scaling/`), y sus datos van en `h4vqe/escalado_recursos.csv`, que se escribe al principio del notebook. La celda de abajo **recalcula en vivo** $N$ = 4, 6 y 8 y los compara con lo guardado, para mostrar que es reproducible. Después ajusta leyes de potencia $y = a\,N^b$ con $N \ge 6$ y extrapola a FeMoco.
