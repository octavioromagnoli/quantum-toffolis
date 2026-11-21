## 6 · Problema 4: ¿por qué confiar en la zona estirada sin mirar FCI?  *(11:42 → 12:58, Joaquín)*

**Observación.** La consigna pide argumentar la confiabilidad **antes** de comparar con FCI. En lo mínimo, la varianza, las simetrías y la comparación entre profundidades se calculaban con el estado del simulador, que en una computadora cuántica no existe.

**Hipótesis.** Todo eso se puede **medir con shots**. Para la varianza, que es una resta de números grandes, un eco de Loschmidt cuenta eventos raros y resuelve valores mucho más chicos.

En tres puntos estirados (1,8, 2,4 y 3,0 Å) comparamos cuatro estados preparados con circuitos y **medidos con shots**: **HF** ($\boldsymbol\theta = 0$), **ADAPT con 12 operadores** (el paso anterior al corte, 169 CZ), **el ansatz de la curva** (ADAPT con 15 operadores, 194 CZ) y **UCCSD ×2** (582 CZ), mucho más profundo. FCI aparece recién en la última tabla, como comprobación.

**1. Principio variacional.** Para cualquier estado, $E(\boldsymbol\theta) \ge E_0$ salvo ruido estadístico (§1.6). Entre estados variacionales, el que da la energía medida más baja es el mejor. CCSD no es variacional: en la zona estirada queda *por debajo* de la exacta, y ahí "más bajo" no quiere decir "mejor".

**2. Límite físico (sin FCI).** A 3,0 Å los átomos casi no interactúan, así que la energía no puede estar muy por debajo de la de 4 átomos de H separados, $4\,E_{\rm H} = -1{,}866327$ Ha. Con un solo electrón, HF es exacto en la base, así que es un número exacto y barato. La celda calcula cuánto queda cada método por debajo de ese límite.

**3. Varianza de la energía medida por grupos.** En un autoestado, $\mathrm{Var}(\hat H) = \langle\hat H^2\rangle - \langle\hat H\rangle^2 = 0$. La medimos como

$$\mathrm{Var}(\hat H) = \big\langle(\hat H - c)^2\big\rangle - \big(\langle\hat H\rangle - c\big)^2,$$

con $\langle(\hat H - c)^2\rangle$ medido con el mismo estimador que la energía: $(\hat H - c)^2$ tiene 1775 cadenas de Pauli, que se miden en 225 grupos QWC. Tomamos $c$ = término constante de $\hat H$ porque reduce la $\sigma$ a la mitad respecto de $c = E$. Con 2 M de shots, $\sigma \approx$ 1,5–2·10⁻³ Ha²: separa HF (0,12–0,18 Ha²) de los ansätze (compatibles con 0), pero no a los ansätze entre sí. Como es una resta, por ruido puede dar negativa; en la tabla se muestra recortada en 0, con su $\sigma$.

**4. Varianza por eco de Loschmidt.** El circuito prepara $\lvert\psi\rangle$, lo evoluciona con $e^{-i\hat Ht}$ y lo des-prepara con $U(\boldsymbol\theta)^\dagger$. La probabilidad de **no** volver a $\lvert 0\cdots 0\rangle$ es

$$P_{\rm eco} = 1 - \big\lvert\langle\psi\rvert e^{-i\hat Ht}\lvert\psi\rangle\big\rvert^2 = \mathrm{Var}(\hat H)\,t^2 + O(t^4).$$

En vez de restar dos números grandes, cuenta **eventos raros**: si $\psi$ es un autoestado, ningún shot vuelve "mal". Con $N$ shots y $k$ eventos, la cota superior del 95 % sale de Clopper-Pearson (con $k = 0$, $\mathrm{Var} \lesssim 3/(Nt^2)$). Usamos $t$ = 1 Ha⁻¹ y Suzuki-Trotter de 2º orden con 4 pasos. El error de Trotter lo verificamos sin FCI, comparando 2 contra 4 pasos: escala como $1/n^4$ y con 4 pasos es 9·10⁻⁸, debajo de la resolución de 10⁶ shots. A orden dominante el eco *subestima* la varianza (Mandelstam-Tamm); a $t$ = 1 la diferencia es ~8 % para UCCSD ×1.
*Costo:* el circuito de eco tiene $U(\boldsymbol\theta)$, la evolución y $U(\boldsymbol\theta)^\dagger$, más de 4000 CZ (columna `echo_cz`). **Es un diagnóstico de simulador sin ruido**: con 0,1 % de error por CZ, el ruido taparía la señal. En hardware, lo equivalente es la varianza por grupos con más shots.
*Qué implica:* por la cota de Weinstein hay un autovalor exacto de $\hat H$ a menos de $\sqrt{\mathrm{Var}}$ de la energía medida. Esa cota es floja; la que importa es la de segundo orden: si hay una brecha $\Delta$ al primer excitado del mismo sector,

$$E - E_0 \;\lesssim\; \frac{\mathrm{Var}(\hat H)}{\Delta}.$$

Con $\mathrm{Var} \sim 10^{-4}$ Ha², ese error queda debajo de la $\sigma$ de los shots (0,5 mHa) si $\Delta \gtrsim$ 0,2 Ha. No medimos $\Delta$; el punto 6 lo chequea de otra forma.

**5. Simetrías medidas.** H₄ tiene un fundamental singulete con $N_\alpha = N_\beta = 2$. Medimos $\langle N_\alpha\rangle$, $\langle N_\beta\rangle$ y $\langle S^2\rangle$ con shots, donde

$$\hat S^2 = \hat S_-\hat S_+ + \hat S_z^2 + \hat S_z, \qquad \hat S_+ = \sum_p a^\dagger_{p\alpha}a_{p\beta},\quad \hat S_z = \tfrac12(\hat N_\alpha - \hat N_\beta).$$

$S^2$ es la prueba que discrimina estados de espín mezclado. Al estirar, el singulete y el triplete se acercan (a 3,0 Å están a 0,95 mHa): una energía "buena" puede ser el estado equivocado. UCCSD ×1 converge a una mezcla, con $\langle S^2\rangle$ medido de 0,64 a 2,4 Å y 0,85 a 3,0 Å, y a 3,0 Å su energía medida no se distingue de la buena. La energía sola no lo delata; la $S^2$ medida sí.

**6. Convergencia con la profundidad.** Si el ansatz de la curva convergió, uno menos profundo tiene que dar una energía medida más alta y una varianza mayor, y uno mucho más profundo no tiene que bajar la energía más allá de la $\sigma$.

**7. Barras de error calibradas** (§4.3): $z$ tiene dispersión 1,01 en 3600 repeticiones.
