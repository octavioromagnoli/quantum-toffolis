## 2 · Lo mínimo que pide la consigna: VQE con UCCSD  *(11:09, Joaquín, Gonzalo y Santino)*

La primera versión que corrió de punta a punta (commit de las 11:09) fue el VQE "de libro". Sus tres piezas siguen siendo las mismas en la versión final; lo que cambia es cómo resolvemos cada una.

### 2.1 El ciclo variacional

$$\lvert \mathrm{HF}\rangle \xrightarrow{\ U(\boldsymbol\theta)\ } \lvert\psi(\boldsymbol\theta)\rangle \xrightarrow{\ \text{shots}\ } E(\boldsymbol\theta) \pm \sigma \xrightarrow{\ \text{optimizador clásico}\ } \boldsymbol\theta' \;\longrightarrow\; \cdots$$

La consigna deja buscar los parámetros como queramos, incluso sin ruido. Por eso separamos dos etapas:

1. **Buscar $\boldsymbol\theta^*$ sin ruido**, con un simulador numpy exacto del mismo circuito y gradiente analítico (L-BFGS-B).
2. **Medir $E(\boldsymbol\theta^*)$ con shots** en `SamplerV2` de `qiskit_aer`. Esa es la energía que se reporta, con su $\sigma$.

### 2.2 El ansatz UCCSD

UCCSD aplica a HF la exponencial de las excitaciones simples y dobles:

$$\lvert\psi(\boldsymbol\theta)\rangle = e^{\hat T(\boldsymbol\theta) - \hat T^\dagger(\boldsymbol\theta)}\lvert \mathrm{HF}\rangle, \qquad \hat T(\boldsymbol\theta) = \sum_{i\in\text{occ},\,a\in\text{vir}} \theta_i^a\, a_a^\dagger a_i \;+\; \sum_{\substack{i<j\,\in\,\text{occ}\\ a<b\,\in\,\text{vir}}} \theta_{ij}^{ab}\, a_a^\dagger a_b^\dagger a_j a_i .$$

Con 2 electrones $\alpha$ y 2 $\beta$ en 4 orbitales, y conservando $S_z$, hay **8 simples** (4 por espín) y **18 dobles** (16 $\alpha\beta$, 1 $\alpha\alpha$ y 1 $\beta\beta$): 26 parámetros. Como las excitaciones no conmutan entre sí, el circuito usa la **factorización de Trotter**:

$$U(\boldsymbol\theta) = \prod_{r=1}^{n_{\rm rep}}\ \prod_{k=1}^{26} e^{-i\theta_{k,r}\hat G_k}, \qquad \hat G_k = i\big(\hat\tau_k - \hat\tau_k^\dagger\big),$$

donde $\hat\tau_k$ es la excitación $k$ y $n_{\rm rep}$ la cantidad de pasos ("UCCSD ×1", "×2"). Cada $\hat G_k$, escrito con Jordan-Wigner, es una suma de cadenas de Pauli **que conmutan entre sí**, así que `PauliEvolutionGate` implementa cada factor **exactamente**, sin error de Trotter adentro del factor. Además $\hat G_k^3 = \hat G_k$, lo que da una forma cerrada para el simulador numpy:

$$e^{-i\theta\hat G} = I + (\cos\theta - 1)\,\hat G^2 - i\sin\theta\ \hat G .$$

El gradiente sale con un barrido hacia atrás ("adjunto"), al costo de unas dos evaluaciones de la energía:

$$\frac{\partial E}{\partial\theta_k} = 2\,\mathrm{Re}\,\big\langle \lambda_k \big\rvert\, (-i\hat G_k)\, \big\lvert\psi_k\big\rangle, \qquad \lvert\psi_k\rangle = U_k\cdots U_1\lvert \mathrm{HF}\rangle,\quad \langle\lambda_k\rvert = \langle\psi\rvert\hat H\, U_n\cdots U_{k+1}.$$

En la curva, cada distancia arranca del $\boldsymbol\theta^*$ de la anterior (*warm start*).

### 2.3 Medir la energía con shots

Un shot no devuelve una energía: devuelve un bitstring. Por eso juntamos las cadenas de Pauli de $\hat H$ en **grupos** que se pueden medir con un mismo circuito:

$$\hat H = c_0\, I + \sum_{g=1}^{G} \hat O_g, \qquad \hat O_g = \sum_{k\in g} c_k P_k .$$

En un grupo QWC (*qubit-wise commuting*), en cada qubit todas las cadenas tienen la misma letra o $I$. Se mide rotando cada qubit a esa base ($H$ para $X$, $S^\dagger H$ para $Y$), **sin CZ extra**. Cada shot $s$ del grupo $g$ da el bitstring $x_s$ y el valor $o_g(x_s) = \sum_{k\in g} c_k\, (-1)^{\,x_s\cdot\mathrm{sop}(P_k)}$. Con $N_g$ shots por grupo:

$$\hat E = c_0 + \sum_{g} \bar o_g, \qquad \sigma^2 = \sum_{g} \frac{s_g^2}{N_g},$$

donde $\bar o_g$ es la media muestral y $s_g^2$ la varianza muestral de $o_g$ (incluye las covarianzas entre las cadenas del grupo). Los grupos usan shots distintos, así que sus errores son independientes y las varianzas se suman. **La σ sale de las propias muestras**, no de una fórmula fija.

En la versión mínima: 46 grupos QWC (el coloreo de `group_commuting` de Qiskit) y **8 M de shots por punto**, repartidos en dos etapas (piloto + Neyman, ver §4).

### 2.4 Lo mínimo, recalculado

La celda de abajo arma UCCSD ×1 y ×2 en cuatro distancias estiradas y calcula, **sin shots y con el estado exacto** (solo para elegir el ansatz), el error contra FCI, $\langle S^2\rangle$ (0 para el singulete) y el costo: CZ con conectividad total y con qubits en línea, y la probabilidad de que ninguna CZ falle con 0,1 % de error por CZ.
