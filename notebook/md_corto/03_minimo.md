## 2 · Lo mínimo que pide la consigna: VQE con UCCSD

Nuestra primera versión completa fue el VQE "de libro". Sus tres piezas siguen en la versión final; lo que cambia es cómo resolvemos cada una.

### 2.1 Dos etapas

1. **Buscar los parámetros $\boldsymbol\theta^*$ sin ruido** (la consigna lo permite), con un simulador exacto del mismo circuito y gradiente analítico.
2. **Medir $E(\boldsymbol\theta^*)$ con shots** en `SamplerV2`. Esa es la energía que se reporta, con su $\sigma$.

### 2.2 El ansatz UCCSD

Parte de Hartree-Fock y aplica todas las excitaciones simples y dobles: **8 simples + 18 dobles = 26 parámetros** por paso. Cada excitación es un factor exacto del circuito, $e^{-i\theta_k \hat G_k}$; "UCCSD ×2" repite la lista dos veces. En la curva, cada distancia arranca de los $\boldsymbol\theta^*$ de la anterior (*warm start*).

### 2.3 Medir la energía con shots

Un shot da un bitstring, no una energía. Por eso las cadenas de Pauli de $\hat H$ se juntan en **grupos** que se miden con un mismo circuito: en cada qubit, todas tienen la misma letra o $I$ (grupos QWC). Se miden rotando cada qubit ($H$ para $X$, $S^\dagger H$ para $Y$), **sin CZ extra**.

$$\hat E = c_0 + \sum_{g} \bar o_g, \qquad \sigma^2 = \sum_{g} \frac{s_g^2}{N_g}$$

$\bar o_g$ y $s_g^2$ son la media y la varianza de las muestras del grupo $g$, y $N_g$ sus shots. **La σ sale de las propias muestras.** En la versión mínima: 46 grupos y 8 M de shots por punto.

### 2.4 Lo mínimo, recalculado

La celda de abajo compara UCCSD ×1 y ×2 en cuatro distancias estiradas: error contra FCI, $\langle S^2\rangle$ (0 = singulete), CZ y probabilidad de que ninguna CZ falle con 0,1 % de error por CZ. Usa el estado exacto, solo para elegir el ansatz.
