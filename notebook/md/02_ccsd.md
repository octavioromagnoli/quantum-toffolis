### 1.6 Por qué CCSD falla al estirar

El gráfico de arriba muestra el error de CCSD contra FCI (izquierda) y el peso del determinante HF en el estado exacto, $w_{\rm HF} = \lvert\langle \mathrm{HF}\vert\Psi_{\rm FCI}\rangle\rvert^2$ (derecha). La zona gris es la zona estirada ($d \ge 1{,}8$ Å).

CCSD escribe el estado como $\lvert\Psi_{\rm CC}\rangle = e^{\hat T}\lvert \mathrm{HF}\rangle$, con $\hat T = \hat T_1 + \hat T_2$ (excitaciones simples y dobles), y saca la energía y las amplitudes de ecuaciones **proyectadas**:

$$E_{\rm CCSD} = \langle \mathrm{HF}\rvert\, e^{-\hat T}\hat H e^{\hat T}\,\lvert \mathrm{HF}\rangle, \qquad \langle \Phi_{i}^{a}\rvert\, e^{-\hat T}\hat H e^{\hat T}\,\lvert \mathrm{HF}\rangle = \langle \Phi_{ij}^{ab}\rvert\, e^{-\hat T}\hat H e^{\hat T}\,\lvert \mathrm{HF}\rangle = 0.$$

- **Cerca del equilibrio** (0,9 Å) el estado exacto es casi HF ($w_{\rm HF}$ = 0,95) y CCSD es prácticamente exacto (+0,007 mHa).
- **Al estirar**, cada electrón se queda en su átomo y el estado exacto pasa a ser una mezcla de muchos determinantes con pesos parecidos: $w_{\rm HF}$ baja a 0,58 a 1,8 Å y a 0,27 a 3,0 Å (§1.4). CCSD solo incluye las excitaciones simples y dobles *conectadas*. Las triples y cuádruples aparecen solo como productos de las dobles ($\tfrac12\hat T_2^2$, …), y para 4 electrones FCI equivale a CCSDTQ: eso deja de alcanzar cuando HF ya no domina.
- **CCSD no es variacional.** Como $e^{-\hat T} \neq (e^{\hat T})^\dagger$, $E_{\rm CCSD}$ no es el valor esperado $\langle\Psi\rvert\hat H\lvert\Psi\rangle/\langle\Psi\vert\Psi\rangle$ de ningún estado, y no es una cota superior. Acá cae **por debajo** de FCI a partir de 1,1 Å: −9,1 mHa a 1,8 Å y −30,6 mHa a 2,6 Å, unas 19 veces la precisión química. Quedar por debajo no es "mejor": el error en valor absoluto es grande, y el signo muestra que el método salió de su rango de validez. CCSD convergió en las 18 distancias, así que no es un problema numérico.

**Por qué VQE le puede ganar.** VQE mide el valor esperado de un estado preparado por un circuito, así que para cualquier $\boldsymbol\theta$

$$E(\boldsymbol\theta) = \langle\psi(\boldsymbol\theta)\rvert\,\hat H\,\lvert\psi(\boldsymbol\theta)\rangle = \sum_n \lvert\langle E_n\vert\psi(\boldsymbol\theta)\rangle\rvert^2\, E_n \;\ge\; E_0 .$$

No puede quedar por debajo de FCI, salvo por el ruido estadístico de los shots, que cubre la barra $\pm\sigma$. Entre dos estados medidos, el de energía más baja es el mejor. Esa garantía es la base del argumento de confiabilidad (§6).
