**Qué dio** (gráfico y tabla de arriba):

- **Cuadrado.** HF ya se equivoca en 124 mHa a 0,7 Å (en el zigzag, en 37 mHa). A 0,7 Å el determinante HF tiene **peso nulo** en el estado exacto (columna `w_hf`), por simetría. CCSD queda 9,4 mHa **por encima** de FCI a 0,7 Å, y 6,5, 18,9 y 14,3 mHa por debajo a 1,2, 1,8 y 2,4 Å. A 3,0 Å no converge, así que su número ahí no sirve (cruz negra en el gráfico).
- **Zigzag.** Se comporta como la cadena recta: CCSD es exacto cerca del equilibrio ($w_{\rm HF}$ = 0,98) y falla al estirar, hasta 31,6 mHa por debajo de FCI a 2,4 Å.
- **VQE.** Sin ruido llega a FCI en todos los puntos de las dos geometrías (columna `err_sin_ruido_mHa`), **incluso donde HF pesa cero**: las excitaciones dobles de UCCSD conectan HF con los determinantes que sí pesan. Con la energía medida (2 M de shots por punto), los errores son del orden de la $\sigma$. En la zona estirada, VQE le gana a CCSD con el criterio de la consigna en todos los puntos donde CCSD converge (columna `beats_ccsd`). En el cuadrado a 3,0 Å, CCSD no converge y su error ahí puede ser chico por casualidad, así que ese punto no cuenta.
- **Costo:** es el mismo circuito que UCCSD ×2 en la cadena (582 CZ, 1648 en línea), porque solo depende de la cantidad de orbitales y de electrones.

**Por qué funcionó.** VQE no necesita que HF sea una buena aproximación: solo necesita que el ansatz pueda llegar al estado exacto desde HF. CCSD, en cambio, se apoya en que HF domine, y falla justo en los casos de correlación estática, sea por estirar o por la simetría del cuadrado.

**Limitación:** con 2 M de shots, $\sigma \approx$ 0,6–1,1 mHa no alcanza para afirmar precisión química punto por punto; para eso harían falta los 8 M de shots de la curva principal.
