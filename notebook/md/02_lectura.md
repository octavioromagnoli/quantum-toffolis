### 1.4 Cómo leer el hamiltoniano de H₄

Los 185 términos se reparten en cinco clases (tabla de arriba), y cada una tiene un significado físico:

$$\hat H = \underbrace{c_0\, I}_{1} \;+\; \underbrace{\sum_{p} c_p\, Z_p}_{8} \;+\; \underbrace{\sum_{p<q} c_{pq}\, Z_p Z_q}_{28} \;+\; \underbrace{\sum c\,\big(X Z\cdots Z X,\ Y Z\cdots Z Y,\ \dots\big)}_{56:\ \text{un electrón salta}} \;+\; \underbrace{\sum c\,\big(XXYY,\ XYYX,\ \dots\big)}_{92:\ \text{dos electrones saltan}}$$

- **$I$, $Z_p$ y $Z_pZ_q$ (37 términos, diagonales):** salen de $\hat n_p = (I - Z_p)/2$ y de $\hat n_p\hat n_q$. Dan la energía de cada determinante (cada bitstring): energía de los orbitales ocupados más la repulsión de Coulomb $J_{pq}$ y el intercambio $K_{pq}$ entre ellos. El determinante de menor energía diagonal es el de Hartree-Fock, `|00110011⟩`.
- **2 X/Y (56 términos):** $a^\dagger_p a_q + a^\dagger_q a_p$, un electrón que salta del orbital $q$ al $p$. Las $Z$ que lo acompañan son la cadena de Jordan-Wigner y, a veces, una ocupación de la que depende el salto ($\hat n_r\, a^\dagger_p a_q$).
- **4 X/Y (92 términos):** $a^\dagger_p a^\dagger_r a_s a_q$ con los cuatro índices distintos: **dos electrones que saltan a la vez**. Son las excitaciones dobles, las que mezclan HF con los determinantes que le siguen.
- Todos los términos conservan $N_\alpha$ y $N_\beta$, así que las $X$/$Y$ aparecen de a pares dentro de cada bloque de espín. Eso hace que $Z_0Z_1Z_2Z_3$ y $Z_4Z_5Z_6Z_7$ conmuten con $\hat H$; lo usamos en §7 para sacar qubits.

**Qué cambia al estirar.** Cerca del equilibrio (0,9 Å) dominan los términos diagonales: la suma de $\lvert c\rvert$ de los $Z_p$ es 1,64 Ha y la brecha HOMO-LUMO, 0,76 Ha. HF es una muy buena primera aproximación (pesa 0,95 en el estado exacto). A 2,4 Å los $Z_p$ caen a 0,31 Ha y la brecha a 0,24 Ha: los orbitales se acercan en energía porque los átomos casi no se ven. En cambio, los saltos dobles (4 X/Y) **no se achican** (2,13 → 2,22 Ha). Muchos determinantes con energías parecidas quedan conectados con la misma fuerza, y el estado exacto pasa a ser una mezcla de todos ellos. Eso es la **correlación estática**, y es lo que rompe a CCSD.

### 1.5 Validación del hamiltoniano

En las 18 distancias de la curva comparamos $\hat H$, ya escrito en qubits, contra PySCF (primera línea de la salida de abajo):

- $\langle \mathrm{HF}\rvert\hat H\lvert \mathrm{HF}\rangle$ coincide con $E_{\rm HF}$ con una diferencia menor que $10^{-14}$ Ha.
- El autovalor más bajo de $\hat H$ **restringido al sector $N_\alpha = N_\beta = 2$** (36 de los 256 estados de la base) coincide con $E_{\rm FCI}$ (CASCI completo) con la misma precisión, y su estado es un singulete ($\langle S^2\rangle = 0$). Restringimos al sector porque $\hat H$ en qubits también describe estados con otra cantidad de electrones (H₄⁺, H₄⁻, …), y el que nos interesa es el de 4 electrones con $S_z = 0$.

Con esto, cualquier diferencia entre VQE y FCI viene del ansatz, de la optimización o de las mediciones, no del hamiltoniano.
