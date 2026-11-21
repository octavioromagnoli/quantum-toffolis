## 1 · El problema: de la molécula a los qubits
### 1.1 Qué hay que calcular

Cuatro átomos de H en línea recta, separados por la misma distancia $d$, en la base STO-3G. Con la aproximación de Born-Oppenheimer los núcleos quedan fijos y, para cada $d$, buscamos el autovalor más chico del hamiltoniano electrónico:

$$\hat H(d)\,\lvert\Psi_0(d)\rangle = E_0(d)\,\lvert\Psi_0(d)\rangle, \qquad E_0(d) = \min_{\Psi}\ \frac{\langle\Psi\rvert\hat H(d)\lvert\Psi\rangle}{\langle\Psi\vert\Psi\rangle}.$$

De la curva $E_0(d)$ salen la distancia de equilibrio y la energía de disociación:

$$R_e = \arg\min_d E_0(d), \qquad D_e = E_0(3{,}0\ \text{Å}) - E_0(R_e).$$

En la zona estirada ($d \ge 1{,}8$ Å) la consigna pide

$$\underbrace{\lvert E_{\rm VQE} - E_{\rm FCI}\rvert + \sigma < \lvert E_{\rm CCSD} - E_{\rm FCI}\rvert}_{\text{ganarle a CCSD}} \qquad\text{y, como meta máxima,}\qquad \underbrace{\lvert E_{\rm VQE} - E_{\rm FCI}\rvert < 1{,}6\ \text{mHa}}_{\text{precisión química}}.$$

FCI es la diagonalización exacta en la base: con 8 qubits se puede calcular, pero **solo sirve para comparar**, no como resultado.

### 1.2 El hamiltoniano en segunda cuantización

STO-3G pone un orbital 1s en cada átomo, así que hay 4 orbitales espaciales. Hartree-Fock restringido (RHF, con PySCF) los combina en 4 orbitales moleculares $\phi_0,\dots,\phi_3$ ordenados por energía, y cada uno admite un electrón $\alpha$ y uno $\beta$: son **8 espín-orbitales para 4 electrones**. En esa base el hamiltoniano electrónico es

$$\boxed{\ \hat H = E_{\rm nuc} + \sum_{\sigma\in\{\alpha,\beta\}}\ \sum_{p,q=0}^{3} h_{pq}\, a^\dagger_{p\sigma} a_{q\sigma} + \frac12 \sum_{\sigma,\tau\in\{\alpha,\beta\}}\ \sum_{p,q,r,s=0}^{3} (pq\vert rs)\, a^\dagger_{p\sigma}\, a^\dagger_{r\tau}\, a_{s\tau}\, a_{q\sigma}\ }$$

con las integrales que calcula PySCF (los orbitales son reales):

$$h_{pq} = \int \phi_p(\mathbf r)\left(-\tfrac12\nabla^2 - \sum_{A=1}^{4} \frac{1}{\lvert\mathbf r - \mathbf R_A\rvert}\right)\phi_q(\mathbf r)\,d\mathbf r, \qquad (pq\vert rs) = \iint \frac{\phi_p(\mathbf r_1)\phi_q(\mathbf r_1)\,\phi_r(\mathbf r_2)\phi_s(\mathbf r_2)}{\lvert\mathbf r_1 - \mathbf r_2\rvert}\,d\mathbf r_1\,d\mathbf r_2,$$

y $E_{\rm nuc} = \sum_{A<B} 1/\lvert\mathbf R_A - \mathbf R_B\rvert$, la repulsión entre núcleos (en unidades atómicas; las energías van en Hartree, 1 Ha ≈ 27,2 eV).

- $h_{pq}$ es la parte de **un electrón**: energía cinética más atracción de los 4 núcleos.
- $(pq\vert rs)$ es la repulsión entre **dos electrones**. Los casos $(pp\vert qq)$ son la repulsión de Coulomb $J_{pq}$, y $(pq\vert qp)$ el intercambio $K_{pq}$.
- $a^\dagger_{p\sigma}$ y $a_{p\sigma}$ crean y destruyen un electrón en el espín-orbital $p\sigma$. Cumplen $\{a_i, a_j^\dagger\} = \delta_{ij}$ y $\{a_i, a_j\} = 0$: el anticonmutador es el principio de exclusión de Pauli.

Es la misma expresión de la consigna, $\sum h_{pq}a_p^\dagger a_q + \tfrac12\sum h_{pqrs}a_p^\dagger a_q^\dagger a_r a_s$, escrita con la notación de químicos de PySCF ($h_{prsq} = (pq\vert rs)$, con los índices de espín explícitos). `hamiltonian.qubit_hamiltonian` la implementa tal cual.

### 1.3 De fermiones a qubits: Jordan-Wigner

Cada espín-orbital es un qubit: $\lvert 1\rangle$ = ocupado y $\lvert 0\rangle$ = vacío. Convención: **qubits 0–3 = orbitales $\alpha$ y 4–7 = orbitales $\beta$** (en cada bloque, el orbital 0 es el de menor energía). El determinante de Hartree-Fock ocupa los dos orbitales más bajos de cada espín: qubits 0, 1, 4 y 5, que en el orden de Qiskit (el qubit 0 a la derecha) es `|00110011⟩`.

La transformación de Jordan-Wigner escribe los operadores fermiónicos con matrices de Pauli:

$$a_j = Z_0 Z_1 \cdots Z_{j-1}\ \frac{X_j + iY_j}{2}, \qquad a_j^\dagger = Z_0 Z_1 \cdots Z_{j-1}\ \frac{X_j - iY_j}{2}.$$

La cadena de $Z$ pone el signo $(-1)^{n_0 + \dots + n_{j-1}}$, que hace que los $a_j$ anticonmuten. De ahí salen los tres bloques con los que se arma $\hat H$:

$$a_j^\dagger a_j = \hat n_j = \frac{I - Z_j}{2}, \qquad a_p^\dagger a_q + a_q^\dagger a_p = \frac12\left(X_p\, Z_{p+1}\cdots Z_{q-1}\, X_q + Y_p\, Z_{p+1}\cdots Z_{q-1}\, Y_q\right) \quad (p<q),$$

y los productos de cuatro operadores, que dejan cadenas con cuatro $X$ o $Y$. El resultado es una suma de cadenas de Pauli de 8 qubits con coeficientes reales que dependen de $d$:

$$\hat H(d) = \sum_{k} c_k(d)\, P_k, \qquad P_k \in \{I, X, Y, Z\}^{\otimes 8}.$$

La celda de abajo arma $\hat H$ a 0,9 Å (cerca del equilibrio) y a 2,4 Å (estirada), lo escribe con sus términos más grandes y los clasifica.
