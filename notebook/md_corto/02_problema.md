## 1 · El problema: de la molécula a los qubits
### 1.1 Qué hay que calcular

Cuatro H en línea separados por $d$, en base STO-3G. Para cada $d$ buscamos la energía más baja del hamiltoniano electrónico, $E_0(d)$. De la curva salen:

$$R_e = \arg\min_d E_0(d), \qquad D_e = E_0(3{,}0\ \text{Å}) - E_0(R_e).$$

En la zona estirada ($d \ge 1{,}8$ Å) la consigna pide **ganarle a CCSD**, $\lvert E_{\rm VQE} - E_{\rm FCI}\rvert + \sigma < \lvert E_{\rm CCSD} - E_{\rm FCI}\rvert$, y como meta, **precisión química**, $\lvert E_{\rm VQE} - E_{\rm FCI}\rvert < 1{,}6$ mHa. FCI (la solución exacta) solo sirve para comparar.

### 1.2 El hamiltoniano

Hay 4 orbitales moleculares (de Hartree-Fock, con PySCF), cada uno con un electrón $\alpha$ y uno $\beta$: **8 espín-orbitales para 4 electrones**.

$$\hat H = E_{\rm nuc} + \sum_{pq,\sigma} h_{pq}\, a^\dagger_{p\sigma} a_{q\sigma} + \tfrac12 \sum_{pqrs,\sigma\tau} (pq\vert rs)\, a^\dagger_{p\sigma} a^\dagger_{r\tau} a_{s\tau} a_{q\sigma}$$

- $h_{pq}$: energía de **un electrón** (cinética + atracción de los núcleos).
- $(pq\vert rs)$: repulsión entre **dos electrones**.
- $a^\dagger$, $a$: crean y destruyen un electrón en un espín-orbital.

### 1.3 De electrones a qubits (Jordan-Wigner)

Cada espín-orbital es un qubit ($\lvert 1\rangle$ = ocupado). Qubits 0–3 = $\alpha$, 4–7 = $\beta$. Hartree-Fock ocupa los dos orbitales más bajos de cada espín: `|00110011⟩`.

Jordan-Wigner escribe $a_j = Z_0\cdots Z_{j-1}\,(X_j + iY_j)/2$; la cadena de $Z$ da el signo fermiónico. Así $\hat H$ queda como una suma de cadenas de Pauli de 8 qubits:

$$\hat H(d) = \sum_k c_k(d)\, P_k, \qquad P_k \in \{I, X, Y, Z\}^{\otimes 8}.$$

La celda de abajo muestra $\hat H$ a 0,9 Å (equilibrio) y a 2,4 Å (estirada).
