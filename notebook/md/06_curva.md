## 5 · Resultado central: la curva, $R_e$ y $D_e$  *(13:41, Joaquín)*

Con las mejoras de §3 y §4, el flujo para cada distancia es:

$$d \;\to\; \hat H(d)\ \text{(PySCF + Jordan-Wigner)} \;\to\; \boldsymbol\theta^*(d) = \arg\min_{\boldsymbol\theta} E(\boldsymbol\theta)\ \text{(sin ruido; warm start y } \boldsymbol\theta = 0\text{)} \;\to\; E(d) \pm \sigma\ \text{(8 M shots, SamplerV2)} .$$

- **Ansatz:** ADAPT con 15 operadores, el mismo circuito en las 18 distancias (194 CZ con conectividad total).
- **Medición:** 37 grupos QWC (*sorted insertion*), piloto + Neyman, semilla propia por punto.
- **Tabla** (lo que pide la consigna): por punto, la energía medida `E_vqe`, su error `sigma_mHa`, los `shots` y las CZ. Hay tres conteos de CZ: todos con todos (`cz_all`), en línea con el transpiler de Qiskit (`cz_line`) y en línea con nuestro ruteo fSWAP (`cz_line_fswap`, §3.3). Después vienen las referencias clásicas (`E_fci`, `E_ccsd`), los errores en mHa y los dos criterios: `beats_ccsd`, $\lvert E_{\rm VQE} - E_{\rm FCI}\rvert + \sigma < \lvert E_{\rm CCSD} - E_{\rm FCI}\rvert$, que solo se evalúa para $d \ge 1{,}8$ Å, y `chem_acc`, error < 1,6 mHa.

La $\sigma$ es **solo estadística** (la dispersión de los shots, calibrada en §4.3). Sin ruido de compuertas no hay otro sesgo que el del ansatz, ≤ 0,38 mHa (§3).
