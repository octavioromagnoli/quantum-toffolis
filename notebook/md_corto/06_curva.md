## 5 · Resultado central: la curva, $R_e$ y $D_e$

Para cada distancia: $\hat H(d)$ → $\boldsymbol\theta^*$ sin ruido → **$E(d) \pm \sigma$ medida con 8 M de shots**.

- **Ansatz:** ADAPT con 15 operadores, el mismo circuito en las 18 distancias (194 CZ).
- **Medición:** 37 grupos QWC, piloto + Neyman, semilla propia por punto.

**Columnas de la tabla:**
- `E_vqe`, `sigma_mHa` y `shots`: la energía medida, su error y los shots.
- `cz_all`, `cz_line` y `cz_line_fswap`: las CZ todos con todos, en línea con el transpiler y en línea con nuestro ruteo.
- `beats_ccsd`: le gana a CCSD contando σ (solo para $d \ge 1{,}8$ Å). `chem_acc`: error < 1,6 mHa.

La σ es solo estadística. Sin ruido, el único sesgo es el del ansatz (≤ 0,38 mHa).
