## 3 · Medir la energía con shots  *(Santino)*

- 185 términos → 46 grupos que conmutan qubit a qubit (QWC). Un circuito por grupo.
- Cada shot da una muestra de O_g = Σ c_k P_k ⇒ E_g = media, σ_g² = var/N_g (covarianzas incluidas).
- **Dos etapas**: 10 % piloto con N_g ∝ Σ|c_k| → estimamos σ_g → resto con Neyman N_g ∝ σ_g.
- Semilla independiente por grupo (`SeedSequence`): con semillas consecutivas de Aer vimos errores correlacionados entre grupos y desvíos de 4σ.
- Calibración: 60 repeticiones a 2,0 Å ⇒ media de z ≈ 0, dispersión de z ≈ 1,16.

> TODO (Santino): agrupamiento general (no QWC) con diagonalización Clifford, bootstrap de σ, estudio de calibración con ≥200 repeticiones, curva σ vs shots, y presupuesto total.
