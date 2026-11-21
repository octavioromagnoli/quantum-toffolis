# Romper una molécula, un qubit a la vez — H₄ con VQE

**Quantum Jam 2026 · equipo quantum-toffolis** (Octavio, Gonzalo, Santino, Joaquín)

Calculamos con **VQE** la curva de energía de la cadena lineal H₄ (STO-3G, 8 qubits con Jordan-Wigner) entre 0,6 y 3,0 Å y la comparamos con CCSD y FCI.

**Resumen de la estrategia**
1. Hamiltoniano por punto con PySCF → Jordan-Wigner (185 cadenas de Pauli).
2. Ansatz químico que conserva N_α, N_β (UCCSD, 2 pasos de Trotter) desde Hartree-Fock.
3. Parámetros buscados sin ruido con gradiente analítico (permitido por la consigna).
4. **Energía reportada = medición con `SamplerV2` de `qiskit_aer`**: grupos QWC, rotación de base, shots repartidos en dos etapas (piloto + Neyman) y σ calculada desde las muestras.
5. R_e por ajuste cuártico en 7 puntos, elegido por su sesgo contra el mínimo FCI verdadero; incertidumbre por Monte Carlo. D_e = E(3,0 Å) − E(R_e).
6. Confiabilidad sin FCI, todo **medido con shots**: varianza de la energía (por grupos y por eco de Loschmidt), simetrías N_α, N_β y S², y convergencia con la profundidad.

**Resultados** (UCCSD ×2, 8M shots por punto, 582 CZ todos con todos):
- Precisión química en las 18 distancias; en la zona estirada (d ≥ 1,8 Å), error ≤ 1,45 mHa contra 9–31 mHa de CCSD.
- R_e = 0,8845 ± 0,0023 Å (FCI: 0,8886 Å); D_e = 313,92 ± 0,85 mHa (FCI: 313,21 mHa).
- Varianza medida por eco ≤ 5·10⁻⁶ Ha² en 1,8, 2,4 y 3,0 Å (HF: 0,11–0,16 Ha²).

> TODO (Joaquín): actualizar al cerrar con el ansatz y los shots finales.
