# Romper una molécula, un qubit a la vez — H₄ con VQE

**Quantum Jam 2026 · equipo quantum-toffolis** (Octavio, Gonzalo, Santino, Joaquín)

Calculamos con **VQE** la curva de energía de la cadena lineal H₄ (STO-3G, 8 qubits con Jordan-Wigner) entre 0,6 y 3,0 Å y la comparamos con CCSD y FCI.

**Resumen de la estrategia**
1. Hamiltoniano por punto con PySCF → Jordan-Wigner (185 cadenas de Pauli).
2. Ansatz químico que conserva N_α y N_β, desde Hartree-Fock: **ADAPT-VQE** fermiónico, con una lista común de 15 excitaciones para toda la curva (194 CZ; 378 con qubits en línea usando un ruteo fSWAP propio).
3. Parámetros buscados sin ruido con gradiente analítico (permitido por la consigna).
4. **Energía reportada = medición con `SamplerV2` de `qiskit_aer`**: grupos QWC, rotación de base, shots repartidos en dos etapas (piloto + Neyman) y σ calculada desde las muestras.
5. R_e por ajuste cuártico en 7 puntos, elegido por su sesgo contra el mínimo FCI verdadero; incertidumbre por Monte Carlo. D_e = E(3,0 Å) − E(R_e).
6. Confiabilidad sin FCI, todo **medido con shots**: varianza de la energía (por grupos y por eco de Loschmidt), simetrías N_α, N_β y S², y convergencia con la profundidad.

**Resultados** (ADAPT 15 operadores, 8M shots por punto, 194 CZ todos con todos):
- Precisión química en las 18 distancias (error máximo 1,15 mHa, a 0,7 Å). En la zona estirada (d ≥ 1,8 Å), error ≤ 0,56 mHa con σ ≈ 0,5 mHa, contra 9–31 mHa de CCSD: le ganamos en 7/7 puntos.
- R_e = 0,8873 ± 0,0018 Å (FCI: 0,8886 Å); D_e = 313,00 ± 0,66 mHa (FCI: 313,21 mHa).
- Varianza medida por eco ≤ 2·10⁻⁴ Ha² en 1,8, 2,4 y 3,0 Å (HF: 0,11–0,16 Ha²), y UCCSD ×2, con 3 veces más CZ, no baja la energía medida más allá de 1,6σ.
- Contra UCCSD ×2 (582 CZ, 1648 en línea) son 3 veces menos CZ con conectividad total y 4,4 veces menos en línea.
