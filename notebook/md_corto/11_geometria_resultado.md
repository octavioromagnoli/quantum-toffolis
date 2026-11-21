**Qué dio** (gráfico y tabla de arriba):

- **Cuadrado:** a 0,7 Å, Hartree-Fock pesa **cero** en el estado exacto. CCSD queda 9,4 mHa por encima de FCI y no converge a 3,0 Å.
- **Zigzag:** igual que la cadena recta; CCSD falla al estirar (hasta −31,6 mHa).
- **VQE:** llega a FCI sin ruido en todos los puntos, **incluso donde HF pesa cero**. Con shots le gana a CCSD en la zona estirada donde CCSD converge.

**Por qué:** VQE no necesita que HF sea buena aproximación, solo que el ansatz pueda llegar al estado exacto. CCSD sí lo necesita.

**Limitación:** con 2 M de shots (σ ≈ 0,6–1,1 mHa) no alcanza para afirmar precisión química punto por punto.
