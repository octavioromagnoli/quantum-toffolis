### 1.6 Por qué CCSD falla al estirar

En el gráfico: el error de CCSD (izquierda) y el peso de Hartree-Fock en el estado exacto (derecha).

- **Cerca del equilibrio** el estado exacto es casi HF (peso 0,95) y CCSD es prácticamente exacto.
- **Al estirar** el estado exacto es una mezcla de muchas configuraciones (peso de HF: 0,58 a 1,8 Å, 0,27 a 3,0 Å). CCSD solo incluye bien las excitaciones simples y dobles a partir de HF, y deja de alcanzar.
- **CCSD no es variacional:** su energía no es la de ningún estado, así que puede quedar *por debajo* de la exacta. Acá llega a −30,6 mHa (2,6 Å), unas 19 veces la precisión química.

**Por qué VQE le puede ganar.** VQE mide la energía de un estado real, así que nunca queda por debajo de la exacta (salvo el ruido de los shots, cubierto por $\pm\sigma$):

$$E(\boldsymbol\theta) = \langle\psi(\boldsymbol\theta)\rvert\hat H\lvert\psi(\boldsymbol\theta)\rangle \ge E_0 .$$

Entonces, entre dos estados medidos, **el de energía más baja es el mejor**. Esa es la base de la confiabilidad (§6).
