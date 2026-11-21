### Extra: otra geometría (H₄ cuadrado y en zigzag)  *(Octavio)*

Repetimos el mismo flujo (ansatz UCCSD×2, búsqueda de θ* sin ruido y energía medida con shots) cambiando solo dónde están los átomos (`h4_geometry(d, shape)`):

- **Cuadrado de lado d.** Por simetría, los dos orbitales del medio son degenerados. RHF tiene que ocupar uno solo y queda con una brecha HOMO-LUMO chica: 0,42 Ha a 0,9 Å, contra 0,76 Ha en la cadena recta. Por eso hay correlación estática fuerte **incluso cerca del equilibrio**.
- **Zigzag:** una cadena con ángulo H-H-H de 120°, una variante de la cadena recta.

Medimos 5 distancias por geometría con 2 M de shots por punto (σ ≈ 1,4 mHa). Resultados (gráfico y tabla de abajo):

- **Cuadrado.** HF ya se equivoca en 124 mHa a 0,7 Å (en zigzag, en 37 mHa). CCSD queda 9,4 mHa **por encima** de FCI a 0,7 Å, 18,9 mHa por debajo a 1,8 Å, y a 3,0 Å no converge, así que su número ahí no sirve. A 0,7 Å el determinante HF tiene peso nulo en el estado exacto, y aun así VQE llega a FCI, porque las excitaciones dobles de UCCSD lo conectan con los determinantes que sí pesan.
- **Zigzag.** Se comporta como la cadena recta: CCSD es exacto cerca del equilibrio y falla al estirar, hasta 31,6 mHa por debajo de FCI a 2,4 Å.
- **VQE.** Sin ruido llega a FCI en todos los puntos de las dos geometrías (error < 0,001 mHa). Con la energía medida, todos los puntos quedan dentro de 2σ de FCI. En la zona estirada (d ≥ 1,8 Å) le gana a CCSD en los 6 puntos con el criterio |E_VQE − E_FCI| + σ < |E_CCSD − E_FCI|. En el cuadrado a 3,0 Å el margen es chico (1,7 contra 1,8 mHa) y CCSD no convergió, así que ese punto no cuenta como victoria.
- **Costo.** Es el mismo circuito que en la cadena recta (582 CZ con conectividad total, 1648 en línea), porque solo depende de la cantidad de orbitales y de electrones.

**Limitación:** con 2 M de shots, σ ≈ 1,4 mHa no alcanza para afirmar precisión química punto por punto. Para eso harían falta los ~8 M de shots de la curva principal.
