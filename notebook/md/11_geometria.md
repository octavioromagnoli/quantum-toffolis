## 10 · Otra geometría: H₄ cuadrado y en zigzag  *(12:07, Octavio)*

**Pregunta.** ¿El método depende de que la molécula sea una cadena recta, o anda igual con otra forma? Repetimos el mismo flujo (UCCSD ×2, $\boldsymbol\theta^*$ sin ruido y energía medida con shots) cambiando **solo** dónde están los átomos (`h4_geometry(d, shape)`):

- **Cuadrado de lado $d$.** Por simetría, los dos orbitales del medio son degenerados. RHF tiene que ocupar uno solo y queda con una brecha HOMO-LUMO chica: 0,42 Ha a 0,9 Å, contra 0,76 Ha en la cadena recta. Por eso hay correlación estática fuerte **incluso cerca del equilibrio**: es el caso difícil para los métodos de un determinante.
- **Zigzag:** una cadena con ángulo H-H-H de 120°, una variante de la recta.

Usamos UCCSD ×2 y no ADAPT porque la lista de ADAPT se eligió para la cadena recta. UCCSD es general y su circuito no depende de la geometría (582 CZ, 1648 en línea). Medimos 5 distancias por forma con 2 M de shots por punto.

**Un problema de reproducibilidad que encontramos al revisar.** En el cuadrado, el SCF de RHF a veces termina en una solución **inestable** (un punto silla), según el entorno: a 1,2 Å queda 76 mHa por encima del RHF estable, y con esa referencia CCSD da +180 mHa en vez de −6,5 mHa. Ahora `hamiltonian.molecule` sigue las inestabilidades internas de RHF (`stability()` de PySCF) hasta llegar a una solución estable. En la cadena recta y en el zigzag RHF ya era estable en las 18 distancias, así que la curva principal no cambia. VQE no tenía ese problema: llega a FCI desde cualquiera de las dos referencias, porque FCI no depende de los orbitales.
