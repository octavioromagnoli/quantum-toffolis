## 10 · Otra geometría: H₄ cuadrado y en zigzag

¿El método depende de que la molécula sea una cadena recta? Repetimos el flujo (UCCSD ×2, 2 M de shots por punto) cambiando solo la forma:

- **Cuadrado:** dos orbitales degenerados; hay correlación fuerte **incluso cerca del equilibrio**. Es el caso difícil para CCSD.
- **Zigzag:** cadena con ángulo de 120°.

Usamos UCCSD y no ADAPT porque la lista de ADAPT se eligió para la cadena recta.

**Un problema que encontramos:** en el cuadrado, Hartree-Fock a veces terminaba en una solución inestable según la máquina, y CCSD daba +180 mHa en vez de −6,5. Ahora `molecule` busca siempre la solución estable. A VQE no le afectaba.
