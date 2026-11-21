## 3 · Problema 1, demasiadas CZ: ADAPT-VQE

**Idea.** No todas las 26 excitaciones aportan lo mismo. Si el circuito lleva solo las que más bajan la energía, alcanza con muchas menos CZ.

**Cómo funciona ADAPT-VQE** (Grimsley *et al.*, 2019). Arranca de Hartree-Fock y agrega de a un operador:
1. Calcula, para cada excitación, cuánto bajaría la energía al agregarla (su gradiente $g_k$).
2. Agrega la de mayor $\lvert g_k\rvert$ y reoptimiza todos los parámetros.
3. Repite.

**Una sola lista para toda la curva.** En cada paso elegimos el operador con mayor gradiente en *cualquiera* de las 18 distancias. Así el circuito (y sus CZ) es el mismo en todos los puntos; solo cambian los $\theta$.

**Dónde cortar sin usar FCI.** La varianza $\mathrm{Var}(\hat H) = \langle\hat H^2\rangle - \langle\hat H\rangle^2$ vale 0 solo en un autoestado: mide qué tan lejos estamos sin conocer la respuesta. La tabla de abajo la muestra al ir agregando operadores (el error contra FCI va solo como comprobación).
