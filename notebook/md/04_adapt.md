## 3 · Problema 1, demasiadas CZ: ADAPT-VQE  *(11:46 → 13:41, Gonzalo)*

**Observación.** UCCSD ×2 llega a FCI, pero con 582 CZ (1648 en línea). Con 0,1 % de error por CZ eso no sirve, y el jurado mira las CZ.

**Hipótesis.** No todas las 26 excitaciones aportan lo mismo. Si agregamos al circuito solo las que más bajan la energía, deberíamos llegar a precisión química con muchas menos CZ.

### 3.1 Cómo funciona ADAPT-VQE

ADAPT-VQE (Grimsley *et al.*, *Nat. Commun.* 10, 3007, 2019) arma el ansatz de a un operador. Arranca de $\lvert\psi^{(0)}\rangle = \lvert \mathrm{HF}\rangle$ y, en el paso $n$, mira cuánto bajaría la energía si agregara cada excitación $\hat A_k = \hat\tau_k - \hat\tau_k^\dagger$ del *pool* (las 26 de UCCSD, que se pueden repetir):

$$g_k = \frac{\partial}{\partial\theta}\,\big\langle\psi^{(n)}\big\rvert\, e^{-\theta\hat A_k}\,\hat H\, e^{\theta\hat A_k}\big\lvert\psi^{(n)}\big\rangle\Big\rvert_{\theta=0} = \big\langle\psi^{(n)}\big\rvert\,[\hat H, \hat A_k]\,\big\lvert\psi^{(n)}\big\rangle .$$

Agrega el de $\lvert g_k\rvert$ máximo, $\lvert\psi^{(n+1)}\rangle = e^{\theta_{n+1}\hat A_{k^*}}\lvert\psi^{(n)}\rangle$, reoptimiza **todos** los $\theta$ y repite. El gradiente es exacto y barato con el simulador (`adapt.py`).

**Una lista común para toda la curva.** En vez de un circuito distinto por distancia, en cada paso elegimos el operador con mayor gradiente en *cualquiera* de las 18 distancias,

$$k^* = \arg\max_k\ \max_{d}\ \lvert g_k(d)\rvert ,$$

y reoptimizamos cada distancia por separado. Así el circuito, y por lo tanto la cantidad de CZ, es el mismo en todos los puntos; solo cambian los $\theta$.

**Dónde cortar sin usar FCI.** La varianza de la energía,

$$\mathrm{Var}(\hat H) = \langle\hat H^2\rangle - \langle\hat H\rangle^2 \;\ge\; 0,$$

vale 0 solo en un autoestado, así que mide qué tan lejos estamos de uno sin conocer $E_0$. La tabla de abajo la calcula con el estado del simulador (solo para elegir el corte; en §6 la **medimos** con shots) al ir agregando operadores en el orden en que ADAPT los eligió. El error contra FCI aparece solo como comprobación.
