## 4 · Problema 2, demasiados shots y una σ dudosa  *(11:46 → 12:40, Santino)*

**Observación.** En lo mínimo, 8 M de shots por punto daban $\sigma \approx$ 0,7 mHa con 46 grupos QWC. Y al repetir 60 veces la medición del mismo circuito, la dispersión de $z = (E - E_{\rm exacta})/\sigma$ dio 1,16 en lugar de 1: o la $\sigma$ subestimaba el error, o era mala suerte.

**Hipótesis.** (a) La varianza total depende de **cómo se agrupan** los términos y de **cómo se reparten** los shots, así que se puede bajar sin tocar el circuito. (b) Con 60 repeticiones, un 1,16 puede ser pura fluctuación estadística.

### 4.1 Cuántos shots hacen falta: $\sigma = \sqrt{K/N}$

Si el grupo $g$ recibe $N_g = f_g N$ shots y su observable tiene varianza $V_g = \langle\hat O_g^2\rangle - \langle\hat O_g\rangle^2$ por shot,

$$\sigma^2 = \sum_g \frac{V_g}{N_g}.$$

Minimizar $\sigma^2$ con $\sum_g N_g = N$ fijo (multiplicadores de Lagrange) da el **reparto de Neyman**,

$$N_g = N\,\frac{\sqrt{V_g}}{\sum_h \sqrt{V_h}} \qquad\Longrightarrow\qquad \sigma^2_{\min} = \frac{\big(\sum_g \sqrt{V_g}\big)^2}{N} \equiv \frac{K}{N}.$$

La $\sigma$ baja como $1/\sqrt N$, y la constante $K$ resume qué tan buena es una estrategia de medición: **para la misma σ, los shots son proporcionales a $K$**. Como las $V_g$ no se conocen de antemano, medimos en **dos etapas**: un piloto con el 10 % de los shots y $N_g \propto \sum_{k\in g}\lvert c_k\rvert$ estima las $V_g$, y el resto se reparte con Neyman. Las muestras de las dos etapas se juntan.

### 4.2 Agrupar mejor: *sorted insertion* y conmutación general

- **Coloreo QWC** (lo mínimo): `group_commuting` de Qiskit arma 46 grupos.
- ***Sorted insertion*** (Crawford *et al.*, *Quantum* 5, 385, 2021): se recorren los términos por $\lvert c_k\rvert$ decreciente, y cada uno va al primer grupo compatible. Quedan 37 grupos, y los términos grandes quedan juntos. Sus covarianzas, muchas veces negativas, se compensan dentro del grupo y bajan $V_g$.
- **Conmutación general (GC)**: términos que conmutan, aunque no qubit a qubit. Quedan 9 grupos. Cada uno se diagonaliza con un Clifford $C$: se toman generadores independientes del grupo (en $\mathrm{GF}(2)$), se arma el estado estabilizador $C\lvert 0\rangle$ que fijan y se mide en $Z$ después de $C^\dagger$, porque $C^\dagger P C = \pm Z\cdots$ para todo $P$ del grupo. **Agrega CZ** al circuito medido.

La celda de abajo calcula $K$ con el estado del ansatz de la curva para las cuatro combinaciones. Es una cuenta para **planificar**: las $\sigma$ que se reportan salen siempre de las muestras.
