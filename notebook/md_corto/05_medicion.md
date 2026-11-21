## 4 · Problema 2, demasiados shots y una σ dudosa

**Qué vimos.** 8 M de shots por punto daban σ ≈ 0,7 mHa. Y al repetir 60 veces la medición, la dispersión de $z = (E - E_{\rm exacta})/\sigma$ dio 1,16 en vez de 1: ¿la σ subestima el error o fue mala suerte?

### 4.1 Cuántos shots hacen falta

Repartiendo bien los shots entre grupos (reparto de Neyman, $N_g \propto \sqrt{V_g}$), la σ queda

$$\sigma = \sqrt{K/N}, \qquad K = \Big(\sum_g \sqrt{V_g}\Big)^2 .$$

$K$ resume qué tan buena es la estrategia: **para la misma σ, los shots son proporcionales a $K$**. Como las varianzas $V_g$ no se conocen de antemano, medimos en dos etapas: un piloto con el 10 % de los shots las estima y el resto se reparte con Neyman.

### 4.2 Agrupar mejor

- **Coloreo QWC** (lo mínimo): 46 grupos.
- ***Sorted insertion***: ordena los términos por tamaño y los mete en el primer grupo compatible. **37 grupos**, sin compuertas extra.
- **Conmutación general (GC)**: 9 grupos, pero cada uno se mide con un circuito Clifford que **agrega CZ**.

La celda de abajo calcula $K$ para cada opción (solo para planificar: las σ reportadas salen de las muestras).
