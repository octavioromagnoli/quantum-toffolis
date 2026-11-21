## 6 · Problema 4: ¿por qué confiar en la zona estirada sin mirar FCI?

La consigna pide argumentar la confiabilidad **antes** de comparar con FCI. Antes lo hacíamos con el estado del simulador; ahora **todo se mide con shots**.

En 1,8, 2,4 y 3,0 Å comparamos cuatro estados: **HF**, **ADAPT con 12 operadores** (uno menos profundo), **la curva** (ADAPT 15) y **UCCSD ×2** (mucho más profundo). Seis pruebas:

1. **Principio variacional:** VQE nunca queda por debajo de la energía exacta. Entre estados medidos, el más bajo es el mejor.
2. **Límite físico:** a 3,0 Å los átomos casi no interactúan, así que la energía no puede estar muy por debajo de la de 4 H separados (−1,866327 Ha, un número exacto y barato).
3. **Varianza por grupos:** en un autoestado, $\mathrm{Var}(\hat H) = 0$. Se mide igual que la energía, pero la resta de números grandes deja σ ≈ 2·10⁻³ Ha²: separa HF del resto, no más.
4. **Varianza por eco de Loschmidt:** preparamos $\psi$, lo evolucionamos con $e^{-i\hat Ht}$ y lo deshacemos. La probabilidad de no volver al inicio es $\approx \mathrm{Var}(\hat H)\,t^2$. Cuenta **eventos raros**, así que resuelve varianzas de ~10⁻⁶ Ha². Usa más de 4000 CZ: es un diagnóstico de simulador, no para hardware.
5. **Simetrías medidas:** el estado buscado tiene $N_\alpha = N_\beta = 2$ y $S^2 = 0$ (singulete). Al estirar, singulete y triplete casi se tocan, así que una energía "buena" puede ser el estado equivocado: la $S^2$ medida lo delata.
6. **Profundidad:** si la curva convergió, un ansatz menos profundo da energía más alta, y uno mucho más profundo no la baja más allá de σ.
