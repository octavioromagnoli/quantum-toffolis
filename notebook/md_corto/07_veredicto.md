### 6.1 Qué dicen las mediciones

| prueba | resultado |
|---|---|
| Límite físico (3,0 Å) | VQE ~1 mHa debajo de 4 H separados (razonable). CCSD **24 mHa** debajo: imposible a esa distancia. |
| Varianza por eco | curva: 2·10⁻⁴ a 5·10⁻⁶ Ha². ADAPT 12: 7–18× más. HF: 0,11–0,16 Ha². |
| Simetrías | $N_\alpha = N_\beta = 2$ en **todos** los shots; $\lvert S^2\rvert \le 0{,}011 \pm 0{,}012$ (singulete) |
| Profundidad ↓ (12 → 15 ops) | la energía baja hasta 2,4 ± 0,8 mHa: los últimos operadores hacen falta |
| Profundidad ↑ (15 → UCCSD ×2) | cambio < 1,6σ con 3× más CZ: no vale la pena |

**Conclusión sin FCI:** la curva usa el estado variacional más bajo que encontramos, con las simetrías correctas y una varianza 500 veces menor que la de HF, y profundizar no lo mejora. **Recién ahora** comparamos con FCI (última tabla de arriba): errores de 0,05–0,56 mHa, contra 9–31 mHa de CCSD.
