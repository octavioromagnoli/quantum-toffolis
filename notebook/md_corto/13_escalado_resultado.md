**Qué dio:**

| | crece como | para FeMoco |
|---|---|---|
| términos, grupos, CZ | ~$N^4$ (polinomial) | ~10⁷ CZ, ~780 mil parámetros |
| shots por evaluación de energía | ~$N^5$ | ~2–3·10¹¹ (y un VQE necesita miles de evaluaciones) |
| circuito sin errores (0,1 %/CZ) | cae exponencialmente | $10^{-4850}$ |
| FCI clásico | exponencial | ~4·10³⁰ determinantes |

- **CCSD empeora con el tamaño:** su error en la zona estirada pasa de 18 mHa (H₄) a 242 mHa (H₁₂).
- **Conclusión:** VQE escala en compuertas, pero las mediciones y el ruido lo vuelven impracticable para FeMoco sin corrección de errores. Las estimaciones serias usan estimación de fase con corrección de errores (~4 millones de qubits físicos, Lee *et al.*, 2021).
- Lo que **sí escala** de nuestro trabajo: argumentar la confiabilidad sin FCI y contar recursos antes de gastar shots.

**Limitación:** H$_N$ no es FeMoco (no tiene metales ni orbitales $d$). La extrapolación sirve para contar recursos, no para predecir la precisión.
