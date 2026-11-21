## 7 · Qué decidimos, qué no funcionó y limitaciones

| decisión | por qué | alternativa descartada |
|---|---|---|
| ADAPT-VQE fermiónico, lista común de 15 operadores | 194 CZ (378 en línea con fSWAP) y error sin ruido ≤ 0,38 mHa; el mismo circuito en toda la curva | UCCSD ×2: exacto, pero 582 CZ (1648 en línea). Compacto de 5 qubits (16 CZ): no conserva N ni S² por construcción y no escala; queda para FakeFez |
| medir el circuito sin rutear con Ĥ original | es el mismo estado que el ruteado; Ĥ reordenado necesita ~1,9× más grupos QWC | medir el circuito ruteado (las 378 CZ en línea se reportan aparte) |
| optimizar sin ruido | permitido; mucho más rápido | SPSA con shots: lento y ruidoso |
| SamplerV2 + σ propia | barras honestas, permite Neyman | EstimatorV2 con `precision`: σ fija |
| R_e con cuártico en 7 puntos (0,7–1,1 Å) | sesgo < 1 mÅ al aplicarlo a FCI en la misma grilla | cuadrático en 5 puntos: sesgo de +8 mÅ, mayor que σ |
| referencia R_e(FCI) minimizando PySCF | compara el VQE contra la verdad, no contra otro ajuste | ajustar FCI con el mismo modelo (esconde el sesgo del ajuste) |
| varianza por eco de Loschmidt | cuenta eventos raros: resuelve Var ~ 10⁻⁶ Ha² con 1M shots | solo ⟨H²⟩ por grupos QWC: σ ≈ 3·10⁻³ Ha² (no separa ×1 de ×2) |

**Qué no funcionó:**
- Semillas consecutivas en el Sampler ⇒ σ subestimada (corregido).
- La misma semilla en las 18 distancias ⇒ ruido de shots correlacionado entre puntos vecinos (detectado por Gonzalo: ~3σ seguidos en 0,7–0,9 Å con 200 k shots). Eso sesga el ajuste de R_e y hace que su σ de Monte Carlo, que supone errores independientes, quede chica. Ahora cada punto y cada diagnóstico usa su propia semilla, separadas por un paso de 1000 para que no choquen con las semillas derivadas del estimador (seed+1 de la etapa de Neyman).
- Medir la varianza como ⟨H²⟩ − ⟨H⟩² con grupos QWC: (H−c)² tiene 1775 Paulis en 225 grupos y cada grupo fluctúa mucho aunque el total sea ≈ 0, así que la σ queda en ~10⁻³ Ha² aunque se usen millones de shots. Restar la constante de H (c) reduce σ a la mitad, pero no alcanza. Lo dejamos como medición directa y sumamos el eco.
- Eco con 2 pasos de Trotter: el error de Trotter (1,3·10⁻⁶) se confundía con la señal. Con 4 pasos baja a 9·10⁻⁸.

> TODO (todos): agregar lo que vayan probando.
