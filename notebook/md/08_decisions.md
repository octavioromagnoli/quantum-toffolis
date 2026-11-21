## 7 · Qué decidimos, qué no funcionó y limitaciones

| decisión | por qué | alternativa descartada |
|---|---|---|
| UCCSD ×2 | exacto en toda la curva | hardware-efficient: difícil de optimizar, rompe N |
| optimizar sin ruido | permitido; mucho más rápido | SPSA con shots: lento y ruidoso |
| SamplerV2 + σ propia | barras honestas, permite Neyman | EstimatorV2 con `precision`: σ fija |

**Qué no funcionó:** semillas consecutivas en el Sampler ⇒ σ subestimada (corregido).

> TODO (todos): agregar lo que vayan probando.
