## 13 · Resumen final: de lo mínimo a nuestro resultado

La celda de abajo **no recalcula nada**: junta los números que salieron de las celdas anteriores en tres tablas y un gráfico.

1. **Nuestros mejores resultados**, uno por fila.
2. **Lo mínimo contra lo final**: el mismo problema resuelto con el VQE "de libro" de las 11:09 (UCCSD ×2, coloreo QWC, ajuste cuadrático, diagnósticos con el estado) y con la versión final. Los valores de la corrida de las 11:09 que ya no se recalculan están marcados como tales; el resto sale de este notebook.
3. **La evolución con ruido**: cuántas CZ usa cada etapa con qubits en línea y qué error deja el ruido a 2,4 Å. Es el sesgo exacto con matriz densidad de §7.1 y §8, así que las etapas se comparan sin ruido de shots.
