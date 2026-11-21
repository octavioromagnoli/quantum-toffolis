## 12 · ¿Y la nitrogenasa? Cuánto costaría la molécula real

La consigna arranca con la nitrogenasa y su cofactor **FeMoco**. Para ver hasta dónde llega el método, medimos cómo crecen los recursos al alargar la cadena (H₄, H₆, …, H₁₄) y extrapolamos a FeMoco (54 electrones en 54 orbitales, 108 qubits).

**Qué contamos** para cada cadena H$_N$ (`scaling.py`): términos de $\hat H$ y grupos, parámetros y CZ de UCCSD, shots para σ = 1,6 mHa, probabilidad de un circuito sin errores, y el error de CCSD contra FCI.

Lo pesado (N = 10 a 14, ~1 h) está precalculado. La celda de abajo recalcula N = 4, 6 y 8 para mostrar que es reproducible, ajusta leyes de potencia y extrapola.
