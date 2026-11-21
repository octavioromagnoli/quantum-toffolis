## 5 · ¿Por qué confiar en la zona estirada?  *(Joaquín)*

La consigna pide argumentar la confiabilidad **antes** de comparar con FCI. Entonces todo lo de esta sección se **mide con shots** (SamplerV2 de Aer) y no usa FCI. FCI aparece recién en la última tabla, como comprobación.

En tres puntos estirados (1,8, 2,4 y 3,0 Å) comparamos cuatro estados preparados con circuitos: **HF** (θ = 0), **ADAPT con 12 operadores** (el paso anterior al corte; 169 CZ), el **ansatz de la curva** (ADAPT, 15 operadores; 194 CZ) y **UCCSD ×2** (52 parámetros; 582 CZ), mucho más profundo.

**1. Principio variacional.** Para cualquier estado, E(θ) ≥ E₀, salvo ruido estadístico (±σ). Por eso, entre estados variacionales, el que da la energía más baja es el mejor. CCSD no es variacional: en la zona estirada queda *por debajo* de la exacta, y ahí "más bajo" no quiere decir "mejor".

**2. Varianza de la energía medida por grupos.** En un autoestado, Var(H) = ⟨H²⟩ − ⟨H⟩² = 0. La medimos como ⟨(H−c)²⟩ − (⟨H⟩−c)², con ⟨(H−c)²⟩ medido con el mismo estimador que la energía: (H−c)² tiene 1775 cadenas de Pauli, que se miden en 225 grupos QWC. Tomamos c = término constante de H porque reduce σ a la mitad respecto de c = E.
*Resultado y limitación:* con 2M de shots, σ ≈ 1,5–2·10⁻³ Ha². HF da 0,12, 0,15 y 0,18 Ha² (1,8, 2,4 y 3,0 Å). Los tres ansätze dan entre 0 y 2·10⁻³ Ha², compatibles con 0 y al menos 50 veces menos que HF. Como es una resta de números grandes, no separa los ansätze entre sí, y por ruido puede dar negativa: en la tabla se muestra recortada en 0, con su σ. Esta medición solo separa HF del resto; la que resuelve es el eco.

**3. Varianza por eco de Loschmidt.** Es la forma fina de medir lo mismo. El circuito prepara ψ, lo evoluciona con e^{−iHt} y lo des-prepara con U(θ)†:

  P(no volver a |0…0⟩) = 1 − |⟨ψ|e^{−iHt}|ψ⟩|² = Var(H)·t² + O(t⁴).

En vez de restar dos números grandes, cuenta **eventos raros**: si ψ es un autoestado, no vuelve ningún shot "mal". Con N shots y 0 eventos, la cota del 95 % es Var ≲ 3/(N t²). Usamos t = 1 Ha⁻¹ y Suzuki-Trotter de 2º orden con 4 pasos. El error de Trotter lo verificamos sin FCI, comparando 2 contra 4 pasos: escala como 1/n⁴ y con 4 pasos es 9·10⁻⁸, debajo de la resolución de 10⁶ shots. A orden dominante, el eco *subestima* Var (Mandelstam-Tamm); a t = 1 la diferencia es ~8 % para UCCSD ×1.
**Resultado:** con 10⁶ shots, el ansatz de la curva da Var = 1,9·10⁻⁴ (1,8 Å), 9·10⁻⁵ (2,4 Å) y 5·10⁻⁶ Ha² (3,0 Å), igual a la calculada con el estado. ADAPT con 12 operadores da entre 7 y 18 veces más (2,9·10⁻³, 6,2·10⁻⁴ y 8,9·10⁻⁵ Ha²), y UCCSD ×2 no tiene ningún evento (Var ≤ 3·10⁻⁶ Ha²). HF da 0,11–0,16 Ha².
*Costo:* el circuito de eco tiene U(θ), la evolución con Trotter y U(θ)†: entre 4100 y 5700 CZ según el estado (columna `echo_cz`). **Es un diagnóstico de simulador sin ruido, no una propuesta para hardware**: con 0,1 % de error por CZ, el ruido taparía la señal. En hardware, el diagnóstico equivalente es la varianza por grupos, con más shots.
*Qué implica:* por la cota de Weinstein, hay un autovalor exacto de Ĥ a menos de √Var de la energía medida: ≤ 15 mHa a 1,8 Å y ≤ 2,3 mHa a 3,0 Å. Esa cota es floja; la que importa es la de segundo orden: si hay un gap Δ al primer excitado del mismo sector, el error es E − E₀ ≈ Var/Δ. Con Var ~10⁻⁴ Ha², ese error queda debajo de la σ de shots (0,5 mHa) si Δ ≳ 0,2 Ha. No medimos Δ; el punto 5 lo chequea de otra forma: un ansatz mucho más profundo no baja la energía. Que ese autovalor sea el fundamental lo sostienen los puntos 1, 4 y 5: es el estado variacional más bajo que encontramos, tiene las simetrías correctas y no baja más al agregar profundidad.

**4. Simetrías medidas.** H₄ tiene un fundamental singlete con N_α = N_β = 2. Medimos ⟨N_α⟩, ⟨N_β⟩ y ⟨S²⟩ con shots:
- **N_α y N_β dan exactamente 2 en todos los shots**: el ansatz de excitaciones conserva el número de partículas.
- **S² es la prueba que discrimina estados de espín mezclado.** Al estirar, el singlete y el triplete se acercan. UCCSD ×1 (un solo paso de Trotter, que no está adaptado al espín) converge a una mezcla, con ⟨S²⟩ medido de 0,64 a 2,4 Å y 0,85 a 3,0 Å, y a 3,0 Å su energía medida no se distingue de la buena: la energía sola no lo delata, la S² medida sí. Por eso la lista de ADAPT completa los pares de espín de las excitaciones simples (sección 2). El ansatz de la curva da |S²| ≤ 0,011 ± 0,012, compatible con 0, igual que HF y UCCSD ×2. Un triplete puro da 2.

**5. Convergencia con la profundidad.** Dos comparaciones, con energías medidas:
- **Hacia abajo** (ADAPT 12 → 15 operadores): la energía baja 2,4 ± 0,8 mHa a 1,8 Å y 1,4 ± 0,7 mHa a 2,4 Å; a 3,0 Å no cambia dentro del ruido (−0,6 ± 0,7 mHa). La varianza por eco cae entre 7 y 18 veces en los tres puntos. Los operadores 13 a 15 hacen falta.
- **Hacia arriba** (15 operadores → UCCSD ×2, con 3 veces más CZ): la energía medida cambia 0,6 ± 0,8, −0,1 ± 0,7 y 1,1 ± 0,7 mHa (1,8, 2,4 y 3,0 Å), todo dentro de 1,6σ. Profundizar mucho más no mejora nada que los shots puedan resolver.

La varianza de UCCSD ×2 sí es menor (≤ 3·10⁻⁶ contra ~10⁻⁴ Ha²): ADAPT no es un autoestado exacto. Pero la energía medida no lo distingue de UCCSD ×2, y por eso elegimos el circuito con menos CZ.

**6. Barras de error calibradas.** Ver la calibración de σ en la sección de medición (Santino).

> La tabla siguiente reúne todo lo medido. La última celda compara recién ahí con FCI.
