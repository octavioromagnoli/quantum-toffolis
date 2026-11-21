## 5 · ¿Por qué confiar en la zona estirada?  *(Joaquín)*

La consigna pide argumentar la confiabilidad **antes** de comparar con FCI. Entonces todo lo de esta sección se **mide con shots** (SamplerV2 de Aer) y no usa FCI. FCI aparece recién en la última tabla, como comprobación.

En tres puntos estirados (1,8, 2,4 y 3,0 Å) comparamos tres estados preparados con circuitos: **HF** (θ = 0), **UCCSD ×1** (26 parámetros) y el **ansatz de la curva** (UCCSD ×2, 52 parámetros).

**1. Principio variacional.** Para cualquier estado, E(θ) ≥ E₀, salvo ruido estadístico (±σ). Por eso, entre estados variacionales, el que da la energía más baja es el mejor. CCSD no es variacional: en la zona estirada queda *por debajo* de la exacta, y ahí "más bajo" no quiere decir "mejor".

**2. Varianza de la energía medida por grupos.** En un autoestado, Var(H) = ⟨H²⟩ − ⟨H⟩² = 0. La medimos como ⟨(H−c)²⟩ − (⟨H⟩−c)², con ⟨(H−c)²⟩ medido con el mismo estimador que la energía: (H−c)² tiene 1775 cadenas de Pauli, que se miden en 225 grupos QWC. Tomamos c = término constante de H porque reduce σ a la mitad respecto de c = E.
*Resultado y limitación:* con 2M de shots, σ ≈ 2·10⁻³ Ha². HF da 0,111, 0,149 y 0,176 ± 0,002 Ha² (1,8, 2,4 y 3,0 Å). El ansatz de la curva da (0,5 / 0,2 / 0,9) · 10⁻³ ± 2·10⁻³ Ha², compatible con 0 y al menos 50 veces menos que HF. Como es una resta de números grandes, no llega a separar UCCSD ×1 (~10⁻⁴–10⁻³) de ×2. Por eso sumamos el eco.

**3. Varianza por eco de Loschmidt.** Es la forma fina de medir lo mismo. El circuito prepara ψ, lo evoluciona con e^{−iHt} y lo des-prepara con U(θ)†:

  P(no volver a |0…0⟩) = 1 − |⟨ψ|e^{−iHt}|ψ⟩|² = Var(H)·t² + O(t⁴).

En vez de restar dos números grandes, cuenta **eventos raros**: si ψ es un autoestado, no vuelve ningún shot "mal". Con N shots y 0 eventos, la cota del 95 % es Var ≲ 3/(N t²). Usamos t = 1 Ha⁻¹ y Suzuki-Trotter de 2º orden con 4 pasos. El error de Trotter lo verificamos sin FCI, comparando 2 contra 4 pasos: escala como 1/n⁴ y con 4 pasos es 9·10⁻⁸, debajo de la resolución de 10⁶ shots. A orden dominante, el eco *subestima* Var (Mandelstam-Tamm); a t = 1 la diferencia es ~8 % para UCCSD ×1.
**Resultado:** con 10⁶ shots, el ansatz de la curva tuvo **1 evento a 1,8 Å y ninguno a 2,4 ni a 3,0 Å**, así que Var ≤ 5·10⁻⁶ Ha² en los tres puntos (95 %). Para comparar: HF da 0,11–0,16 Ha² y UCCSD ×1 da 2,8·10⁻³ (1,8 Å), 3,1·10⁻⁴ (2,4 Å) y 1,0·10⁻⁴ Ha² (3,0 Å).
*Qué implica:* por la cota de Weinstein, hay un autovalor exacto de Ĥ a menos de √Var ≤ 2,2 mHa de la energía medida. Si además hay un gap Δ al primer excitado del mismo sector, el error es E − E₀ ≈ Var/Δ, mucho menor que la σ de shots. Que ese autovalor sea el fundamental lo sostienen los puntos 1, 4 y 5: es el estado variacional más bajo que encontramos, tiene las simetrías correctas y no baja más al agregar profundidad.

**4. Simetrías medidas.** H₄ tiene un fundamental singlete con N_α = N_β = 2. Medimos ⟨N_α⟩, ⟨N_β⟩ y ⟨S²⟩ con shots:
- **N_α y N_β dan exactamente 2 en todos los shots**: el ansatz de excitaciones conserva el número de partículas.
- **S² es la prueba que discrimina.** Al estirar, el singlete y el triplete se acercan, y UCCSD ×1 (un solo paso de Trotter, que no está adaptado al espín) converge a una mezcla. Su ⟨S²⟩ medido da 0,64 ± 0,01 a 2,4 Å y 0,85 ± 0,01 a 3,0 Å. A 3,0 Å su energía medida es **indistinguible** de la del ansatz de la curva (−0,2 ± 0,9 mHa): la energía sola no lo delata, la S² medida sí. El ansatz de la curva da |S²| ≤ 0,04 con σ = 0,016 (compatible con 0). Un singlete da 0 y un triplete puro da 2; HF da entre −0,02 y 0,00 ± 0,011, así que la σ alcanza de sobra para separar los casos.

**5. Convergencia con la profundidad.** Al pasar de UCCSD ×1 a ×2, la energía medida baja 3,3 ± 1,0 mHa a 1,8 Å y 1,6 ± 0,9 mHa a 2,4 Å, y no cambia dentro del ruido a 3,0 Å (−0,2 ± 0,9 mHa). La varianza medida por eco cae de 10⁻⁴–10⁻³ Ha² a ≤ 5·10⁻⁶ Ha² en los tres. Por eso la energía sola no alcanza como criterio de convergencia: a 3,0 Å es la varianza (y la S²) la que separa los dos estados. Una varianza ≈ 0 es la firma de un autoestado: si se agregan más parámetros, la energía ya no puede bajar por más que ~Var/Δ.

**6. Barras de error calibradas.** Ver la calibración de σ en la sección de medición (Santino).

> La tabla siguiente reúne todo lo medido. La última celda compara recién ahí con FCI.
