**Qué dio la reducción.**

- **Validación:** en las 18 distancias, $\hat H_6$ y $\hat H_5$ dan exactamente la energía FCI (error < 10⁻¹³ Ha).
- **Shots:** de 185 a 132 términos, pero casi los mismos grupos (37 → 35). La reducción **no ahorra shots, ahorra compuertas**.
- **Ansatz:** el compacto llega a FCI sin ruido con **16 CZ, 100 veces menos** que UCCSD ×2 en línea.

**Por qué funcionó:** las simetrías son información que ya conocemos. Sacándola, el espacio pasa de 256 estados a 32 y un circuito simple alcanza.

**Qué no funcionó o cuesta:** con menos de 7 capas aparecen mínimos locales; poner todas las CZ en cada capa traba el ansatz (fidelidad 0,96); y el ansatz en ladrillo **no escala** a moléculas grandes: es una solución para H₄, no un método general.
