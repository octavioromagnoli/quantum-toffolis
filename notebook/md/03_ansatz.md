## 2 · El ansatz  *(Gonzalo)*

Partimos de UCCSD con Trotter: 8 simples + 18 dobles = 26 excitaciones por paso. Cada excitación es exp(−iθG) con G = i(T − T†). Sus términos de Pauli conmutan, así que el circuito es exactamente Π exp(−iθ_k G_k). Eso nos deja optimizar con un simulador numpy (gradiente analítico) y medir **el mismo circuito**.

El problema: UCCSD×2 llega a FCI, pero cuesta ~584 CZ (~1600 en línea). Con 0,1 % de error por CZ eso no sirve.

### ADAPT-VQE: solo las excitaciones que hacen falta

ADAPT (Grimsley et al., 2019) arranca de |HF⟩ y agrega de a una la excitación del pool con gradiente máximo, g_k = ⟨ψ|i[G_k, H]|ψ⟩. Después de cada agregado reoptimiza todos los parámetros. Pool: las 26 excitaciones de UCCSD, que se pueden repetir (`h4vqe/adapt.py`).

**Una lista común para toda la curva.** En cada paso elegimos el operador con mayor gradiente en *cualquiera* de las 18 distancias y reoptimizamos cada una por separado. Así el circuito (y la cantidad de CZ) es el mismo en todos los puntos, y solo cambian los θ.

**Dónde cortar sin usar FCI.** Usamos la varianza de la energía ⟨H²⟩ − ⟨H⟩², que vale 0 solo en un autoestado. Para elegir el corte la calculamos con el estado del simulador: al pasar de 12 a 13 operadores, la varianza máxima sobre la curva cae de 3,2·10⁻³ a 2,6·10⁻⁴ Ha² (y el error máximo, de 3,7 a 0,6 mHa). La misma varianza **medida con shots** por eco de Loschmidt (sección 4) lo confirma: 2,8·10⁻⁴ a 1,8 Å, 7,3·10⁻⁵ a 2,4 Å y 5·10⁻⁶ a 3,0 Å, contra 2,6·10⁻⁴, 6,5·10⁻⁵ y 6,7·10⁻⁶ calculadas.

**Completar los pares de espín.** ADAPT eligió las simples α 0→2 y β 5→7, pero no sus pares β 4→6 y α 1→3, así que el estado rompe un poco el espín (S² ≤ 0,013). Agregamos esas dos simples **al lado de su pareja** (`ADAPT_FERMIONIC`, 15 operadores): el error máximo baja de 0,59 a **0,38 mHa** en toda la curva y S² baja a ≤ 6,5·10⁻³, por 9 CZ más.

**Dos arranques en la optimización.** `optimize_params` prueba siempre el warm start (θ* de la distancia anterior) **y** θ = 0 (|HF⟩), y se queda con la energía más baja (es variacional). Solo con warm start, el optimizador caía en mínimos malos en algunos puntos: 2 mHa a 2,8 Å con la lista de 15, 0,84 mHa a 1,4 Å con la de 13 y 9,8 mHa a 2,0 Å con QEB.

| ansatz | parámetros | CZ (todos con todos) | CZ (línea) | error máx. vs FCI (sin ruido) | S² máx. |
|---|---|---|---|---|---|
| UCCSD ×1 | 26 | 292 | 806 | 4,5 mHa (2,0 Å) | hasta 0,85 |
| UCCSD ×2 | 52 | ~584 | ~1600 | < 0,001 mHa | ≈ 0 |
| ADAPT, 13 ops (`ADAPT_FERMIONIC_13`) | 13 | 185 | 526 (transpiler) / 370 (ruteo fSWAP) | 0,59 mHa (2,0 Å) | 0,013 |
| **ADAPT + pares de espín, 15 ops (`ADAPT_FERMIONIC`)** | **15** | **194** | **530 (transpiler) / 378 (ruteo fSWAP)** | **0,38 mHa** | **6,5·10⁻³** |
| ADAPT-QEB, 15 ops (extra) | 15 | 159 | ~451 (transpiler) | 0,81 mHa | 0,074 |

Las CZ en línea con ruteo fSWAP salen de `routing.line_cz_count(routed_ansatz(m, lista, LINE_ORDER)[0])` (ver abajo). `cz_count(circuito, "line")` usa el transpiler común y da el número más alto.

El pool QEB usa excitaciones de *qubit* (sin cadenas Z de Jordan-Wigner). Conservan N_α y N_β, tocan solo 2 o 4 qubits y cuestan menos CZ.

### Qubits en línea: ruteo propio con fSWAP

Con SWAP comunes, el transpiler multiplica las CZ ~2,8× en línea. Nuestro ruteo (`h4vqe/routing.py`), inspirado en Co-ADAPT (Ramôa et al., 2026) y en las redes fSWAP (Kivlichan et al., 2018), hace tres cosas:

1. Antes de cada excitación junta sus modos en qubits **contiguos** con **fSWAP** (2 CZ cada uno; intercambian dos modos fermiónicos *con* el signo de Jordan-Wigner).
2. **No deshace** la permutación: el orden de los modos se sigue de forma clásica.
3. Aplica cada doble con un circuito pensado para la línea (18 CZ, casi todas entre vecinos) y mide la energía con Ĥ escrito en el **orden final** de modos (`hamiltonian_in_order`).

Los θ son los mismos (las fSWAP solo reetiquetan modos). Verificamos que la energía coincide hasta 1e-15 (`tests/test_routing.py`) y que la medición con shots no tiene sesgo (16 semillas: media de z = 0,14, desvío 0,91). El orden inicial se eligió entre 72 candidatos: [α3 β3 α1 β1 α2 β2 α0 β0]. Resultado: **530 → 378 CZ en línea** (−29 %).

Costo: Ĥ en el orden final se agrupa mucho peor (69 grupos contra 37 con el agrupamiento de la sección 3), así que necesita más shots para la misma σ. Por eso la curva central se mide con el circuito sin rutear y Ĥ original (es el mismo estado); el circuito ruteado se usa para las CZ en línea y para la simulación con ruido.

### Qué no funcionó

- **Umbral de gradiente 1e-4:** ADAPT se frena en mesetas a 1,6, 2,8 y 3,0 Å con 0,6-2 mHa de error, y después el gradiente vuelve a crecer. Con 1e-5 llega a FCI en toda la curva.
- **Agregar las simples de espín al final de la lista:** con el warm start cae en un mínimo malo (2,5 mHa, S² = 0,52). Tienen que ir al lado de su pareja.
- **Optimizar solo con warm start entre distancias:** cae en mínimos locales malos (hasta 9,8 mHa con QEB a 2,0 Å, S² = 0,78). Se arregló probando también θ = 0. QEB es más barato, pero rompe más el espín (S² ≤ 0,074), así que queda como extra.
- **Reordenar los operadores** para necesitar menos fSWAP (314 CZ en línea): el error sube a 5,2 mHa. El orden de los operadores cambia el ansatz.
- **Co-ADAPT** (elegir cada operador por gradiente / costo en línea): necesita 16-19 operadores para precisión química y termina en ~460-470 CZ en línea. En H4 casi todas las dobles hacen falta, así que penalizar por costo solo agrega operadores.
- **Circuitos explícitos de excitación sin ruteo propio** (escalera de CNOT + RY multicontrolada): exactos, pero el transpiler los rutea peor que `PauliEvolutionGate` (705 contra 526 CZ en línea).
- **Probar layouts y semillas del transpiler** para conectividad en línea: el mejor resultado fue ~515 CZ, contra ~526 por defecto.

### Extra: comparación de ansätze

**Sin ruido, en las 18 distancias** (`scripts/compare_ansatze.py`). Los θ se optimizan como en la curva. FCI se usa solo para reportar el error. S² y las CZ en línea varían un poco entre corridas (el optimizador puede llegar a θ distintos con la misma energía, y el transpiler usa semillas). CZ y grupos a 2,4 Å (el circuito es el mismo en toda la curva).

| ansatz | qubits | θ | CZ todos con todos | CZ en línea | grupos | error máx. | error máx. (d ≥ 1,8 Å) | S² máx. |
|---|---|---|---|---|---|---|---|---|
| UCCSD ×1 | 8 | 26 | 292 | 804 | 37 | 4,46 mHa | 4,46 mHa | 0,85 |
| UCCSD ×2 | 8 | 52 | 582 | 1648 | 37 | 0,001 mHa | < 0,001 mHa | ≈ 0 |
| ADAPT 13 | 8 | 13 | 185 | 370 (fSWAP) | 37 | 0,59 mHa | 0,59 mHa | 0,013 |
| **ADAPT 15 (curva central)** | 8 | 15 | 194 | **378 (fSWAP)** | 37 | **0,38 mHa** | 0,38 mHa | ≤ 7·10⁻³ |
| ADAPT 15 reducido | 6 | 15 | 194 | ~510 | 37 | 0,38 mHa | 0,38 mHa | ≤ 7·10⁻³ |
| ADAPT 15 reducido | 5 | 15 | 172 | 427 | 35 | 0,38 mHa | 0,38 mHa | ≤ 7·10⁻³ |
| ADAPT-QEB 15 | 8 | 15 | 159 | 453 | 37 | 0,81 mHa | 0,81 mHa | 0,07-0,14 |
| compacto (Octavio), L = 5 | 5 | 30 | 10 | 10 | 35 | 56 mHa | 50 mHa | 0,025 |
| compacto, L = 6 | 5 | 35 | 12 | 12 | 35 | 3,8 mHa | 3,8 mHa | 3·10⁻³ |
| **compacto, L = 7** | 5 | 40 | **14** | **14** | 35 | **< 0,001 mHa** | < 0,001 mHa | ≈ 0 |
| compacto, L = 8 (default) | 5 | 45 | 16 | 16 | 35 | < 0,001 mHa | < 0,001 mHa | ≈ 0 |

Qué muestra la tabla:

- **Ansatz químico:** ADAPT 15 llega a 0,38 mHa con 3,5 veces menos parámetros y 3 veces menos CZ que UCCSD ×2. QEB ahorra CZ con conectividad total, pero en línea queda peor que el ruteo fSWAP y rompe más el espín.
- **Reducir qubits no le sirve a ADAPT:** la energía es idéntica en 8, 6 y 5 qubits (la reducción está bien hecha), pero las CZ casi no bajan. En línea quedan peor que el ruteo fSWAP de 8 qubits, porque la reducción mezcla los modos y las excitaciones dejan de ser locales.
- **La reducción rinde con el compacto:** en 5 qubits, una brickwork de RY + CZ entre vecinos (sin SWAP) llega a FCI con **14 CZ** (L = 7). Con 12 CZ ya falla (3,8 mHa). No conserva N ni S² por construcción (solo por la penalización de S² en la optimización): es "hardware-efficient", no químico.

**Con ruido, a 2,4 Å** (`scripts/noise_ansatze.py`; modelo de `noise.py`, qubits en línea; 500 k shots sin mitigar, 2 M con ZNE de plegado global con λ = 1, 3). CCSD − FCI = −30,6 mHa.

| ansatz | CZ (línea) | compuertas de 1 qubit | error sin ruido | con ruido, sin mitigar | con ZNE |
|---|---|---|---|---|---|
| compacto 5q L = 8 | 16 | 141 | 0,0 mHa | +12,9 ± 2,2 mHa | **+0,2 ± 2,2 mHa** |
| ADAPT 15 + fSWAP | 313 | 1377 | 0,2 mHa | +178 ± 3 mHa | +70 ± 3 mHa |
| UCCSD ×2 | 1546 | 3801 | 0,0 mHa | +548 ± 2 mHa | — |

**Con ruido, la cantidad de compuertas decide todo.** Los tres ansätze son igual de buenos sin ruido, pero con 0,1 % por CZ solo el de 16 CZ le gana a CCSD; ADAPT y UCCSD ×2 quedan a cientos de mHa. La ZNE no alcanza a salvar 313 CZ.

Con muchos shots, el compacto con mitigación llega a **precisión química con ruido, contando la σ** (sección de mitigación de Santino, `scripts/mitigation_shots.py`: ZNE lineal + post-selección por simetría + corrección de lectura, 54 M de shots, con 2 % de error de lectura): −0,48 ± 0,54 mHa a 1,8 Å, −0,09 ± 0,50 a 2,4 Å y +0,31 ± 0,48 a 3,0 Å.

**Cuál usar:** sin ruido (la curva central), ADAPT 15. Es químico (conserva N exacto y casi S²), tiene pocos parámetros con significado físico, y se puede justificar sin FCI con la varianza. Con ruido (extras), el compacto de 5 qubits.

### Limitaciones y próximos pasos

- S² ≤ 6,5·10⁻³ (no es exactamente 0), porque el pool no está adaptado al espín: los pares de dobles αβ tampoco comparten parámetro.
- Para bajar más las CZ: operadores CEO de 9 CNOT, podar operadores y orbitales localizados u optimizados (rotarlos clásicamente no cuesta CZ).
