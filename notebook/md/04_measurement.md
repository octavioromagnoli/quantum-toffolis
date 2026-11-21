## 3 · Medir la energía con shots  *(Santino)*

### Cómo medimos

- H tiene 184 términos de Pauli más una constante. Los agrupamos en conjuntos que se miden con **un mismo circuito**. Cada shot de un grupo da una muestra de O_g = Σ c_k P_k. Entonces E_g es la media muestral, σ_g² = var/N_g (con las covarianzas dentro del grupo incluidas), E = c₀ + Σ_g E_g y σ² = Σ_g σ_g², porque los grupos usan shots distintos.
- **Agrupamiento QWC (default)**: términos que conmutan qubit a qubit. Se miden rotando cada qubit (H para X, S†H para Y), **sin CZ extra**. Armamos los grupos por *sorted insertion* (Crawford et al., *Quantum* 5, 385 (2021)): por |c| decreciente, cada término va al primer grupo compatible. Así quedan **37 grupos**, contra los 46 del coloreo de `group_commuting`, y se necesitan **1,6–1,8× menos shots** para la misma σ (tabla de abajo).
- **Agrupamiento GC (opcional, `grouping="gc"`)**: términos que conmutan en general, **9 grupos**. Cada grupo se diagonaliza con un Clifford: tomamos generadores independientes (sobre GF(2)), construimos el estado estabilizador C|0⟩ que fijan y medimos en Z después de C†, porque C†PC = ±Z… para todo P del grupo. Necesita **2–2,8× menos shots que QWC**, pero agrega **4–23 CZ por grupo** con conectividad total (**4–77 con qubits en línea**).
- **Reparto de shots en dos etapas**: primero un piloto con el 10 % del presupuesto y N_g ∝ Σ|c_k|, para estimar σ_g. Después el resto con la asignación de Neyman, N_g ∝ σ_g, que es la que minimiza σ para un N total fijo. Las muestras de las dos etapas se juntan.
- **Semillas**: cada grupo y cada etapa tiene su propia semilla, derivada con `SeedSequence`. Al principio usábamos semillas consecutivas de Aer y vimos errores correlacionados entre grupos, con desvíos de 4σ.
- **Velocidad**: los conteos se arman con `np.bincount` sobre el `BitArray` del `SamplerV2`. Con 8 M de shots, `get_counts()` (un string por shot) se comía la mayor parte del tiempo. Además, los circuitos se transpilan una sola vez.

### Cuántos shots hacen falta

Para un estado |ψ⟩ y fracciones de shots f_g, σ² = Σ_g V_g/(f_g N) con V_g = ⟨O_g²⟩ − ⟨O_g⟩². Por eso **σ = √(K/N)**: la σ baja como 1/√N, y alcanza con calcular K una vez para cada estrategia. Las V_g se calculan con el estado FCI, al que el VQE llega con un error menor que 1 µHa. Esta cuenta **solo sirve para planificar**; la energía reportada sale siempre de los shots (`scripts/shots_study.py`).

Shots necesarios para **σ = 0,5 mHa** (reparto en dos etapas):

| d (Å) | QWC con coloreo (antes) | **QWC con sorted insertion** | GC con coloreo | **GC con sorted insertion** |
|---|---|---|---|---|
| 0,9 | 17,8 M | 11,1 M | 5,3 M | 4,3 M |
| 1,8 | 16,5 M | 9,1 M | 7,2 M | 4,4 M |
| 2,0 | 15,5 M | 8,6 M | 6,8 M | 4,1 M |
| 2,4 | 14,2 M | 8,0 M | 5,9 M | 3,6 M |
| 3,0 | 13,6 M | 8,0 M | 4,9 M | 2,9 M |

**Reparto de shots** (QWC con sorted insertion, 2,0 Å, K en mHa²·M shots):

| reparto | Σ\|c\| | Neyman con varianzas de HF | dos etapas (piloto + Neyman) | Neyman ideal (V_g exactas) |
|---|---|---|---|---|
| K | 3,69 | 24,6 | 2,15 | 2,13 |

Las dos etapas quedan a menos de 1 % del óptimo. Usar las varianzas del estado HF es **11× peor**: HF es un estado de la base computacional, así que da varianza cero a todos los grupos que miden solo Z y les asigna muy pocos shots.

**En 5 qubits** (reducción de Octavio, ansatz compacto): con el coloreo, la reducción costaba ~3× más shots que 8 qubits. Con sorted insertion ese costo desaparece: queda igual que en 8 qubits.

| d (Å) | QWC con coloreo | **QWC con sorted insertion** | GC con sorted insertion |
|---|---|---|---|
| 1,8 | 31 M | 11 M | 4,3 M |
| 2,4 | 28 M | 8,9 M | 3,5 M |
| 3,0 | 25 M | 7,9 M | 3,0 M |

(Shots para σ = 0,5 mHa, sin ruido; `scripts/mitigation_study.py`.)

**Presupuesto elegido**: con QWC, **8 M de shots por punto dan σ = 0,50–0,53 mHa** en la zona estirada (~0,6 mHa cerca del mínimo). Con GC, 4 M dan lo mismo.

### Calibración de σ: ¿la barra de error es honesta?

Repetimos la medición del **mismo circuito** con semillas independientes y miramos z = (E − E_exacta)/σ, donde E_exacta = ⟨ψ|H|ψ⟩ del circuito medido. Si σ está bien, z ~ N(0,1). Hicimos 500 repeticiones por configuración (300 para Σ|c|), a 2,0 y 3,0 Å, con un total de **17,8 mil millones de shots** (`scripts/calib_sigma.py`, ~35 min en 12 núcleos). También calculamos σ por **bootstrap**: re-muestrear los shots de cada grupo y recalcular E.

| configuración | d (Å) | σ media (mHa) | σ predicha √(K/N) | media de z | std de z [IC 95 %] | σ_boot/σ | \|z\|>2 |
|---|---|---|---|---|---|---|---|
| QWC, 8 M, dos etapas | 2,0 | 0,519 | 0,519 | −0,02 ± 0,05 | 1,01 [0,95; 1,07] | 0,997 | 5,0 % |
| QWC, 8 M, dos etapas | 3,0 | 0,500 | 0,500 | +0,04 ± 0,05 | 1,00 [0,94; 1,06] | 0,998 | 4,4 % |
| QWC, 1 M, dos etapas | 2,0 | 1,467 | 1,467 | +0,04 ± 0,05 | 1,03 [0,97; 1,09] | 0,999 | 4,8 % |
| QWC, 1 M, dos etapas | 3,0 | 1,413 | 1,413 | −0,02 ± 0,05 | 1,04 [0,98; 1,11] | 0,999 | 6,2 % |
| QWC, 8 M, Σ\|c\| | 2,0 | 0,679 | 0,679 | −0,16 ± 0,06 | 0,98 [0,91; 1,06] | 0,997 | 5,3 % |
| QWC, 8 M, Σ\|c\| | 3,0 | 0,681 | 0,681 | −0,02 ± 0,06 | 1,03 [0,95; 1,12] | 1,001 | 5,0 % |
| GC, 4 M, dos etapas | 2,0 | 0,507 | 0,507 | +0,01 ± 0,04 | 0,96 [0,91; 1,03] | 0,996 | 3,4 % |
| GC, 4 M, dos etapas | 3,0 | 0,424 | 0,424 | +0,01 ± 0,05 | 1,04 [0,98; 1,11] | 0,994 | 5,8 % |

- **σ está calibrada**: la dispersión de z es compatible con 1 en las 8 configuraciones. Con las 3600 repeticiones juntas: media de z = −0,006 ± 0,017, dispersión = 1,011 [0,988; 1,035], |z|>2 en 5,0 % (normal: 4,6 %), 11 casos con |z|>3 (esperados ~10) y prueba KS contra N(0,1) con p = 0,78. No hay correlación entre repeticiones seguidas (|r| < 0,12). Histogramas: `results/calib_sigma/hist_z.png`. La fórmula √(K/N) predice la σ media hasta el tercer dígito, y de 1 M a 8 M la σ baja en √8 = 2,83 (medido: 2,83).
- **El 1,16 de antes era ruido estadístico**: con 60 repeticiones, aunque σ sea exacta, la dispersión de z cae en [0,82; 1,18] el 95 % de las veces.
- **El bootstrap coincide** con la σ analítica dentro de 0,7 %: la σ no depende de suponer normalidad.
- La media de z = −0,16 ± 0,06 (Σ|c| a 2,0 Å) es un desvío de 2,7σ entre unas 18 pruebas, algo que por azar pasa ~10 % de las veces. No es un sesgo: el mismo circuito medido con dos etapas da −0,02 ± 0,05.

### Qué no funcionó o descartamos

- **Neyman con varianzas de HF**: 11× peor que el piloto (ver arriba).
- **Reusar shots entre grupos** (cada base QWC también estima los Paulis compatibles de otros grupos, promediados por shots): sobre el coloreo gana 1,4–2,3×, pero sobre *sorted insertion* casi no suma. Solo gana 1,5× cerca del mínimo, y en la zona estirada empeora un poco porque el promedio por shots ignora las covarianzas. No lo implementamos.
- **Calibrar con 60 repeticiones**: no alcanza para distinguir 1,16 de 1.

### Limitaciones

- La calibración es **sin ruido**: verifica la estadística de los shots, no el sesgo del ruido. Con ruido, σ sigue midiendo solo la dispersión estadística, y el sesgo se trata en el notebook de extras (mitigación).
- GC baja los shots pero **suma CZ** al circuito medido: con 0,1 % de error por CZ, ~50 CZ más (qubits en línea) son ~5 % más de error. Por eso proponemos QWC para las corridas con ruido: además permite post-seleccionar por N_α y N_β en el grupo que mide solo Z.
- K se planificó con el estado FCI. Si el ansatz cambia, el piloto vuelve a estimar las σ_g en cada punto, así que el reparto se ajusta solo.
