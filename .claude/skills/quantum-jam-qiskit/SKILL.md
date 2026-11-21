---
name: quantum-jam-qiskit
description: Usar para cualquier tarea de computación cuántica con Qiskit: circuitos, algoritmos (Grover, DJ, BV, Simon, QPE, Shor), BB84, teleportación, VQE, SQD y química con PySCF.
---

# Quantum Jam 2026: computación cuántica con Qiskit

Esta skill resume el curso de la Quantum Jam 2026 / Qiskit Fall Fest 2026 del ITBA (repo `Computer-Society-ITBA/quantumjam-workshops`): 4 clases, 2 partes extra y 2 notebooks de química cuántica. Sirve para resolver desafíos de hackathon con Qiskit 2.x.

El usuario no programa a mano: dirige y valida. Por eso, antes de escribir código, **explicá el plan**: qué algoritmo, por qué, cuántos qubits, qué resultado se espera y cómo lo vas a validar. Después implementá en scripts o notebooks chicos y reportá siempre un chequeo de cordura.

---

## 0. Reglas de trabajo (siempre)

### Stack e imports
```python
import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.circuit import Parameter, ParameterVector
from qiskit.circuit.library import PauliEvolutionGate, efficient_su2, QFT
from qiskit.quantum_info import Statevector, Operator, Pauli, SparsePauliOp
from qiskit.visualization import plot_histogram, plot_bloch_multivector
from qiskit.transpiler import generate_preset_pass_manager
from qiskit_aer import AerSimulator
from qiskit_ibm_runtime import SamplerV2 as Sampler, EstimatorV2 as Estimator, QiskitRuntimeService
from qiskit_ibm_runtime.fake_provider import FakeVigoV2, FakeFez
from qiskit.primitives import StatevectorEstimator          # exacto, sin ruido
from qiskit_aer.primitives import EstimatorV2 as AerEstimator, SamplerV2 as AerSampler
```
Instalación: `pip install qiskit[visualization] qiskit-ibm-runtime qiskit-aer qiskit_qasm3_import` (+ `pyscf` para química). En Colab hay que reinstalar en cada sesión nueva.

### Trampas que hay que evitar
1. **Little-endian.** En Qiskit, q0 es el carácter de más a la DERECHA en bitstrings, en `Statevector.from_label` y en strings de Pauli (`'ZYI'` = Z en q2, Y en q1, I en q0). Los libros usan el orden inverso, así que las matrices de más de 2×2 no coinciden con las de los libros. Si un resultado sale espejado, es esto. Comentalo en el código cada vez que importe.
2. **`execute` no existe** en Qiskit 1.x/2.x. Usá las primitivas V2 o `backend.run(isa_circuit)`.
3. **Transpilá siempre** antes de ejecutar en un backend: `pm = generate_preset_pass_manager(backend=backend, optimization_level=1); isa = pm.run(qc)`.
4. **Estimator:** el circuito va SIN mediciones y el observable con `obs.apply_layout(isa.layout)` después de transpilar.
5. **Sampler:** el circuito va CON mediciones. Con `measure_all()` los counts están en `result[0].data.meas.get_counts()`; con `QuantumCircuit(n, n)` + `measure`, en `result[0].data.c.get_counts()`.
6. `Statevector.from_instruction(qc)` no acepta circuitos con mediciones.
7. **Validá en escalera:** primero exacto (`Statevector` / `StatevectorEstimator`), después con shots (`AerSimulator`), después con ruido (`FakeVigoV2`, `FakeFez`, `AerEstimator` con `default_precision`) y recién al final hardware real.

### Chequeos de cordura que hay que reportar
- Las probabilidades suman 1 y el estado es válido (`sv.is_valid()`).
- VQE/SQD: la energía tiene que estar **por encima** de la exacta (principio variacional). Si da por debajo, hay un bug.
- Grover: P(solución) alta, con r = ⌊π/4·√(N/M)⌋ iteraciones.
- DJ/BV: resultado determinista (100% en un solo bitstring en simulación ideal).
- Profundidad (`qc.depth()`) y CNOTs (`transpile(qc, basis_gates=['cx','rz','sx','x']).count_ops()['cx']`): las CNOTs son el principal costo en hardware.

### Qué algoritmo para qué problema
| Problema | Herramienta |
|---|---|
| ¿f es constante o balanceada? | Deutsch / Deutsch-Jozsa |
| f(x) = x·s, encontrar s | Bernstein-Vazirani |
| f dos a uno con f(x) = f(x⊕b), encontrar b | Simon |
| Búsqueda no estructurada / satisfacer una condición | Grover |
| Autovalor de un unitario | QPE |
| Factorizar, período de a^x mod N | Shor (QPE + fracciones continuas) |
| Mandar un qubit con canal clásico + EPR | Teleportación |
| 2 bits clásicos con 1 qubit + EPR | Codificación superdensa |
| Distribuir una clave secreta y detectar espías | BB84 |
| Energía del estado fundamental (química, Hamiltonianos) | VQE o SQD |

Patrón común de los algoritmos: **superposición (H) → marcar con fase (oráculo o kickback) → interferencia (H, QFT o difusor) → medir.**

---

## 1. Fundamentos (Clase 1)

**Qubit:** |ψ⟩ = α|0⟩ + β|1⟩, con |α|² + |β|² = 1. Al medir se obtiene 0 con probabilidad |α|² y 1 con |β|² (regla de Born), y el estado colapsa.

**Notación de Dirac:** ket = columna; bra = ⟨v| = |v⟩† (transpuesto y conjugado; las constantes salen conjugadas); braket ⟨w|v⟩ = Σ wᵢ* vᵢ, con ⟨w|v⟩ = ⟨v|w⟩*; ket-bra |v⟩⟨w| = matriz (no conmuta). Norma: √⟨v|v⟩. Probabilidad de obtener |φ⟩ al medir |ψ⟩: |⟨φ|ψ⟩|².

**Bases de un qubit:**
- Z: |0⟩, |1⟩ (etiquetas `'0'`, `'1'`).
- X (Hadamard): |±⟩ = (|0⟩ ± |1⟩)/√2 (etiquetas `'+'`, `'-'`).
- Y (circular): |±i⟩ = (|0⟩ ± i|1⟩)/√2 (etiquetas `'r'`, `'l'`).

**Muchos qubits:** 2ⁿ amplitudes. Superposición uniforme (1/√2ⁿ)Σ|i⟩. GHZ = (|0…0⟩ + |1…1⟩)/√2. W = (1/√N)(|10…0⟩ + … + |0…01⟩), con un único 1.

**Operadores:** Â = Σ aᵢⱼ|uᵢ⟩⟨uⱼ|. Adjunto A† = (A*)ᵀ; (AB)† = B†A†. Traza tr(A) = Σ⟨k|A|k⟩.

Paulis: X = [[0,1],[1,0]], Y = [[0,−i],[i,0]], Z = [[1,0],[0,−1]]; H = (1/√2)[[1,1],[1,−1]].
- X: |0⟩↔|1⟩; |−⟩ → −|−⟩.
- Z: |1⟩ → −|1⟩; |+⟩↔|−⟩.
- Y = iXZ: |0⟩ → i|1⟩, |1⟩ → −i|0⟩.
- H: |0⟩ → |+⟩, |1⟩ → |−⟩ y viceversa (intercambia las bases Z y X).

**Medición:**
- Observable = operador hermítico (A† = A): autovalores reales, autovectores ortogonales. Los resultados posibles son los autovalores.
- Proyector P = |φ⟩⟨φ|, con P² = P y P† = P. Clausura: Σ Pᵢ = I.
- P(λᵢ) = |⟨φᵢ|ψ⟩|²; si el autovalor es degenerado, P(λ) = ⟨ψ|P_m|ψ⟩.
- Estado después de medir: P_m|ψ⟩/√P(λ).
- Descomposición espectral: A = Σ λᵢ|φᵢ⟩⟨φᵢ|.
- Valor esperado: ⟨A⟩ = ⟨ψ|A|ψ⟩ = Σ λᵢ P(λᵢ); por ejemplo, ⟨Z⟩ = P(0) − P(1).

**Evolución:** los unitarios cumplen UU† = I: conservan la norma, son reversibles y el producto de unitarios es unitario. Schrödinger: iħ d|ψ⟩/dt = Ĥ|ψ⟩, cuya solución es U = e^{−iĤt/ħ}.

**Fases y Bloch:** la fase global e^{iθ}|ψ⟩ es irrelevante; la fase relativa entre α y β sí importa. Parametrización: |ψ⟩ = cos(θ/2)|0⟩ + e^{iφ} sin(θ/2)|1⟩, con coordenadas (sinθcosφ, sinθsinφ, cosθ). Los Paulis son giros de 180° en su eje. Dos estados ortogonales se distinguen con certeza midiendo en la base adecuada. Para medir en X: H y medir. Para medir en Y: S†, H y medir.

```python
psi = Statevector([(1+1j)/2, (1-1j)/2])
psi.probabilities_dict(); psi.inner(Statevector.from_label('+'))   # ⟨psi|+⟩
np.conj(psi.data)                                  # el bra
Operator.from_label('H'); Statevector.from_label('0').evolve(Operator.from_label('H'))
psi.expectation_value(SparsePauliOp('Z'))

def ghz_sv(n):
    v = np.zeros(2**n); v[0] = v[-1] = 1/np.sqrt(2); return Statevector(v)

def w_sv(n):
    v = np.zeros(2**n)
    for k in range(n): v[1 << k] = 1/np.sqrt(n)
    return Statevector(v)
```

---

## 2. Circuitos y Qiskit (Clase 2)

**Compuertas de un qubit:** son unitarias. La columna j de la matriz es la imagen de |j⟩.

| Compuerta | Efecto | Qiskit |
|---|---|---|
| X, Y, Z | Giros de 180° en x, y, z | `qc.x(q)` etc. |
| R(θ) | diag(1, e^{iθ}) | `qc.p(θ,q)` |
| S | R(π/2) | `qc.s(q)` |
| T | R(π/4) | `qc.t(q)` |
| H | Interferencia: (α+β)/√2 · \|0⟩ + (α−β)/√2 · \|1⟩ | `qc.h(q)` |
| RX, RY, RZ | Rotaciones de ángulo θ | `qc.ry(θ,q)` |

Dato útil: RY(θ)|0⟩ = cos(θ/2)|0⟩ + sin(θ/2)|1⟩.

**Reglas de los circuitos:**
- Se leen de izquierda a derecha y todos los qubits arrancan en |0⟩.
- Serie: U₃U₂U₁|ψ⟩. Paralelo: ⊗.
- Prohibido: lazos, FANIN y FANOUT. No clonación: no se puede copiar un estado desconocido.

**Varios qubits:**
- |v⟩⊗|w⟩, con ⟨v₁w₁|v₂w₂⟩ = ⟨v₁|v₂⟩⟨w₁|w₂⟩.
- Un estado producto es separable; un estado entrelazado no se puede factorizar.
- Bell: Φ± = (|00⟩ ± |11⟩)/√2 y Ψ± = (|01⟩ ± |10⟩)/√2. Se obtienen con H(q0) + CX(0,1), más X o Z antes.
- CNOT: |x,y⟩ → |x, x⊕y⟩ (`qc.cx(c,t)`).
- Toffoli: |a,b,c⟩ → |a,b,c⊕ab⟩ (`qc.ccx`). También `qc.swap`, `qc.cz` y `qc.mcx(ctrls, t)`.
- Transformada de Hadamard: H^{⊗n}|0ⁿ⟩ = (1/√2ⁿ)Σₓ|x⟩ y H^{⊗n}|x⟩ = (1/√2ⁿ)Σ_y (−1)^{x·y}|y⟩.

**Qiskit Patterns:**
1. **Map:** problema → circuitos y observables.
2. **Optimize:** transpilar a un circuito ISA (compuertas nativas + conectividad). `optimization_level` va de 0 a 3: más nivel da menos profundidad a cambio de más tiempo de compilación.
3. **Execute:** Sampler (distribución de bitstrings) o Estimator (valores esperados).
4. **Post-process:** análisis clásico.

```python
backend = AerSimulator()
pm = generate_preset_pass_manager(backend=backend, optimization_level=1)

# Sampler
bell = QuantumCircuit(2); bell.h(0); bell.cx(0, 1); bell.measure_all()
counts = Sampler(mode=backend).run([pm.run(bell)], shots=1000).result()[0].data.meas.get_counts()
plot_histogram(counts, sort='value_desc')

# Estimator
qc = QuantumCircuit(2); qc.h(0); qc.cx(0, 1)
isa = pm.run(qc)
obs = [SparsePauliOp(l).apply_layout(isa.layout) for l in ["ZZ", "XX", "ZI"]]
evs = Estimator(mode=backend).run([(isa, obs)]).result()[0].data.evs      # [1, 1, ~0]

# Circuito dinámico
dyn = QuantumCircuit(2, 1); dyn.h(0); dyn.measure(0, 0)
with dyn.if_test((dyn.clbits[0], 1)):
    dyn.x(1)

# Parámetros
theta = Parameter('θ'); pc = QuantumCircuit(1); pc.rx(theta, 0)
bound = pc.assign_parameters({theta: np.pi/2})

# Transpilación y profundidad
isa_ghz = generate_preset_pass_manager(optimization_level=3, backend=FakeVigoV2()).run(ghz_circuit)
print(ghz_circuit.depth(), isa_ghz.depth())

# OpenQASM 3
import qiskit.qasm3
s = qiskit.qasm3.dumps(qc); qc2 = qiskit.qasm3.loads(s)
```

**Conceptos de runtime y ruido:**
- Modos de ejecución: *job* (trabajo suelto), *session* (iterativo con QPU reservado: el indicado para VQE) y *batch* (trabajos independientes agrupados).
- Supresión y mitigación: Dynamical Decoupling (pulsos en tiempos muertos contra la decoherencia), readout/TREX (errores de medición) y ZNE (corre con varios niveles de ruido y extrapola a cero).
- OpenQASM 3 agrega flujo de control clásico en tiempo real (if/while/for) respecto de OpenQASM 2.
- Sin Python: REST API de Qiskit Runtime, con API key de IBM Cloud (→ token IAM) y CRN de la instancia.

Hardware real:
```python
QiskitRuntimeService.save_account(channel="ibm_quantum_platform", token=API_KEY, instance=CRN, name="fallfest-2026")
service = QiskitRuntimeService(name="fallfest-2026")
backend = service.least_busy(operational=True, simulator=False)
```

---

## 3. Comunicaciones y oráculos (Clase 3)

Actores: Alice (emisora), Bob (receptor), canal clásico o cuántico, ruido y Eve (espía).

### Teleportación
Alice y Bob comparten Φ+. Pasos: CNOT(ψ → EPR de Alice), H(ψ), Alice mide M1 y M2 y los manda; Bob aplica Z^{M1} X^{M2}. No clona (el original se destruye) ni comunica más rápido que la luz (hacen falta los 2 bits).
```python
qc = QuantumCircuit(3, 3)
qc.ry(1.2, 0)                         # estado a enviar
qc.h(1); qc.cx(1, 2)                  # EPR: q1 de Alice, q2 de Bob
qc.cx(0, 1); qc.h(0)
qc.measure(0, 0); qc.measure(1, 1)
with qc.if_test((qc.clbits[1], 1)): qc.x(2)
with qc.if_test((qc.clbits[0], 1)): qc.z(2)
qc.measure(2, 2)                      # P(1) ≈ sin²(0.6)
```

### Codificación superdensa
Alice aplica X^{M2} Z^{M1} sobre su mitad del EPR (00 → I, 01 → X, 10 → Z, 11 → ZX) y manda el qubit. Bob aplica CNOT(Alice → Bob) y H(Alice) y mide. Qiskit imprime 'q1q0' = 'M2M1'.

### Oráculo
U_f|x,y⟩ = |x, y⊕f(x)⟩. Con el auxiliar en |−⟩: U_f|x⟩|−⟩ = (−1)^{f(x)}|x⟩|−⟩ (phase kickback + paralelismo cuántico).

- **Deutsch** (1 bit) y **Deutsch-Jozsa** (n bits): |0ⁿ⟩|1⟩ → H en todos → U_f → H en la entrada → medir. Todo 0 → constante; otra cosa → balanceada. Clásicamente hacen falta 2^{n−1}+1 consultas; cuánticamente, 1.
  - Oráculo de paridad: CNOT de cada entrada al auxiliar (balanceada).
  - f = 1: X en el auxiliar (constante).
- **Bernstein-Vazirani:** mismo circuito con f(x) = x·s; el oráculo son CNOTs desde los qubits donde sᵢ = 1. Medís s en 1 consulta (contra n clásicas).
- **Simon:** U_f|x⟩|0ⁿ⟩ = |x⟩|f(x)⟩ con f(x) = f(x⊕b). Aplicás H, U_f, H y medís la entrada: cada y cumple y·b = 0 mod 2. Con ~n corridas, eliminación gaussiana mod 2 da b. Es el primer algoritmo con ventaja exponencial.

```python
def bernstein_vazirani(s):
    n = len(s); q = QuantumCircuit(n + 1, n)
    q.x(n); q.h(range(n + 1))
    for i, bit in enumerate(reversed(s)):     # little-endian
        if bit == '1': q.cx(i, n)
    q.h(range(n)); q.measure(range(n), range(n))
    return q
```

### Extra: el camino a Shor
- **Phase kickback:** si U|v⟩ = e^{iθ}|v⟩, U controlada con el control en |+⟩ da (|0⟩ + e^{iθ}|1⟩)|v⟩/√2.
- **QFT:** |j⟩ → (1/√2ⁿ) Σₖ e^{2πijk/2ⁿ}|k⟩. Se arma con H + rotaciones controladas R_k = diag(1, e^{2πi/2^k}) + SWAPs finales; usa O(n²) compuertas. En Qiskit: `QFT(n)` (inversa: `QFT(n, inverse=True)`).
- **QPE:** m qubits de conteo en |+⟩, U^{2^k} controlada desde el qubit k, QFT inversa y medir. Da la fase j de U|u⟩ = e^{2πij}|u⟩, con precisión que crece con m.
- **Shor:**
  1. Elegir a con mcd(a, N) = 1 (si no, ya hay un factor).
  2. U|x⟩ = |ax mod N⟩ tiene autovalores e^{2πis/r}; QPE sobre |1⟩ (que es la suma de todos los |u_s⟩) da s/r.
  3. Fracciones continuas → r. Si r es par y a^{r/2} ≢ −1 (mod N), los factores son mcd(a^{r/2} ± 1, N). Si no, probar otro a.

  Ejemplo: N = 15, a = 7 → r = 4 → mcd(48, 15) = 3 y mcd(50, 15) = 5. Romper RSA requiere millones de qubits lógicos.

---

## 4. BB84 y Grover (Clase 4)

### Criptografía y BB84
Simétrica: misma clave, el problema es distribuirla. Asimétrica: RSA, vulnerable a Shor. BB84 usa un canal clásico público autenticado y un canal cuántico (fotones: |0⟩ = →, |1⟩ = ↑, |+⟩ = ↗, |−⟩ = ↖).

1. Alice genera bits al azar.
2. Codifica cada uno en una base al azar (Z: 0 → |0⟩, 1 → |1⟩; X: 0 → |+⟩, 1 → |−⟩).
3. Bob mide en una base al azar.
4. Sifting: comparan bases (no bits) por el canal público.
5. Se quedan con las coincidencias: ~n/2 bits.

Después sacrifican m bits para estimar el QBER; si supera el umbral, hubo espía (siguen reconciliación y amplificación de privacidad).

Ataque interceptar y reenviar: Eve acierta el bit con P = 0,75, lo que genera ~25% de QBER. P(detectar) = 1 − (3/4)^m; m ≈ 72 da 0,999999999. Otros ataques: MitM y photon number splitting.

Simulación en Qiskit: por cada bit, `QuantumCircuit(1,1)`; `x` si bit = 1, `h` si la base es X; Bob aplica `h` si mide en X y mide. Eve mide y re-prepara en el medio.

### Grover (extra)
Búsqueda no estructurada en O(√N) (contra O(N) clásico). Es óptima (BBBV).
- |s⟩ = H^{⊗n}|0ⁿ⟩ = cos(θ/2)|w⊥⟩ + sin(θ/2)|w⟩, con sin(θ/2) = 1/√N.
- Oráculo U_w: invierte el signo de |w⟩.
- Difusor U_s = 2|s⟩⟨s| − I (inversión sobre la media) = H X (Z multicontrolada) X H, salvo fase global.
- Cada G = U_s U_w rota θ. Iteraciones: r ≈ π/4 · √(N/M). Iterar de más EMPEORA el resultado.

```python
def oraculo(w):                       # w en orden Qiskit, ej. '101'
    n = len(w); o = QuantumCircuit(n)
    ceros = [i for i, b in enumerate(reversed(w)) if b == '0']
    if ceros: o.x(ceros)
    o.h(n-1); o.mcx(list(range(n-1)), n-1); o.h(n-1)
    if ceros: o.x(ceros)
    return o

def difusor(n):
    d = QuantumCircuit(n); d.h(range(n)); d.x(range(n))
    d.h(n-1); d.mcx(list(range(n-1)), n-1); d.h(n-1)
    d.x(range(n)); d.h(range(n)); return d

n, w = 3, '101'; r = int(np.floor(np.pi/4*np.sqrt(2**n)))
qc = QuantumCircuit(n); qc.h(range(n))
for _ in range(r): qc.compose(oraculo(w), inplace=True); qc.compose(difusor(n), inplace=True)
qc.measure_all()                      # ~95% '101'
```
Con 2 qubits y w = 11: H H · CZ · H H · X X · CZ · X X · H H da 100% en una iteración.

---

## 5. Química cuántica (material extra, partes 1 y 2)

Objetivo: el **autovalor más chico de Ĥ** (energía del estado fundamental).

**Conceptos:**
- Spin-orbital (α/β) = un bit, por el principio de Pauli.
- Born-Oppenheimer: los núcleos quedan fijos y se calcula la energía a cada geometría. H₂ tiene su mínimo en ~0,74 Å.
- Unidades: Å y Hartree (≈27,2 eV). **Precisión química: 1,6 mHa.**
- Base de orbitales: STO-3G < 6-31G < cc-pVDZ. Cada orbital cuesta 2 qubits.
- Hartree-Fock: campo medio, una sola configuración. FCI: exacto, costo exponencial. E_corr = E_FCI − E_HF (H₂ STO-3G a 0,735 Å: HF −1,11700, FCI −1,13731, correlación −20,3 mHa; a 2,5 Å: −233 mHa).
- Correlación fuerte: aparece al estirar o romper enlaces; ahí HF falla.

**Configuraciones como bitstrings** (convención IBM): α a la derecha, β a la izquierda, y en cada mitad los orbitales de menor energía a la derecha. El HF de H₂ es `0101`. El estado exacto es una superposición de configuraciones. Configuraciones válidas: C(n, n_α)·C(n, n_β), una fracción mínima de los 2^{2n} bitstrings.

**Hamiltoniano:** Ĥ = E_nuc + Σ h_pq a†_p a_q + ½ Σ g_pqrs a†_p a†_r a_s a_q.
- Jordan-Wigner: a = (X + iY)/2, a_p = Z₀…Z_{p−1}(X_p + iY_p)/2, a†a = (I − Z)/2.
- Paridad: el qubit k guarda la paridad de 0..k. Se reducen 2 qubits (paridad α y paridad total ya conocidas → se reemplazan esas Z por ±1).
- Congelar orbitales de core (por ejemplo, O 1s): −2 qubits cada uno, con error pequeño. Congelar orbitales de valencia rompe el resultado.
- Agua STO-3G: 14 qubits → 10 (congelar 1 + reducir 2). Con simetrías se puede más (qubit tapering).

**Costo de medir:** ⟨H⟩ = Σ c_k⟨P_k⟩. Los términos que conmutan qubit a qubit se miden juntos: `H.group_commuting(qubit_wise=True)`.

### Funciones del curso (copiadas del notebook, PySCF + Qiskit)
```python
from pyscf import gto, scf, fci, mcscf, ao2mo, cc
from scipy.optimize import minimize
from math import comb

def molecula(atomos, base="sto-3g", n_congelados=0):
    mol = gto.M(atom=atomos, basis=base, verbose=0)
    hf = scf.RHF(mol).run()
    norb = mol.nao_nr() - n_congelados
    n_electrones = mol.nelectron - 2 * n_congelados
    activo = mcscf.CASCI(hf, norb, n_electrones)
    h, constante = activo.get_h1cas()
    g = ao2mo.restore(1, activo.get_h2cas(), norb)
    return {"hf": hf, "norb": norb, "nelec": (n_electrones // 2, n_electrones // 2),
            "h": h, "g": g, "constante": constante,
            "e_hf": hf.e_tot, "e_fci": activo.kernel()[0]}

def destruccion_jordan_wigner(p, n_qubits):
    return SparsePauliOp.from_sparse_list(
        [("Z" * p + "X", list(range(p + 1)), 0.5),
         ("Z" * p + "Y", list(range(p + 1)), 0.5j)], n_qubits)

def hamiltoniano_qubits(m, destruccion=destruccion_jordan_wigner):
    norb = m["norb"]; n_qubits = 2 * norb; h, g = m["h"], m["g"]
    a = [destruccion(p, n_qubits) for p in range(n_qubits)]
    a_dag = [op.adjoint() for op in a]
    salto = [[[a_dag[s*norb+p].dot(a[s*norb+q]).simplify() for q in range(norb)]
              for p in range(norb)] for s in range(2)]
    terminos = [SparsePauliOp("I" * n_qubits, m["constante"])]
    for s in range(2):
        for p in range(norb):
            for q in range(norb):
                terminos.append(h[p, q] * salto[s][p][q])
    for s in range(2):
        for t in range(2):
            for p in range(norb):
                for q in range(norb):
                    segundo = SparsePauliOp.sum([g[p, q, r, u] * salto[t][r][u]
                                                 for r in range(norb) for u in range(norb)])
                    terminos.append(0.5 * salto[s][p][q].dot(segundo))
                    if s == t:
                        for u in range(norb):
                            terminos.append(-0.5 * g[p, q, q, u] * salto[s][p][u])
    return SparsePauliOp.sum(terminos).simplify(atol=1e-10)

def destruccion_paridad(p, n_qubits):
    resto = list(range(p + 1, n_qubits))
    local = [("X", [0], 0.5), ("Y", [0], 0.5j)] if p == 0 else [("ZX", [p-1, p], 0.5), ("Y", [p], 0.5j)]
    return SparsePauliOp.from_sparse_list(
        [(l + "X" * len(resto), qs + resto, c) for l, qs, c in local], n_qubits)

def reducir_dos_qubits(H, m):
    norb = m["norb"]; n_alfa, n_beta = m["nelec"]
    conocidos = {norb - 1: (-1) ** n_alfa, 2 * norb - 1: (-1) ** (n_alfa + n_beta)}
    quedan = [k for k in range(2 * norb) if k not in conocidos]
    etiquetas, coefs = [], []
    for etiqueta, coef in H.to_list():
        letras = etiqueta[::-1]
        for k, v in conocidos.items():
            if letras[k] == "Z": coef = coef * v
        etiquetas.append("".join(letras[k] for k in quedan)[::-1]); coefs.append(coef)
    return SparsePauliOp(etiquetas, coefs).simplify()
```
Uso y verificaciones:
```python
m = molecula("H 0 0 0; H 0 0 0.735")
H = hamiltoniano_qubits(m)                                  # 4 qubits
np.linalg.eigvalsh(H.to_matrix())[0]                         # == m["e_fci"]
Statevector.from_label("0101").expectation_value(H).real     # == m["e_hf"]
H_red = reducir_dos_qubits(hamiltoniano_qubits(m, destruccion_paridad), m)   # 2 qubits
```

### VQE
Principio variacional: ⟨ψ|H|ψ⟩ ≥ E₀, así que toda energía obtenida es una cota superior.

Ciclo: ansatz |ψ(θ)⟩ → Estimator mide E(θ) → optimizador clásico propone nuevos θ → repetir.

```python
def vqe(ansatz, H, metodo="COBYLA", x0=None, estimator=None, maxiter=3000):
    estimator = estimator or StatevectorEstimator()
    x0 = np.zeros(ansatz.num_parameters) if x0 is None else x0
    historial = []
    def energia(p):
        v = float(estimator.run([(ansatz, H, p)]).result()[0].data.evs)
        historial.append(v); return v
    res = minimize(energia, x0, method=metodo, options={"maxiter": maxiter})
    return res, historial            # res.fun, res.x, res.nfev

def contar_cnots(c):
    return transpile(c, basis_gates=["cx", "rz", "sx", "x"], optimization_level=1).count_ops().get("cx", 0)
```

**Ansatz:**
- Hardware-efficient: `efficient_su2(n, reps=k)`. Circuito corto, pero difícil de optimizar; arrancar con parámetros chicos aleatorios (`np.random.default_rng(0).uniform(-0.1, 0.1, n_params)`).
- Químico UCCSD: X en los spin-orbitales ocupados (estado HF) y después, por cada salto simple o doble ocupado → vacío, `PauliEvolutionGate(i(T − T†), time=θ)`. Con θ = 0 es HF. Pocos parámetros con sentido físico, pero muchas CNOTs. Transpilar una sola vez al construirlo.

```python
def saltos_posibles(m):
    norb = m["norb"]; na, nb = m["nelec"]
    occ = [list(range(na)), [norb + i for i in range(nb)]]
    vac = [list(range(na, norb)), [norb + i for i in range(nb, norb)]]
    S = [((i,), (a,)) for s in range(2) for i in occ[s] for a in vac[s]]
    S += [((i, j), (a, b)) for i in occ[0] for j in occ[1] for a in vac[0] for b in vac[1]]
    for s in range(2):
        for x, i in enumerate(occ[s]):
            for j in occ[s][x+1:]:
                for y, a in enumerate(vac[s]):
                    for b in vac[s][y+1:]:
                        S.append(((i, j), (a, b)))
    return S

def ansatz_ucc(m):
    norb = m["norb"]; na, nb = m["nelec"]; nq = 2 * norb
    a = [destruccion_jordan_wigner(p, nq) for p in range(nq)]
    saltos = saltos_posibles(m); theta = ParameterVector("θ", len(saltos))
    c = QuantumCircuit(nq)
    for i in range(na): c.x(i)
    for i in range(nb): c.x(norb + i)
    for ang, (desde, hacia) in zip(theta, saltos):
        T = SparsePauliOp("I" * nq)
        for p in hacia: T = T.dot(a[p].adjoint())
        for p in reversed(desde): T = T.dot(a[p])
        G = (1j * (T - T.adjoint())).simplify()
        c.append(PauliEvolutionGate(G, time=ang), range(nq))
    return transpile(c, basis_gates=["cx", "rz", "sx", "x"], optimization_level=1)
```

**Optimizadores:** COBYLA y Powell usan solo valores; L-BFGS-B usa gradientes. Con ruido de shots, los de gradiente por diferencias finitas se pierden; preferir COBYLA o Powell. Ruido simulado: `AerEstimator(options={"default_precision": 0.01, "run_options": {"seed": 5}})` (≈10⁴ shots). Para comparar de forma justa, reevaluar `res.x` con el `StatevectorEstimator`.

**Curvas de energía:** usar los θ óptimos de la distancia anterior como `x0` (warm start).

### SQD: muestrear y diagonalizar
Los problemas de VQE son ansatz largos, optimización difícil y miles de evaluaciones. En SQD, el Sampler solo **señala configuraciones importantes**; después se diagonaliza clásicamente en ese subespacio (sigue siendo variacional). IBM lo usó en 2024 con hasta 77 qubits.

```python
def energia_en_subespacio(H, bitstrings, devolver_estado=False):
    M = H.to_matrix(sparse=True)
    idx = sorted({int(b, 2) for b in bitstrings})
    chica = M[idx][:, idx].toarray()
    w, v = np.linalg.eigh(chica)
    if not devolver_estado: return w[0]
    est = np.zeros(2 ** H.num_qubits, dtype=complex); est[idx] = v[:, 0]
    return w[0], est

def es_valida(b, m):
    norb = m["norb"]; na, nb = m["nelec"]
    return b[norb:].count("1") == na and b[:norb].count("1") == nb

def combinar_alfa_beta(bitstrings, m):
    n = m["norb"]
    alfas = {b[n:] for b in bitstrings}; betas = {b[:n] for b in bitstrings}
    return [be + al for be in betas for al in alfas]

# Flujo: circuito aproximado (por ejemplo, un VQE de 40 pasos) → medir → filtrar → combinar → diagonalizar
circ = ansatz.assign_parameters(res.x); circ.measure_all()
cuentas = AerSampler(seed=11).run([circ], shots=1000).result()[0].data.meas.get_counts()
validas = [b for b in cuentas if es_valida(b, m)]
E = energia_en_subespacio(H, combinar_alfa_beta(validas, m))
```
- Ruido de lectura: la fracción de bitstrings sin errores es (1 − p)^n (p = 0,02 y n = 100 → ~13%). En la práctica se usa *recuperación de configuraciones* (librería `qiskit-addon-sqd`).
- `get_bitstrings()` devuelve la lista completa de muestras, útil para inyectar ruido.

### Validar sin FCI
- **Cota variacional:** entre métodos variacionales, gana el de energía más baja. CCSD (`cc.CCSD(m["hf"]).run().e_tot`) NO es variacional y puede quedar por debajo de la exacta.
- **Convergencia:** el resultado se estabiliza al subir shots, subespacio, capas o base.
- **Varianza:** ⟨H²⟩ − ⟨H⟩² → 0 cerca de un autovector. Con ψ y M: ⟨H⟩ = ψ†Mψ y ⟨H²⟩ = ‖Mψ‖².
- **Comparar** con HF y CCSD y verificar el número de electrones.

---

## 6. Cómo presentar resultados
Para cada resultado, reportar:
- qué algoritmo, por qué y con qué alternativas se comparó;
- qubits, parámetros, profundidad y CNOTs del circuito transpilado;
- backend (exacto, shots o ruido) y cantidad de shots/evaluaciones;
- métrica contra la referencia (error en mHa contra 1,6 mHa, P(solución), QBER) y su chequeo de cordura;
- un gráfico (histograma, curva de convergencia o curva de energía) cuando ayude al pitch.
