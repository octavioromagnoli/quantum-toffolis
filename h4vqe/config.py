"""Configuración compartida. Cambiar acá, no en los módulos.

Dueño: todos (avisar en el grupo antes de tocar).
"""

BASIS = "sto-3g"

# Distancias H-H (Å). La consigna pide >= 10 puntos entre 0,6 y 3,0 Å.
# Hay puntos densos cerca del mínimo (~0,9 Å) para ajustar R_e.
DISTANCES = [0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.4, 1.6,
             1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0]

STRETCHED_FROM = 1.8          # zona donde hay que ganarle a CCSD
CHEM_ACC = 1.6e-3             # Ha

# Presupuesto de shots por punto de la curva (estimación final, no optimización).
SHOTS_PER_POINT = 8_000_000   # σ ≈ 0,5 mHa en la zona estirada (QWC con sorted insertion + Neyman)

# Compuertas nativas que fija la consigna para el extra con ruido: cz, sx, x y rz.
# (ruido.py no se distribuye; el modelo de ruido propio está en noise.py.)
BASIS_GATES = ["cz", "rz", "sx", "x"]

SEED = 1234

# Ansatz de la curva principal:
#   "adapt": ADAPT-VQE fermiónico, lista común de 15 operadores (adapt.ADAPT_FERMIONIC).
#            194 CZ todos con todos, 378 en línea con el ruteo fSWAP de routing.py.
#            Error sin ruido ≤ 0,38 mHa en toda la curva.
#   "uccsd": UCCSD con ANSATZ_REPS pasos de Trotter (582 CZ con 2 pasos; < 1 µHa sin ruido).
ANSATZ = "adapt"

# Pasos de Trotter de UCCSD. Con ANSATZ = "adapt", UCCSD ×2 queda como referencia más
# profunda (prueba de profundidad, comparación de ansätze, geometría y FakeFez).
# Con 1 paso el error máximo es ~4,5 mHa (a 2,0 Å): le gana a CCSD pero sin precisión química.
ANSATZ_REPS = 2
