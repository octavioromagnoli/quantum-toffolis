"""Calibración de σ del estimador (issue #3, Santino).

Repite measure_energy muchas veces con semillas independientes sobre el MISMO
circuito y mira z = (E − E_exacta)/σ. Si σ está bien calibrada, z ~ N(0, 1).
E_exacta = <ψ|H|ψ> del circuito medido (statevector), solo para calibrar.

Uso:
    python scripts/calib_sigma.py --quick        # prueba de ~1 min
    python scripts/calib_sigma.py                # corrida completa (PC de escritorio)
    python scripts/calib_analyze.py              # histogramas y tabla

Cada repetición se guarda apenas termina en results/calib_sigma/runs.jsonl.
Si se corta, relanzar el mismo comando: retoma donde quedó.
"""
from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")  # un hilo por proceso; paralelizamos con procesos

import argparse
import json
import pickle
import sys
import time
import zlib
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")  # consola/log de Windows
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
OUT = os.path.join(ROOT, "results", "calib_sigma")

# (nombre, shots, pilot_fraction, repeticiones, agrupamiento). Cada config en cada distancia.
CONFIGS = [
    ("neyman_8M", 8_000_000, 0.1, 500, "qwc"),  # lo que usa el pipeline hoy
    ("sumc_8M", 8_000_000, 1.0, 300, "qwc"),    # una sola etapa N_g ∝ Σ|c| (¿el piloto sesga σ?)
    ("neyman_1M", 1_000_000, 0.1, 500, "qwc"),  # z no debería depender de los shots
    ("gc_4M", 4_000_000, 0.1, 500, "gc"),       # conmutación general + Clifford
]
DISTANCES = [2.0, 3.0]
N_BOOT = 200


def prep(d: float) -> dict:
    """Hamiltoniano, θ* y circuito ligado para la distancia d (cacheado en disco)."""
    path = os.path.join(OUT, f"prep_d{d:.2f}.pkl")
    if os.path.exists(path):
        with open(path, "rb") as f:
            return pickle.load(f)
    from qiskit.quantum_info import Statevector
    from h4vqe.ansatz import build_ansatz, cz_count
    from h4vqe.config import ANSATZ_REPS
    from h4vqe.hamiltonian import molecule, qubit_hamiltonian
    from h4vqe.optimize import optimize_params
    m = molecule(d)
    H = qubit_hamiltonian(m)
    ans = build_ansatz(m, reps=ANSATZ_REPS)
    opt = optimize_params(ans, H, n_restarts=2)
    bound = ans.circuit.assign_parameters(opt["x"])
    p = {"d": d, "circuit": bound, "H": H,
         "E_exact": float(Statevector(bound).expectation_value(H).real),
         "E_fci": m["e_fci"], "cz": cz_count(ans.circuit), "ansatz_reps": ANSATZ_REPS}
    with open(path, "wb") as f:
        pickle.dump(p, f)
    return p


_PLANS = {}


def _init_worker(preps, groupings):
    from h4vqe.estimator import prepare_measurement
    for d, p in preps.items():
        for gr in groupings:
            _PLANS[d, gr] = (p, prepare_measurement(p["circuit"], p["H"], grouping=gr))


def _job(cfg, d, shots, pilot, grouping, rep, seed):
    from h4vqe.estimator import measure_energy
    p, plan = _PLANS[d, grouping]
    t0 = time.time()
    r = measure_energy(None, None, shots, seed=seed, pilot_fraction=pilot, plan=plan, n_boot=N_BOOT)
    E0 = p["E_exact"]
    return {"cfg": cfg, "d": d, "rep": rep, "seed": seed, "shots_req": shots, "pilot": pilot,
            "grouping": grouping,
            "E": r["energy"], "sigma": r["sigma"], "sigma_boot": r["sigma_boot"],
            "shots": r["shots"], "n_groups": r["n_groups"], "E_exact": E0,
            "z": (r["energy"] - E0) / r["sigma"], "z_boot": (r["energy"] - E0) / r["sigma_boot"],
            "seconds": time.time() - t0}


def job_seed(cfg: str, d: float, rep: int) -> int:
    """Semilla de 62 bits, independiente para cada (config, d, rep). measure_energy la
    expande con SeedSequence en flujos independientes por etapa y por grupo."""
    ss = np.random.SeedSequence([20261121, zlib.crc32(cfg.encode()), int(round(d * 1000)), rep])
    return int(ss.generate_state(1, np.uint64)[0] >> np.uint64(2))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) // 2),
                    help="procesos en paralelo (default: núcleos físicos aprox.)")
    ap.add_argument("--quick", action="store_true", help="prueba corta (pocas reps, pocos shots)")
    ap.add_argument("--distances", type=float, nargs="+", default=DISTANCES)
    ap.add_argument("--configs", nargs="+", default=None, help="subconjunto de nombres de CONFIGS")
    ap.add_argument("--reps-scale", type=float, default=1.0, help="multiplica las repeticiones")
    args = ap.parse_args()

    os.makedirs(OUT, exist_ok=True)
    configs = [c for c in CONFIGS if args.configs is None or c[0] in args.configs]
    runs_path = os.path.join(OUT, "runs.jsonl")
    if args.quick:
        configs = [("quick_200k", 200_000, 0.1, 4, "qwc"), ("quick_gc_200k", 200_000, 0.1, 4, "gc")]
        runs_path = os.path.join(OUT, "runs_quick.jsonl")

    t0 = time.time()
    preps = {d: prep(d) for d in args.distances}
    for d, p in preps.items():
        print(f"d={d:.2f} Å  E_exacta={p['E_exact']:.8f}  E_exacta−FCI={1e3*(p['E_exact']-p['E_fci']):.4f} mHa"
              f"  CZ={p['cz']}", flush=True)
    print(f"preparación: {time.time()-t0:.0f} s", flush=True)

    done = set()
    if os.path.exists(runs_path):
        with open(runs_path) as f:
            for line in f:
                try:
                    r = json.loads(line)
                    done.add((r["cfg"], r["d"], r["rep"]))
                except json.JSONDecodeError:
                    pass  # línea cortada por una interrupción
    # Orden intercalado: si se corta antes, todas las configs tienen ~la misma cantidad.
    nmax = max(int(c[3] * args.reps_scale) for c in configs)
    jobs = [(name, d, shots, pilot, grouping, rep, job_seed(name, d, rep))
            for rep in range(nmax) for name, shots, pilot, reps, grouping in configs
            if rep < int(reps * args.reps_scale) for d in args.distances
            if (name, d, rep) not in done]
    print(f"{len(done)} hechas, {len(jobs)} pendientes, {args.workers} procesos", flush=True)
    if not jobs:
        return

    t0, n = time.time(), 0
    with ProcessPoolExecutor(args.workers, initializer=_init_worker,
                            initargs=(preps, sorted({c[4] for c in configs}))) as ex, \
            open(runs_path, "a") as fout:
        futs = [ex.submit(_job, *j) for j in jobs]
        for fut in as_completed(futs):
            r = fut.result()
            fout.write(json.dumps(r) + "\n")
            fout.flush()
            n += 1
            if n % 10 == 0 or n == len(jobs):
                el = time.time() - t0
                print(f"[{time.strftime('%H:%M:%S')}] {n}/{len(jobs)}  "
                      f"{el/60:.1f} min  faltan ~{el/n*(len(jobs)-n)/60:.0f} min  "
                      f"(último: {r['cfg']} d={r['d']} z={r['z']:+.2f})", flush=True)
    print("listo. Ahora: python scripts/calib_analyze.py", flush=True)


if __name__ == "__main__":
    main()
