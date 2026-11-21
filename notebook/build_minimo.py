"""Genera notebook/H4_VQE_minimo.ipynb: solo lo que pide la consigna, sin extras.

    python notebook/build_minimo.py            # genera
    python notebook/build_minimo.py --run      # genera y ejecuta en una carpeta vacía (como un Colab nuevo)
    python notebook/build_minimo.py --text-only  # actualiza solo el texto del notebook ya ejecutado

Usa las mismas celdas de código y los mismos .md que build_notebook.py (no los modifica).
Diferencias con el notebook completo:
  - sin la sección "Extras" (07_extras.md); la sección de decisiones pasa a ser la 6;
  - sin las líneas "> TODO" de los .md;
  - sin la referencia a la sección de mitigación (no está en el mínimo);
  - con un párrafo de limitaciones al final de la sección de decisiones.
"""
from __future__ import annotations

import pathlib
import shutil
import sys
import tempfile
import time

import nbformat as nbf

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import build_notebook as full  # noqa: E402

OUT = full.ROOT / "notebook" / "H4_VQE_minimo.ipynb"

LIMITS = """
**Limitaciones**
- **Compuertas:** UCCSD ×2 usa 582 CZ con conectividad total y 1648 con qubits en línea (columnas `cz_all` y `cz_line` de la tabla). Sin ruido no importa, pero con el 0,1 % de error por CZ del modelo de la consigna, la probabilidad de que no falle ninguna CZ sería 0,999^582 ≈ 56 % (todos con todos) y 0,999^1648 ≈ 19 % (en línea). En hardware actual este circuito no serviría sin bajar las CZ o mitigar.
- **Optimización sin ruido:** los parámetros se buscan con el simulador exacto, como permite la consigna. La curva mide qué tan expresivo es el ansatz más el ruido de los shots, no la dificultad de optimizar con mediciones.
- **σ es solo estadística:** cubre la dispersión de los shots (calibrada en la sección 3), no sesgos de ruido de compuertas o de lectura, que en esta simulación no hay.
- **Diagnósticos de la sección 5:** la varianza, las simetrías y la comparación con 1 paso de Trotter se calculan con el estado del simulador. En hardware habría que medirlos con shots.
- **Costo en shots:** 8 M por punto (144 M en las 18 distancias) para σ ≈ 0,5–0,6 mHa.
"""

# (archivo, texto original, reemplazo). Si el texto original no está, falla: hay que revisar el .md.
REPLACE = [
    ("04_measurement.md", "y el sesgo se trata en la sección de mitigación.",
     "y no incluye el sesgo que agregaría el ruido."),
    ("08_decisions.md", "## 7 · ", "## 6 · "),
]


def md_text(name: str) -> str:
    text = (full.MD / name).read_text(encoding="utf-8")
    for fname, old, new in REPLACE:
        if fname == name:
            if old not in text:
                raise ValueError(f"{name}: no encontré {old!r}")
            text = text.replace(old, new)
    lines = [ln for ln in text.splitlines() if not ln.lstrip().startswith("> TODO")]
    text = "\n".join(lines).strip() + "\n"
    if name == "08_decisions.md":
        text += LIMITS
    return text


def md(name: str):
    return nbf.v4.new_markdown_cell(md_text(name))


def build() -> nbf.NotebookNode:
    nb = nbf.v4.new_notebook()
    cells = [md("00_intro.md"), md("01_setup.md"), nbf.v4.new_code_cell(full.INSTALL)]
    for mod in full.MODULES:
        src = (full.PKG / f"{mod}.py").read_text(encoding="utf-8")
        cells.append(nbf.v4.new_code_cell(f"%%writefile h4vqe/{mod}.py\n{src}"))
    cells += [md("02_hamiltonian.md"), md("03_ansatz.md"), md("04_measurement.md"),
              nbf.v4.new_code_cell(full.CHECKS),
              md("05_results.md"), nbf.v4.new_code_cell(full.RUN), nbf.v4.new_code_cell(full.TABLE),
              nbf.v4.new_code_cell(full.PLOT), nbf.v4.new_code_cell(full.RE_DE),
              md("06_reliability.md"), nbf.v4.new_code_cell(full.RELIABILITY),
              md("08_decisions.md")]
    nb.cells = cells
    nb.metadata["kernelspec"] = {"name": "python3", "display_name": "Python 3", "language": "python"}
    nb.metadata["colab"] = {"provenance": []}
    leftover = [c.source for c in cells if c.cell_type == "markdown" and "TODO" in c.source]
    assert not leftover, f"quedó un TODO en el texto: {leftover}"
    return nb


def run(nb: nbf.NotebookNode, kernel: str = "python3") -> nbf.NotebookNode:
    """Ejecuta en una carpeta temporal vacía: %%writefile crea h4vqe/ ahí, como en Colab."""
    from nbclient import NotebookClient

    workdir = pathlib.Path(tempfile.mkdtemp(prefix="h4vqe_minimo_"))
    try:
        t0 = time.time()
        NotebookClient(nb, timeout=None, kernel_name=kernel,
                       resources={"metadata": {"path": str(workdir)}}).execute()
        print(f"ejecutado OK en {time.time() - t0:.0f} s")
    finally:
        shutil.rmtree(workdir, ignore_errors=True)
    return nb


def update_text(path: pathlib.Path = OUT) -> None:
    """Reemplaza solo las celdas de texto del notebook ya ejecutado (las salidas quedan).
    Sirve para corregir la redacción sin volver a correr todo."""
    old, new = nbf.read(path, as_version=4), build()
    kinds = lambda nb: [c.cell_type for c in nb.cells]  # noqa: E731
    if kinds(old) != kinds(new):
        raise ValueError("cambió la estructura de celdas: hay que volver a ejecutar con --run")
    for c_old, c_new in zip(old.cells, new.cells):
        if c_new.cell_type == "markdown":
            c_old.source = c_new.source
        elif c_old.source != c_new.source:
            raise ValueError("cambió una celda de código: hay que volver a ejecutar con --run")
    nbf.write(old, path)
    print("texto actualizado en", path)


if __name__ == "__main__":
    if "--text-only" in sys.argv:
        update_text()
        sys.exit()
    nb = build()
    nbf.write(nb, OUT)
    print("escrito", OUT)
    if "--run" in sys.argv:
        nbf.write(run(nb), OUT)
        print("escrito (ejecutado)", OUT)
