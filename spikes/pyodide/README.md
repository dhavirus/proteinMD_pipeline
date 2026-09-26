# ADR-0002 spike: simprep in the browser via Pyodide

Measures whether the simprep detector core runs under Pyodide, using Node as the
WebAssembly host (same Pyodide build as the browser; timings indicative only).

```bash
cd spikes/pyodide
npm install                                   # pyodide 314.0.7 from the npm registry
python export_inputs.py ../../tests/panel/5FQL.cif.gz /tmp/5fql.json   # CPython + gemmi
node run_spike.mjs /tmp/5fql.json /tmp/5fql_pyodide.json               # Pyodide
python check_parity.py ../../tests/panel/5FQL.cif.gz /tmp/5fql_pyodide.json
```

`export_inputs.py` parses with gemmi in CPython and writes the parser-independent
structure model plus the knowledge base as JSON. `run_spike.mjs` loads Pyodide, copies
`src/simprep` into its filesystem, rebuilds the model, and runs `run_audit`.
`check_parity.py` confirms the Pyodide findings equal a CPython run.

Results and the recommendation are in `docs/decisions/ADR-0002-pyodide-spike.md`.
