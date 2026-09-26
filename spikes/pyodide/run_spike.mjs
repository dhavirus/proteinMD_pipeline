// ADR-0002 spike driver: load Pyodide in Node, install the simprep pure core from the
// working tree into Pyodide's filesystem, and time an audit of exported inputs.
// Usage: node run_spike.mjs INPUTS.json [FINDINGS_OUT.json]
// (INPUTS.json from export_inputs.py; FINDINGS_OUT.json for the parity check)
import { readFileSync, readdirSync, statSync, writeFileSync } from "node:fs";
import { dirname, join, relative } from "node:path";
import { fileURLToPath } from "node:url";
import { loadPyodide } from "pyodide";

const here = dirname(fileURLToPath(import.meta.url));
const srcRoot = join(here, "..", "..", "src");
const mb = (bytes) => Math.round(bytes / 1e6);

function pythonFiles(dir) {
  return readdirSync(dir).flatMap((name) => {
    const path = join(dir, name);
    if (statSync(path).isDirectory()) return pythonFiles(path);
    return path.endsWith(".py") ? [path] : [];
  });
}

const t0 = performance.now();
const pyodide = await loadPyodide();
const t1 = performance.now();
for (const path of pythonFiles(join(srcRoot, "simprep"))) {
  const target = "/lib/" + relative(srcRoot, path);
  pyodide.FS.mkdirTree(dirname(target));
  pyodide.FS.writeFile(target, readFileSync(path, "utf8"));
}
pyodide.FS.writeFile("/lib/run_in_pyodide.py", readFileSync(join(here, "run_in_pyodide.py"), "utf8"));
pyodide.runPython("import sys; sys.path.insert(0, '/lib')");
const t2 = performance.now();
pyodide.runPython("import run_in_pyodide");
const t3 = performance.now();
const inputs = readFileSync(process.argv[2], "utf8");
pyodide.globals.set("inputs_json", inputs);
const result = JSON.parse(pyodide.runPython("run_in_pyodide.main(inputs_json)"));
const t4 = performance.now();
if (process.argv[3]) writeFileSync(process.argv[3], JSON.stringify(result.findings));
const memory = process.memoryUsage();
console.log(JSON.stringify({
  pyodide_version: pyodide.version,
  load_pyodide_ms: Math.round(t1 - t0),
  install_sources_ms: Math.round(t2 - t1),
  import_simprep_ms: Math.round(t3 - t2),
  rebuild_ms: result.rebuild_ms,
  audit_ms: result.audit_ms,
  end_to_end_ms: Math.round(t4 - t0),
  inputs_mb: mb(inputs.length),
  atoms: result.atoms,
  findings: result.findings.length,
  node_rss_mb: mb(memory.rss),
  wasm_heap_mb: mb(pyodide._module.HEAPU8.length),
}, null, 2));
