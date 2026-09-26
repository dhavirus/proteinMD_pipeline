// Residue names per chain:number from a structure's coordinate records (pure; no Mol*).
// mmCIF: the _atom_site loop (model 1; author chain, number, insertion code and component).
// PDB: ATOM/HETATM fixed columns. Used to check typed mutations against the loaded file.

const PDB_RECORD = /^(ATOM  |HETATM)/;

/** Map "A:468" (with insertion code appended, e.g. "A:52A") -> component id. */
export function residueNames(text, format) {
  return format === "pdb" ? pdbResidues(text) : mmcifResidues(text);
}

function pdbResidues(text) {
  const names = new Map();
  for (const line of text.split(/\r?\n/)) {
    if (line.startsWith("ENDMDL")) break;
    if (!PDB_RECORD.test(line)) continue;
    const key = `${line.slice(21, 22).trim()}:${line.slice(22, 26).trim()}${line.slice(26, 27).trim()}`;
    if (!names.has(key)) names.set(key, line.slice(17, 20).trim());
  }
  return names;
}

function mmcifResidues(text) {
  const lines = text.split(/\r?\n/);
  const start = lines.findIndex((line, i) => line.trim() === "loop_" && lines[i + 1]?.startsWith("_atom_site."));
  if (start < 0) return new Map();
  const columns = [];
  let row = start + 1;
  while (lines[row]?.startsWith("_atom_site.")) columns.push(lines[row++].trim().slice("_atom_site.".length));
  const at = (name) => columns.indexOf(name);
  const index = {
    chain: at("auth_asym_id"), seq: at("auth_seq_id"), icode: at("pdbx_PDB_ins_code"),
    comp: at("auth_comp_id") >= 0 ? at("auth_comp_id") : at("label_comp_id"), model: at("pdbx_PDB_model_num"),
  };
  return collect(lines, row, columns.length, index);
}

function collect(lines, row, width, index) {
  const names = new Map();
  let firstModel = null;
  for (let i = row; i < lines.length; i += 1) {
    const line = lines[i];
    if (line.startsWith("_") || line.startsWith("loop_") || line.startsWith("#")) break;
    const values = tokens(line);
    if (values.length < width) continue;
    const model = index.model >= 0 ? values[index.model] : "1";
    firstModel ??= model;
    if (model !== firstModel) break;
    const icode = index.icode >= 0 && !["?", "."].includes(values[index.icode]) ? values[index.icode] : "";
    const key = `${values[index.chain]}:${values[index.seq]}${icode}`;
    if (!names.has(key)) names.set(key, values[index.comp]);
  }
  return names;
}

/** Whitespace-separated CIF values; quoted values may contain spaces. */
export function tokens(line) {
  const values = [];
  const pattern = /'([^']*)'(?=\s|$)|"([^"]*)"(?=\s|$)|(\S+)/g;
  for (const match of line.matchAll(pattern)) values.push(match[1] ?? match[2] ?? match[3]);
  return values;
}
