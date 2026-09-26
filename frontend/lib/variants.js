// Variants as recorded in manifest.schema.json (variant, mutation), from typed text such
// as "A:468 R>Q" or "A:468 ARG>GLN" (pure). simprep variants re-checks everything.

export const THREE_LETTER = {
  A: "ALA", R: "ARG", N: "ASN", D: "ASP", C: "CYS", Q: "GLN", E: "GLU", G: "GLY", H: "HIS", I: "ILE",
  L: "LEU", K: "LYS", M: "MET", F: "PHE", P: "PRO", S: "SER", T: "THR", W: "TRP", Y: "TYR", V: "VAL",
};
const ONE_LETTER = Object.fromEntries(Object.entries(THREE_LETTER).map(([one, three]) => [three, one]));
const MUTATION = /^([A-Za-z0-9]+):(-?\d+)([A-Za-z]?)\s+([A-Za-z]{1,3})\s*>\s*([A-Za-z]{1,3})$/;

const residueName = (code) => {
  const upper = code.toUpperCase();
  return upper.length === 1 ? THREE_LETTER[upper] : (ONE_LETTER[upper] ? upper : null);
};

/** {mutation, error}: one mutation from "A:468 R>Q"; error explains what is wrong. */
export function parseMutation(text) {
  const match = MUTATION.exec(text.trim());
  if (!match) return { error: `"${text.trim()}" is not a mutation like "A:468 R>Q".` };
  const [, chain, number, icode, fromCode, toCode] = match;
  const from = residueName(fromCode);
  const to = residueName(toCode);
  if (!from || !to) return { error: `"${text.trim()}": use one-letter or three-letter standard amino-acid codes.` };
  if (from === to) return { error: `"${text.trim()}": from and to are the same residue.` };
  return { mutation: { chain, seq_num: Number(number), ins_code: icode.toUpperCase(), from, to } };
}

/** Conventional name, e.g. R468Q (joined with "_" for several mutations). */
export function variantName(mutations) {
  return mutations.map((m) => `${ONE_LETTER[m.from]}${m.seq_num}${m.ins_code}${ONE_LETTER[m.to]}`).join("_");
}

/**
 * {variant, errors} from a draft {mutations: "A:468 R>Q, ...", rationale}; residues is
 * the Map from residueNames() of the loaded structure; otherNames are the names in use.
 */
export function variantFromDraft(draft, residues, otherNames) {
  const parsed = draft.mutations.split(/[,;\n]/).map((t) => t.trim()).filter(Boolean).map(parseMutation);
  const errors = parsed.filter((p) => p.error).map((p) => p.error);
  const mutations = parsed.filter((p) => p.mutation).map((p) => p.mutation);
  if (!parsed.length) errors.push("Enter at least one mutation, e.g. A:468 R>Q.");
  for (const m of mutations) errors.push(...residueProblems(m, residues));
  const name = variantName(mutations);
  if (mutations.length && otherNames.includes(name)) errors.push(`A variant named ${name} already exists.`);
  if (!(draft.rationale || "").trim()) errors.push("Rationale: say why this variant is part of the study.");
  if (errors.length) return { errors };
  return { errors: [], variant: { name, mutations, rationale: draft.rationale.trim(), references: [] } };
}

function residueProblems(mutation, residues) {
  const key = `${mutation.chain}:${mutation.seq_num}${mutation.ins_code}`;
  const found = residues.get(key);
  if (!found) return [`${key} is not a residue of the loaded structure.`];
  if (found !== mutation.from) return [`${key} is ${found} in the loaded structure, not ${mutation.from}.`];
  return [];
}

export const formatMutation = (m) => `${m.chain}:${m.seq_num}${m.ins_code} ${m.from}>${m.to}`;
