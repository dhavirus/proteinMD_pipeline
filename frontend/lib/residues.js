// Residue addresses as simprep writes them: chain:number[insertion code], author numbering.

const TOKEN = /^([A-Za-z0-9]{1,4}):(-?\d+)([A-Za-z]?)(?:-(-?\d+))?$/;
export const MAX_RANGE_LENGTH = 2000;

export function residueLabel(ref) {
  return `${ref.chain}:${ref.seq_num}${ref.ins_code || ""}`;
}

export function sameResidue(a, b) {
  return a.chain === b.chain && a.seq_num === b.seq_num && (a.ins_code || "") === (b.ins_code || "");
}

/**
 * Parse "A:45, A:52A, A:444-453" into residue refs. Returns {residues, errors};
 * duplicates are dropped, order of first appearance is kept.
 */
export function parseResidueList(text) {
  const residues = [];
  const errors = [];
  for (const token of text.split(/[\s,;]+/).filter(Boolean)) {
    const parsed = parseToken(token);
    if (parsed.error) errors.push(parsed.error);
    for (const ref of parsed.residues || []) {
      if (!residues.some((seen) => sameResidue(seen, ref))) residues.push(ref);
    }
  }
  return { residues, errors };
}

function parseToken(token) {
  const match = TOKEN.exec(token);
  if (!match) return { error: `"${token}" is not a residue (expected e.g. A:45, A:52A or A:444-453)` };
  const [, chain, start, insCode, end] = match;
  if (end === undefined) return { residues: [{ chain, seq_num: Number(start), ins_code: insCode }] };
  if (insCode) return { error: `"${token}": ranges cannot carry an insertion code` };
  const [first, last] = [Number(start), Number(end)];
  if (last < first) return { error: `"${token}": range end is before its start` };
  if (last - first + 1 > MAX_RANGE_LENGTH) return { error: `"${token}": range longer than ${MAX_RANGE_LENGTH} residues` };
  const residues = [];
  for (let n = first; n <= last; n += 1) residues.push({ chain, seq_num: n, ins_code: "" });
  return { residues };
}

export function formatResidueList(residues) {
  return residues.map(residueLabel).join(", ");
}
