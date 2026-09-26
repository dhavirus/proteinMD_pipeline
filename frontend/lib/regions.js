// Regions of interest: named residue sets as in manifest.schema.json (region).
import { parseResidueList, residueLabel } from "./residues.js";

const NAME = /^[A-Za-z0-9_-]+$/;

/** Validate a region form; returns {region, errors}. */
export function regionFromDraft(draft, otherNames) {
  const errors = [];
  const name = (draft.name || "").trim();
  const description = (draft.description || "").trim();
  if (!NAME.test(name)) errors.push("Name: letters, digits, _ and - only, no spaces.");
  if (otherNames.includes(name)) errors.push(`Name "${name}" is already used by another region.`);
  if (!description) errors.push("Description: say what the region is and why it matters.");
  const { residues, errors: residueErrors } = parseResidueList(draft.residues || "");
  errors.push(...residueErrors);
  if (!residues.length && !residueErrors.length) errors.push("Residues: list at least one residue.");
  return { region: errors.length ? null : { name, description, residues }, errors };
}

/** Order-insensitive comparison of two region lists (used for the stale-severity flag). */
export function sameRegions(a, b) {
  return normalize(a) === normalize(b);
}

function normalize(regions) {
  return JSON.stringify(
    regions
      .map((r) => [r.name, r.residues.map(residueLabel).sort()])
      .sort(([x], [y]) => (x < y ? -1 : x > y ? 1 : 0)),
  );
}
