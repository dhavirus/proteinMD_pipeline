// Presentation helpers for findings.json: grouping, filtering, evidence formatting.

export const SEVERITIES = ["blocking", "warn", "info"];
export const FAMILY_LABELS = {
  metals: "Metal sites",
  nonstandard_residues: "Non-standard residues",
  altlocs: "Alternative locations",
  missing_residues: "Missing residues",
  covalent_contacts: "Unannotated covalent contacts",
  unrecognized: "Unrecognized chemistry",
};
// Evidence items whose atoms are drawn as context around the selected finding.
const CONTEXT_EVIDENCE = ["ligand_distance", "attachment", "link", "contact_distance"];
const UNIT_SUFFIX = { angstrom: " Å", angstrom_squared: " Å²", degree: "°", fraction: "" };

export function severityCounts(findings) {
  const counts = Object.fromEntries(SEVERITIES.map((s) => [s, 0]));
  for (const finding of findings) counts[finding.effective_severity] += 1;
  return counts;
}

/**
 * Findings grouped by family (in report order), each group sorted by effective severity.
 * filter: {severities: string[] | null, family: string | null, undecidedOnly, decidedIds: Set}
 */
export function groupFindings(findings, filter) {
  const groups = new Map();
  for (const finding of findings) {
    if (!passes(finding, filter)) continue;
    if (!groups.has(finding.rule_family)) groups.set(finding.rule_family, []);
    groups.get(finding.rule_family).push(finding);
  }
  const rank = (f) => SEVERITIES.indexOf(f.effective_severity);
  return [...groups].map(([family, items]) => ({
    family,
    label: FAMILY_LABELS[family] || family,
    items: [...items].sort((a, b) => rank(a) - rank(b)),
  }));
}

function passes(finding, filter) {
  if (filter.severities && !filter.severities.includes(finding.effective_severity)) return false;
  if (filter.family && finding.rule_family !== filter.family) return false;
  return !(filter.undecidedOnly && filter.decidedIds.has(finding.id));
}

export function formatEvidenceValue(item) {
  if (item.type === "source_record") {
    return Object.entries(item.fields)
      .filter(([, value]) => value !== null)
      .map(([key, value]) => `${key}: ${value}`)
      .join(" · ");
  }
  if (item.type === "flag") return item.value ? "yes" : "no";
  const unit = item.unit in UNIT_SUFFIX ? UNIT_SUFFIX[item.unit] : item.unit ? ` ${item.unit}` : "";
  return `${item.value}${unit}`;
}

export function formatAtom(atom) {
  const altloc = atom.altloc ? `[${atom.altloc}]` : "";
  return `${atom.res_name} ${atom.chain}:${atom.seq_num}${atom.ins_code || ""} ${atom.atom_name}${altloc}`;
}

/** Altloc identifiers a "specific_altloc" decision may choose from (from the evidence). */
export function altlocChoices(finding) {
  const item = finding.evidence.find((e) => e.key === "altloc_ids");
  return item ? item.value.split(",").filter(Boolean) : [];
}

/** Residues to show for a finding: its anchors, plus context atoms named in the evidence. */
export function residuesToShow(finding) {
  const focus = finding.anchor_residues.length ? finding.anchor_residues : finding.locus.extent;
  const context = finding.evidence
    .filter((e) => CONTEXT_EVIDENCE.includes(e.key))
    .flatMap((e) => e.atoms || [])
    .map((a) => ({ chain: a.chain, seq_num: a.seq_num, ins_code: a.ins_code }));
  return { focus, context };
}
