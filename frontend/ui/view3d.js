// 3D view: Mol* driven through MolViewSpec. Only displays; never changes coordinates.
import { residuesToShow } from "../lib/findings.js";

const COLORS = { protein: "#9aa7a3", other: "#7d8b86", focus: "#e0a33a", context: "#d4dad8", region: "#4fb3bb" };
const BACKGROUND = "#10161a";
// Highlights are drawn slightly larger than the generic ligand/ion layer so they are not
// hidden inside it; the camera keeps FOCUS_MARGIN_ANGSTROM of context around them.
const SIZE = { generic: 0.8, context: 0.95, focus: 1.25 };
const FOCUS_MARGIN_ANGSTROM = 8;

const selector = (refs) =>
  refs.map((r) => ({ auth_asym_id: r.chain, auth_seq_id: r.seq_num, ...(r.ins_code ? { pdbx_PDB_ins_code: r.ins_code } : {}) }));

export async function createView(container) {
  const molstar = globalThis.molstar;
  if (!molstar || !molstar.Viewer) throw new Error("The 3D viewer (Mol*) did not load; findings can still be reviewed.");
  const viewer = await molstar.Viewer.create(container, {
    layoutIsExpanded: false, layoutShowControls: false, layoutShowRemoteState: false,
    layoutShowSequence: false, layoutShowLog: false, layoutShowLeftPanel: false,
    viewportShowExpand: true, viewportShowSelectionMode: false, viewportShowAnimation: false,
    viewportBackgroundColor: BACKGROUND,
  });
  let structureUrl = null;
  let format = "mmcif";
  return {
    setStructure(text, structureFormat) {
      if (structureUrl) URL.revokeObjectURL(structureUrl);
      structureUrl = URL.createObjectURL(new Blob([text], { type: "text/plain" }));
      format = structureFormat;
    },
    async show(finding, regions) {
      if (!structureUrl) return;
      const builder = molstar.lib.extensions.mvs.createBuilder();
      const structure = builder.download({ url: structureUrl }).parse({ format }).modelStructure({});
      structure.component({ selector: "polymer" }).representation({ type: "cartoon" }).color({ color: COLORS.protein });
      for (const kind of ["ligand", "branched", "ion"]) {
        structure.component({ selector: kind }).representation({ type: "ball_and_stick", size_factor: SIZE.generic }).color({ color: COLORS.other });
      }
      const regionResidues = regions.flatMap((r) => r.residues);
      if (regionResidues.length) {
        structure.component({ selector: selector(regionResidues) }).representation({ type: "ball_and_stick", size_factor: SIZE.context }).color({ color: COLORS.region });
      }
      if (finding) addFinding(structure, finding);
      builder.canvas({ background_color: BACKGROUND });
      await viewer.loadMvsData(builder.getState(), "mvsj", { replaceExisting: true });
    },
  };
}

function addFinding(structure, finding) {
  const { focus, context } = residuesToShow(finding);
  if (context.length) {
    structure.component({ selector: selector(context) }).representation({ type: "ball_and_stick", size_factor: SIZE.context }).color({ color: COLORS.context });
  }
  if (focus.length) {
    const component = structure.component({ selector: selector(focus) });
    component.representation({ type: "ball_and_stick", size_factor: SIZE.focus }).color({ color: COLORS.focus });
    component.focus({ radius_extent: FOCUS_MARGIN_ANGSTROM });
  }
}
