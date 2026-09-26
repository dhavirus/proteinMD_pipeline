import { readFileSync } from "node:fs";
import { fileURLToPath } from "node:url";

const root = new URL("../../", import.meta.url);
export const readJson = (relative) => JSON.parse(readFileSync(fileURLToPath(new URL(relative, root)), "utf8"));
export const miniReport = () => readJson("schema/examples/findings_report/valid_mini.json");
export const exampleManifest = () => readJson("schema/examples/manifest/valid_with_decision.json");
