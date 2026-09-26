// Schema validation in the browser with ajv (pinned CDN build), against the repository's
// own schema/*.schema.json files, the same files simprep validates with in Python.

const SCHEMA_NAMES = ["finding", "rule", "manifest", "findings_report"];

export async function createValidator(schemaBaseUrl) {
  const library = globalThis.ajv2020;
  if (!library) throw new Error("The schema validator (ajv) did not load; check your connection and reload.");
  const Ajv = library.default || library;
  // Formats are not asserted, matching the Python validator (jsonschema without a format checker).
  const ajv = new Ajv({ strict: false, validateFormats: false, allErrors: true });
  for (const name of SCHEMA_NAMES) {
    const response = await fetch(new URL(`${name}.schema.json`, schemaBaseUrl));
    if (!response.ok) throw new Error(`Could not load schema ${name} (HTTP ${response.status}).`);
    ajv.addSchema(await response.json());
  }
  return (name, document_) => {
    const valid = ajv.validate(`urn:simprep:schema:${name}`, document_);
    return valid ? [] : ajv.errors.slice(0, 20).map((e) => `${e.instancePath || "<root>"}: ${e.message}`);
  };
}
