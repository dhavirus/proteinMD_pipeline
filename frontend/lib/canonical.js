// RFC 8785 (JCS) canonical JSON and SHA-256, identical to simprep.canonical in Python.
// JSON.stringify already writes numbers and strings the way JCS requires; only object
// keys need sorting (by UTF-16 code units, which is JavaScript's default string order).

export function canonicalJson(value) {
  if (value === null || typeof value !== "object") {
    if (typeof value === "number" && !Number.isFinite(value)) {
      throw new TypeError(`${value} is not representable in JSON`);
    }
    return JSON.stringify(value);
  }
  if (Array.isArray(value)) return `[${value.map(canonicalJson).join(",")}]`;
  const keys = Object.keys(value).sort();
  return `{${keys.map((key) => `${JSON.stringify(key)}:${canonicalJson(value[key])}`).join(",")}}`;
}

export async function sha256Hex(data) {
  const bytes = typeof data === "string" ? new TextEncoder().encode(data) : data;
  const digest = await globalThis.crypto.subtle.digest("SHA-256", bytes);
  return [...new Uint8Array(digest)].map((b) => b.toString(16).padStart(2, "0")).join("");
}

export const sha256Canonical = (value) => sha256Hex(canonicalJson(value));
