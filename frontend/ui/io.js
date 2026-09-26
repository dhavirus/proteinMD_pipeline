// Local file access: nothing is uploaded; files are read in the browser.
import { sha256Hex } from "../lib/canonical.js";

const GZIP_MAGIC = [0x1f, 0x8b];

export async function readBytes(file) {
  return new Uint8Array(await file.arrayBuffer());
}

export async function fetchBytes(url) {
  const response = await fetch(url);
  if (!response.ok) throw new Error(`Could not load ${url} (HTTP ${response.status}).`);
  return new Uint8Array(await response.arrayBuffer());
}

export const isGzip = (bytes) => bytes[0] === GZIP_MAGIC[0] && bytes[1] === GZIP_MAGIC[1];

export async function gunzip(bytes) {
  const stream = new Blob([bytes]).stream().pipeThrough(new DecompressionStream("gzip"));
  return new Uint8Array(await new Response(stream).arrayBuffer());
}

export function parseJson(bytes, label) {
  try {
    return JSON.parse(new TextDecoder().decode(bytes));
  } catch (error) {
    throw new Error(`${label} is not valid JSON: ${error.message}`);
  }
}

/** Structure file: SHA-256 of the bytes as opened (like simprep), plain text, format. */
export async function readStructure(name, bytes) {
  const plain = isGzip(bytes) ? await gunzip(bytes) : bytes;
  const base = name.toLowerCase().replace(/\.gz$/, "");
  const format = /\.(pdb|ent)$/.test(base) ? "pdb" : "mmcif";
  return { name, sha256: await sha256Hex(bytes), text: new TextDecoder().decode(plain), format };
}

/** Offer a JSON document as a file download. */
export function download(document_, fileName) {
  const blob = new Blob([`${JSON.stringify(document_, null, 2)}\n`], { type: "application/json" });
  const url = URL.createObjectURL(blob);
  const link = Object.assign(document.createElement("a"), { href: url, download: fileName });
  document.body.append(link);
  link.click();
  link.remove();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}
