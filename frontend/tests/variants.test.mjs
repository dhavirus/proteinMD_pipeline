import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import { gunzipSync } from "node:zlib";
import { fileURLToPath } from "node:url";
import { test } from "node:test";
import { residueNames, tokens } from "../lib/atom_site.js";
import { parseMutation, variantFromDraft, variantName } from "../lib/variants.js";

const panel = (name) => fileURLToPath(new URL(`../../tests/panel/${name}`, import.meta.url));
const fiveFql = () => residueNames(gunzipSync(readFileSync(panel("5FQL.cif.gz"))).toString(), "mmcif");

test("mmCIF atom_site gives author residue names (5FQL A:468 is ARG)", () => {
  const names = fiveFql();
  assert.equal(names.get("A:468"), "ARG");
  assert.equal(names.get("A:84"), "ALS");
  assert.equal(names.get("A:1551"), "CA");
});

test("PDB records give residue names too", () => {
  const pdb = "ATOM      1  N   ARG A 468      1.000   2.000   3.000  1.00 60.00           N\nEND\n";
  assert.deepEqual([...residueNames(pdb, "pdb")], [["A:468", "ARG"]]);
});

test("CIF tokens keep quoted values with spaces together", () => {
  assert.deepEqual(tokens(`HETATM 1 O "O5'" 'a b' x`), ["HETATM", "1", "O", "O5'", "a b", "x"]);
});

test("mutations parse from one- or three-letter codes", () => {
  assert.deepEqual(parseMutation("A:468 R>Q").mutation, { chain: "A", seq_num: 468, ins_code: "", from: "ARG", to: "GLN" });
  assert.equal(parseMutation("A:468 arg > trp").mutation.to, "TRP");
  assert.match(parseMutation("A468 R>Q").error, /not a mutation/);
  assert.match(parseMutation("A:468 R>R").error, /same residue/);
  assert.match(parseMutation("A:468 R>X").error, /standard amino-acid/);
});

test("a variant is checked against the loaded structure and named conventionally", () => {
  const names = fiveFql();
  const ok = variantFromDraft({ mutations: "A:468 R>Q", rationale: "study" }, names, []);
  assert.deepEqual(ok.errors, []);
  assert.equal(ok.variant.name, "R468Q");
  const wrong = variantFromDraft({ mutations: "A:468 K>Q, A:9999 R>W", rationale: "" }, names, []);
  assert.equal(wrong.errors.length, 3);
  assert.match(wrong.errors[0], /is ARG in the loaded structure, not LYS/);
  assert.match(variantFromDraft({ mutations: "A:468 R>Q", rationale: "x" }, names, ["R468Q"]).errors[0], /already exists/);
  assert.equal(variantName([ok.variant.mutations[0], { ...ok.variant.mutations[0], seq_num: 470, to: "TRP" }]), "R468Q_R470W");
});
