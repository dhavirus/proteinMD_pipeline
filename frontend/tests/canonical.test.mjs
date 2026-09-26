import assert from "node:assert/strict";
import { test } from "node:test";
import { canonicalJson, sha256Canonical, sha256Hex } from "../lib/canonical.js";
import { exampleManifest } from "./fixtures.mjs";

test("keys are sorted and whitespace removed", () => {
  assert.equal(canonicalJson({ b: 1, a: [true, null, "x"] }), '{"a":[true,null,"x"],"b":1}');
});

test("numbers follow ECMAScript formatting", () => {
  assert.equal(canonicalJson([1.0, 0.29, 1e21, 1e-7, -0]), "[1,0.29,1e+21,1e-7,0]");
});

test("non-finite numbers are rejected", () => {
  assert.throws(() => canonicalJson(Number.NaN), TypeError);
});

test("sha256 of known text", async () => {
  assert.equal(await sha256Hex("abc"), "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
});

test("snapshot hash matches the one Python recorded", async () => {
  const snapshot = exampleManifest().findings_snapshot;
  assert.equal(await sha256Canonical(snapshot.findings), snapshot.findings_sha256);
});
