import assert from "node:assert/strict";
import { test } from "node:test";
import { formatResidueList, parseResidueList } from "../lib/residues.js";
import { regionFromDraft, sameRegions } from "../lib/regions.js";

test("parses single residues, insertion codes, negatives and ranges", () => {
  const { residues, errors } = parseResidueList("A:45, A:52A;B:-3  A:444-446");
  assert.deepEqual(errors, []);
  assert.equal(formatResidueList(residues), "A:45, A:52A, B:-3, A:444, A:445, A:446");
});

test("negative ranges and duplicates", () => {
  const { residues } = parseResidueList("A:-2--1 A:-1");
  assert.equal(formatResidueList(residues), "A:-2, A:-1");
});

test("reports malformed tokens and bad ranges", () => {
  const { errors } = parseResidueList("45 A:10-5 A:3A-5 A:1-5000");
  assert.equal(errors.length, 4);
});

test("region drafts are validated", () => {
  const ok = regionFromDraft({ name: "site", description: "Ca site", residues: "A:45 A:46" }, []);
  assert.deepEqual(ok.errors, []);
  assert.equal(ok.region.residues.length, 2);
  const bad = regionFromDraft({ name: "bad name", description: "", residues: "" }, []);
  assert.equal(bad.region, null);
  assert.equal(bad.errors.length, 3);
  const taken = regionFromDraft({ name: "site", description: "x", residues: "A:1" }, ["site"]);
  assert.match(taken.errors[0], /already used/);
});

test("region comparison ignores order and descriptions", () => {
  const a = [{ name: "x", description: "one", residues: [{ chain: "A", seq_num: 2, ins_code: "" }, { chain: "A", seq_num: 1, ins_code: "" }] }];
  const b = [{ name: "x", description: "two", residues: [{ chain: "A", seq_num: 1, ins_code: "" }, { chain: "A", seq_num: 2, ins_code: "" }] }];
  assert.ok(sameRegions(a, b));
  assert.ok(!sameRegions(a, []));
});
