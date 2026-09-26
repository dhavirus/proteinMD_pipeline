# ADR-0003: Front-end architecture

- Status: accepted (TASK-002)
- Date: 2026-09-26

## Context

TASK-002 adds a static review page (`frontend/`, GitHub Pages) where a user opens a
structure, its `findings.json` and optionally a manifest, records a decision per finding,
edits regions of interest and exports `manifest.json`. CLAUDE.md requires plain HTML/JS
without a build step, Mol* and ajv from a CDN at pinned exact versions, and a front end
that records decisions but never transforms coordinates. ADR-0002 decided that the page
consumes `findings.json` and does not re-run detection or severity weighting.

## Options

1. **Framework + bundler** (React/Vite or similar). Rich components, but a build step,
   a `node_modules` tree and a toolchain to keep current; CLAUDE.md requires an ADR to
   justify that, and nothing in v0.1 needs it.
2. **One large inline script.** No tooling, but logic and DOM code tangle, and nothing can
   be unit-tested without a browser.
3. **Native ES modules, split into pure logic and a thin UI layer.** No build step; the
   browser loads `app.js` as a module. The logic modules have no DOM, Mol* or ajv
   dependency and are tested with Node's built-in test runner.

## Decision

Option 3.

- `frontend/lib/`: pure logic, no DOM, no network, no dependencies:
  `canonical.js` (RFC 8785 JSON + SHA-256), `residues.js`, `regions.js`, `findings.js`,
  `decisions.js`, `state.js`, `manifest.js`. Tested by `node --test "tests/*.test.mjs"`.
- `frontend/ui/`: the edges: `io.js` (file reading, gunzip, download), `validate.js`
  (ajv), `view3d.js` (Mol* through MolViewSpec), `draft.js` (localStorage).
- `frontend/app.js`: wires them to `index.html`.
- **State**: one plain object (`lib/state.js`) replaced on every change, never mutated;
  every change re-renders the page from it.
- **I/O**: files are opened with file inputs and read in the browser; nothing is uploaded
  and there is no server. The structure is hashed (SHA-256 over the bytes as opened, as
  `simprep` hashes them) and refused if it differs from `findings.json` `input.sha256`.
  The export is a download (plus copy to clipboard).
- **Validation**: ajv validates `findings.json`, a loaded manifest, and the manifest before
  export, using the repository's `schema/*.schema.json` files (served beside the page, not
  copied). ajv runs with `strict: false, validateFormats: false`, matching the Python
  validator, which does not assert formats. Parity was checked on every schema example
  and rule file (19 documents, 0 disagreements).
- **Snapshot hash**: `findings_sha256` is computed over RFC 8785 canonical JSON in both
  languages (`simprep.canonical`, `frontend/lib/canonical.js`), because Python's
  `json.dumps` and JavaScript's `JSON.stringify` format numbers differently (`1.0` vs
  `1`). `simprep` verifies the hash when it loads a manifest.
- **Drafts**: unsent decisions and regions are kept in `localStorage`, keyed by the input
  SHA-256, labelled on the page as "kept in this browser only"; the exported manifest is
  the record.
- **Severity after region edits** is not recomputed in the browser (ADR-0002; TASK-002
  decision 3): the page marks severities as stale until `simprep audit --manifest` is
  re-run.
- **Pinned CDN builds with Subresource Integrity**:

  | library | version | file | integrity |
  |---|---|---|---|
  | Mol* | 5.11.0 | `molstar@5.11.0/build/viewer/molstar.js` | `sha384-5Mfx4eL50NkWPky+mcH//qY0sbml4il0CLFFmrMp8uv/saB3Z6uZMHn2dUpAnH92` |
  | Mol* | 5.11.0 | `molstar@5.11.0/build/viewer/molstar.css` | `sha384-RIontCdJN53gEl2fmiHN+4bscIBvaUaOiCeeGktXqmFqdEBF+COnSdt9O4IKFSvq` |
  | ajv (2020 dialect) | 8.17.1 | `ajv-dist@8.17.1/dist/ajv2020.min.js` | `sha384-vsMiC2R9fcX/SqeA9xJ+Lpbbe1fRi7RDpgHenZH7A8yXAoX+9T4m8MfuwrafZAAL` |

  Hashes were computed from the npm tarballs of those versions (jsDelivr serves the npm
  files unchanged). `ajv` itself ships no browser bundle; `ajv-dist` is its UMD build
  (global `ajv2020`). To bump a version: download the new tarball, recompute the SHA-384,
  update `index.html` and this table, and run the browser smoke test.
- **Hosting**: GitHub Pages from `main` (`.github/workflows/pages.yml`). The site keeps
  the repository layout (`frontend/`, `schema/`, `tests/panel/5FQL.cif.gz`) so the page's
  relative paths are the same locally and deployed.

## Consequences

- No build step and no runtime dependencies to install; a contributor can serve the
  repository root with any static server and open `/frontend/`.
- Browser behaviour (Mol*, ajv, downloads) is covered by a Playwright smoke test that needs
  the CDN and therefore is not in CI (`frontend/tests/smoke/smoke.mjs`; it can map the CDN
  to local npm copies with `SIMPREP_CDN_DIR`). The logic, and the manifests the logic
  produces, are covered offline in CI (`node --test`, `tests/test_frontend_contract.py`).
- Offline use of the deployed page is not possible: Mol* and ajv come from the CDN.
- Changing the snapshot hash to RFC 8785 changed `findings_sha256` for manifests written
  by TASK-001's `simprep`; such manifests must be re-audited once (none exist outside
  tests).
