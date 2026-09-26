// Browser smoke test for the review page (TASK-002). Needs Playwright and, normally, the
// network (Mol* and ajv load from the pinned CDN), so it is not part of CI.
//
//   node frontend/tests/smoke/smoke.mjs                # uses the CDN
//   SIMPREP_CDN_DIR=/path node frontend/tests/smoke/smoke.mjs
//
// With SIMPREP_CDN_DIR set, requests to https://cdn.jsdelivr.net/npm/<pkg>@<ver>/<file> are
// answered from <dir>/<pkg>@<ver>/<file> (unpacked npm tarballs of the same versions);
// Subresource Integrity still checks that the bytes match the pinned hashes.
// Writes the exported manifest to $SIMPREP_SMOKE_OUT (default: a temp dir) and exits 1 on
// any failed check.
import { execFileSync } from "node:child_process";
import { createReadStream, existsSync, mkdtempSync, statSync } from "node:fs";
import { createServer } from "node:http";
import { tmpdir } from "node:os";
import { extname, join, normalize } from "node:path";
import { fileURLToPath } from "node:url";

const playwrightPath = process.env.PLAYWRIGHT_MODULE || "playwright";
const { chromium } = await import(playwrightPath);
const ROOT = fileURLToPath(new URL("../../../", import.meta.url));
const CDN = "https://cdn.jsdelivr.net/npm/";
const TYPES = { ".html": "text/html", ".js": "text/javascript", ".mjs": "text/javascript", ".css": "text/css", ".json": "application/json", ".gz": "application/gzip" };
const results = [];
const check = (name, ok, detail = "") => {
  results.push({ name, ok });
  console.log(`${ok ? "PASS" : "FAIL"}  ${name}${detail ? `  (${detail})` : ""}`);
};

function serve(root) {
  const server = createServer((request, response) => {
    const path = normalize(join(root, decodeURIComponent(new URL(request.url, "http://x").pathname)));
    const file = existsSync(path) && statSync(path).isDirectory() ? join(path, "index.html") : path;
    if (!file.startsWith(root) || !existsSync(file)) {
      response.writeHead(404).end();
      return;
    }
    response.writeHead(200, { "content-type": TYPES[extname(file)] || "application/octet-stream" });
    createReadStream(file).pipe(response);
  });
  return new Promise((resolve) => server.listen(0, "127.0.0.1", () => resolve(server)));
}

async function routeCdn(page, dir) {
  await page.route(`${CDN}**`, (route) => {
    const file = join(dir, route.request().url().slice(CDN.length));
    route.fulfill(existsSync(file) ? { path: file, headers: { "access-control-allow-origin": "*" } } : { status: 404 });
  });
  await page.route("https://fonts.googleapis.com/**", (route) => route.abort());
}

async function decideAll(page) {
  const rows = page.locator(".row");
  const count = await rows.count();
  for (let index = 0; index < count; index += 1) {
    await rows.nth(index).click();
    // Keep the recommended option when it is final, otherwise take the first final one.
    if (!(await page.locator('input[name="option"]:checked[data-final="true"]').count())) {
      await page.locator('input[name="option"][data-final="true"]').first().check();
    }
    await page.fill("#decision-rationale", "smoke test: recommended or first final option");
    await page.fill("#decision-by", "smoke-test");
    await page.click('#decision-form button[type="submit"]');
    await page.waitForSelector("#decision-status.ok");
  }
  return count;
}

async function run(page, base, outDir) {
  const consoleErrors = [];
  page.on("pageerror", (error) => consoleErrors.push(error.message));
  await page.goto(`${base}/frontend/`);
  await page.click("#load-demo");
  await page.waitForSelector("#review:not([hidden])", { timeout: 30000 });
  check("demo shows 16 findings", (await page.locator(".row").count()) === 16);

  await page.locator(".row", { hasText: "metals/A:1551" }).click();
  const status = await page.waitForFunction(
    () => /Showing|failed|did not load/.test(document.getElementById("viewer-status").textContent), null, { timeout: 60000 },
  ).then(() => page.textContent("#viewer-status"));
  check("Mol* focuses the Ca2+ site", status.includes("Showing metals/A:1551"), status);

  await page.locator(".row", { hasText: "nonstandard_residues/A:84" }).click();
  await page.check("#option-revert_to_parent");
  await page.fill("#decision-rationale", "smoke test");
  await page.fill("#decision-by", "smoke-test");
  await page.click('#decision-form button[type="submit"]');
  check("explicit-choice option needs confirmation", /deliberately/.test(await page.textContent("#decision-status")));

  const decided = await decideAll(page);
  check("every finding decided", (await page.textContent("#summary")).includes(`${decided} of ${decided} decided`));

  await page.fill("#region-name", "active_site");
  await page.fill("#region-description", "Ca2+ site and catalytic FGly");
  await page.fill("#region-residues", "A:45, A:46, A:84, A:334-335, A:1551");
  await page.click('#region-form button[type="submit"]');
  check("region edit marks severities stale", await page.locator(".notice.stale").isVisible());

  const [download] = await Promise.all([page.waitForEvent("download"), page.click("#export-build")]);
  const manifestPath = join(outDir, "manifest.json");
  await download.saveAs(manifestPath);
  check("manifest exported", existsSync(manifestPath));
  check("no page errors", consoleErrors.length === 0, consoleErrors.join("; "));
  return manifestPath;
}

async function checkManifestWithDraft(page, base, manifestPath) {
  // The browser still holds a draft from run(); opening the exported manifest must show
  // the manifest and offer the draft, not apply it.
  await page.goto(`${base}/frontend/`);
  await page.setInputFiles("#file-structure", join(ROOT, "tests/panel/5FQL.cif.gz"));
  await page.setInputFiles("#file-findings", join(ROOT, "frontend/demo/5FQL.findings.json"));
  await page.setInputFiles("#file-manifest", manifestPath);
  await page.click('#open-form button[type="submit"]');
  await page.waitForSelector("#review:not([hidden])", { timeout: 30000 });
  const offered = page.locator(".notice", { hasText: "unsaved work" });
  check("opened manifest is shown and the browser draft only offered", (await offered.isVisible())
    && (await page.textContent("#summary")).includes("16 of 16 decided"));
  check("expert review is shown as not final", await expertReviewIsNotFinal(page));
  await offered.getByRole("button", { name: "Discard it" }).click();
  check("discarding the draft removes the offer", !(await page.locator(".notice", { hasText: "unsaved work" }).count()));
}

async function expertReviewIsNotFinal(page) {
  await page.locator(".row", { hasText: "unrecognized/group/A:1567" }).click();
  await page.check("#option-expert_review");
  const note = await page.locator("#decision-extra").textContent();
  await page.fill("#decision-rationale", "smoke test: still open");
  await page.click('#decision-form button[type="submit"]');
  const chip = await page.locator(".row", { hasText: "unrecognized/group/A:1567" }).textContent();
  return note.includes("stays undecided") && chip.includes("not final")
    && (await page.textContent("#summary")).includes("15 of 16 decided");
}

async function checkMismatch(page, base) {
  await page.goto(`${base}/frontend/`);
  await page.setInputFiles("#file-structure", join(ROOT, "tests/panel/3KS3.cif.gz"));
  await page.setInputFiles("#file-findings", join(ROOT, "frontend/demo/5FQL.findings.json"));
  await page.click('#open-form button[type="submit"]');
  await page.waitForSelector("#open-errors li", { timeout: 30000 });
  const text = await page.textContent("#open-errors");
  check("wrong structure is refused with both hashes", (text.match(/SHA-256 [0-9a-f]{64}/g) || []).length === 2);
}

function checkWithSimprep(manifestPath, outDir) {
  const simprep = (...args) => {
    try {
      execFileSync("simprep", args, { cwd: ROOT, stdio: "pipe" });
      return 0;
    } catch (error) {
      return error.status;
    }
  };
  check("simprep manifest status exits 0", simprep("manifest", "status", manifestPath, "--structure", "tests/panel/5FQL.cif.gz") === 0);
  check("re-audit with the exported manifest succeeds",
    simprep("audit", "tests/panel/5FQL.cif.gz", "--manifest", manifestPath, "--out", join(outDir, "reaudit")) === 0);
}

const outDir = process.env.SIMPREP_SMOKE_OUT || mkdtempSync(join(tmpdir(), "simprep-smoke-"));
const server = await serve(ROOT);
const base = `http://127.0.0.1:${server.address().port}`;
const browser = await chromium.launch(process.env.CHROMIUM_PATH ? { executablePath: process.env.CHROMIUM_PATH } : {});
try {
  const context = await browser.newContext({ acceptDownloads: true });
  const page = await context.newPage();
  if (process.env.SIMPREP_CDN_DIR) await routeCdn(page, process.env.SIMPREP_CDN_DIR);
  const manifestPath = await run(page, base, outDir);
  await checkManifestWithDraft(page, base, manifestPath);
  await checkMismatch(page, base);
  checkWithSimprep(manifestPath, outDir);
} finally {
  await browser.close();
  server.close();
}
const failed = results.filter((r) => !r.ok).length;
console.log(`${results.length - failed}/${results.length} checks passed; output in ${outDir}`);
process.exit(failed ? 1 : 0);
