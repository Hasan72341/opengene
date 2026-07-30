import assert from "node:assert/strict";
import { readFile } from "node:fs/promises";
import test from "node:test";

const root = new URL("../../", import.meta.url);

test("browser workstation is wired as a local module application", async () => {
  const [html, app, serviceWorker] = await Promise.all([
    readFile(new URL("index.html", root), "utf8"),
    readFile(new URL("src/app.js", root), "utf8"),
    readFile(new URL("service-worker.js", root), "utf8"),
  ]);
  assert.match(html, /src="src\/app\.js"/);
  assert.match(app, /function renderVerificationTool/);
  assert.match(app, /function renderAnnotationTool/);
  assert.match(app, /function createAssemblyResult/);
  assert.match(serviceWorker, /opengene-studio-v3/);
});
