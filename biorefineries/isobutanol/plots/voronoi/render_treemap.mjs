#!/usr/bin/env node
// Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
// Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
// UIUC open-source license -- see the biosteam LICENSE.txt.
//
// Stage 2 of the proteome-allocation Voronoi-treemap figure: render the
// Stage-1 allocation JSON into a grid of amCharts5 VoronoiTreemap charts
// and export PNG/SVG/PDF via headless Chromium (Puppeteer).
//
//   node render_treemap.mjs --doc alloc.json --out /path/stem [--cols 3] [--scale 2]

import { readFileSync } from "node:fs";
import { fileURLToPath, pathToFileURL } from "node:url";
import { dirname, resolve } from "node:path";
import puppeteer from "puppeteer";

function arg(name, def) {
  const i = process.argv.indexOf("--" + name);
  return i >= 0 && i + 1 < process.argv.length ? process.argv[i + 1] : def;
}

const here = dirname(fileURLToPath(import.meta.url));
const docPath = resolve(arg("doc"));
const outStem = resolve(arg("out"));
const cols = parseInt(arg("cols", "0"), 10) || 0;
const scale = parseFloat(arg("scale", "2"));

if (!arg("doc") || !arg("out")) {
  console.error("usage: render_treemap.mjs --doc <json> --out <stem> "
    + "[--cols N] [--scale S]");
  process.exit(2);
}

const doc = JSON.parse(readFileSync(docPath, "utf-8"));
const htmlUrl = pathToFileURL(resolve(here, "treemap.html")).href;

const browser = await puppeteer.launch({
  headless: "new",
  args: ["--allow-file-access-from-files", "--no-sandbox"],
});
try {
  const page = await browser.newPage();
  await page.setViewport({ width: 1400, height: 1000, deviceScaleFactor: scale });
  page.on("console", (m) => { if (m.type() === "error") console.error("[page]", m.text()); });
  page.on("pageerror", (e) => console.error("[pageerror]", e.message));

  // inject the doc + column override BEFORE the page's own script runs
  await page.evaluateOnNewDocument((d, c) => {
    window.TREEMAP_DOC = d;
    if (c > 0) window.__COLS__ = c;
  }, doc, cols);

  await page.goto(htmlUrl, { waitUntil: "networkidle0" });
  await page.waitForFunction(() => window.__RENDER_DONE__ === true,
    { timeout: 60000 });
  // amCharts animates polygons in; give the tessellation a beat to settle
  await new Promise((r) => setTimeout(r, 800));

  const figure = await page.$("#figure");
  const box = await figure.boundingBox();

  // PNG (raster, high-DPI)
  await figure.screenshot({ path: outStem + ".png" });

  // PDF (vector page sized to the figure)
  await page.pdf({
    path: outStem + ".pdf",
    printBackground: true,
    width: Math.ceil(box.width + box.x * 2) + "px",
    height: Math.ceil(box.height + box.y * 2) + "px",
    pageRanges: "1",
  });

  // SVG (serialize amCharts' rendered SVG per chart, wrapped in one root)
  const svg = await page.evaluate(() => {
    const svgs = Array.from(document.querySelectorAll("#grid .chart svg"));
    if (!svgs.length) return null;
    const pad = 14, cellW = 300, titleH = 40;
    const per = svgs.map((s, i) => {
      const clone = s.cloneNode(true);
      clone.removeAttribute("style");
      const cols = getComputedStyle(document.getElementById("grid"))
        .gridTemplateColumn.split(" ").length;
      const r = Math.floor(i / cols), c = i % cols;
      const x = c * (cellW + pad), y = r * (cellW + titleH + pad);
      return '<g transform="translate(' + x + ',' + y + ')">'
        + clone.outerHTML + "</g>";
    }).join("");
    const cols = getComputedStyle(document.getElementById("grid"))
      .gridTemplateColumns.split(" ").length;
    const rows = Math.ceil(svgs.length / cols);
    const W = cols * (cellW + pad), H = rows * (cellW + titleH + pad);
    return '<svg xmlns="http://www.w3.org/2000/svg" width="' + W
      + '" height="' + H + '" viewBox="0 0 ' + W + " " + H + '">'
      + per + "</svg>";
  });
  if (svg) {
    const { writeFileSync } = await import("node:fs");
    writeFileSync(outStem + ".svg", svg, "utf-8");
  } else {
    console.error("warning: no chart SVGs found; skipped .svg export");
  }

  console.log("wrote " + outStem + ".png / .pdf" + (svg ? " / .svg" : ""));
} finally {
  await browser.close();
}
