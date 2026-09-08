#!/usr/bin/env node
// Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
// Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
// UIUC open-source license -- see the biosteam LICENSE.txt.
//
// Build step for the Voronoi-treemap renderer: bundle amCharts5 (core +
// hierarchy) from node_modules into a single browser IIFE that exposes the
// `am5` and `am5hierarchy` globals, so treemap.html can load amCharts with a
// plain <script> tag offline. amCharts5 distributes ES modules for bundlers
// only -- no UMD/script build -- so this esbuild pass is required.
//
// Runs automatically on `npm install` (postinstall) and via `npm run build`.
// Output (vendor/am5.bundle.js) is gitignored and regenerated from the pinned
// package-lock, never committed.
import { build } from "esbuild";
import { mkdirSync } from "node:fs";
import { dirname, resolve } from "node:path";
import { fileURLToPath } from "node:url";

const here = dirname(fileURLToPath(import.meta.url));
const entry = resolve(here, "amcharts_entry.js");
const outfile = resolve(here, "vendor", "am5.bundle.js");

mkdirSync(dirname(outfile), { recursive: true });
await build({
  entryPoints: [entry],
  bundle: true,
  format: "iife",
  outfile,
  legalComments: "none",
  logLevel: "info",
});
console.log("built " + outfile);
