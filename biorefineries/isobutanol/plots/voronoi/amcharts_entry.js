// Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
// Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
// UIUC open-source license -- see the biosteam LICENSE.txt.
//
// Bundle entry for build_amcharts.mjs. amCharts5 ships ES modules (for
// bundlers), not browser UMD/script bundles, so a classic <script> tag on
// node_modules/@amcharts/amcharts5/index.js cannot expose the `am5` global.
// esbuild bundles this entry into an IIFE (vendor/am5.bundle.js) that
// attaches the amCharts core and hierarchy namespaces to `window`, which is
// what treemap.html loads.
import * as am5 from "@amcharts/amcharts5";
import * as am5hierarchy from "@amcharts/amcharts5/hierarchy";

window.am5 = am5;
window.am5hierarchy = am5hierarchy;
