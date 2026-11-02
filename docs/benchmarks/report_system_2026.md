# Report system prototype

Prototype run 2026-11-01 to 2026-11-02 (roadmap Phase 11c, items P1 and P2), on dev `c8115c2`, the first commit with the Phase 11b summary tables. It answers two questions before Phase 11d rebuilds the summary's documents and figures:

- **(P1) Report engine.** The same two sections, an overview (run, samples, read QC, community) and an ARG section (read-level method comparison, SARG, contig ARGs), plus a methods section built from `software_versions.yml` with references, were built two ways: **(a)** Quarto with its Typst PDF backend, one `.qmd` per section giving a multi-page HTML site and a PDF; **(b)** Typst templates reading the data as JSON for the PDF, plus a Python/Jinja multi-page HTML site.
- **(P2) Figures.** One shared, paper-ready style applied to three figures (the taxonomy heatmap, the QC read-fate plot and the SARG heatmap) at 4 and at 120 samples, drawn two ways: **(i)** matplotlib for the static figures plus Vega-Lite (Altair) for the interactive ones; **(ii)** one Vega-Lite spec per figure for both, with static export through vl-convert (no browser).

Nothing is wired into the pipeline. Only this document is committed; the prototype code and renders stay outside the repository.

**Decision (owner):** pending.

**Outcome (recommendations):**

- **(P1) Engine (b): Typst for the PDF, Python/Jinja for the HTML site.** Its PDF was correct on the first render: long tables break across pages with the header repeated, numbers are right-aligned, wide tables go on landscape pages, figures and tables are numbered and cross-referenced, and the bibliography uses Typst's built-in Nature style. A whole report took 4-7 s (Typst itself 0.7-1.9 s). Quarto (a) needed four workarounds before its PDF was usable, and a fifth problem was still open at the end:
  1. Tables are wrapped in figures of kind `quarto-float-tbl` that do not break across pages, so the rows past the page end were printed on top of each other. The fix is a raw Typst rule, `#show figure: set block(breakable: true)`.
  2. The book format's bundled `orange-book` template ignores `mainfont` (body text in Libertinus Serif), inserts blank pages before chapters, and fails to compile without an `author`.
  3. pandoc gives the table columns widths in proportion to the text, so short columns wrap into four lines and the wide SARG table overprints even on a landscape page. This one was still open.
  4. The conda package works only with its activation variables set (`QUARTO_*`, one of them pointing at a path from the package's own build). Nextflow does not source activation scripts, so the module would have to set them.
  5. Under Apptainer, Quarto needs `XDG_CACHE_HOME`, `XDG_DATA_HOME` and `DENO_DIR` in the task directory, or it stops with `Read-only file system: mkdir '~/.cache/quarto'`.

  Quarto's HTML site is also not offline as configured. Its theme `@import`s Google Fonts. Opened from disk (`file://`), the browser blocks its module scripts (`quarto.js`, `tabsets.js`), so search and the page tools stop working. It took 39-49 s per run (188-219 s at 120 samples), and its image is 645 MB larger than today's report image. Engine (b)'s site works from `file://` with no failed request. The cost of (b) is maintenance: each section exists twice, as a Jinja template and a Typst file (about 20 lines each here). The HTML reference list and the figure and table numbers come from Python, not from a citation processor. Phase 11d can limit the duplication by keeping every number, caption and table in the shared Python layer, as the prototype does, so that only the short prose is written twice.
- **Typst's own HTML export is not usable for (b)'s HTML yet.** Typst 0.15.0 (2026-06-15) added an experimental "bundle" target for multi-file sites, and "the rest of HTML export remains feature-flagged". Typst 0.15.1 (2026-07-17) changed nothing in HTML export. The HTML reference page still says it is "very incomplete and only available for experimentation behind a feature flag. Do not use this feature for production use cases", and that Typst "does not output CSS style sheets". Compiling engine (b)'s sources with `--features html` gave one 2.4 MB page that ignored every page, alignment and page-break rule.
- **(P2) (i): matplotlib for the canonical static figures, Vega-Lite (Altair) for the interactive versions of the large figures only.** Both backends met the basic style: Inter embedded in PDF and SVG, the Okabe-Ito and viridis palettes, PNG at 600 dpi, a source CSV beside each figure, and widths of exactly 89 and 183 mm. The static paper needs are what separates them:
  - Vega-Lite has no fill patterns. A low-depth heatmap column can be striped with one diagonal rule per cell, but a bar cannot be hatched (grey only).
  - Vega-Lite sizes the plot area, not the figure. An exact width took a measure-and-refit loop of 2 to 10 renders, and the height is not controlled: at 120 samples the heatmaps came out 252-260 mm tall, more than the 247 mm page.
  - Legend placement is coarser. A legend beside a strip widens that sub-chart, and a facet of a concatenated chart is not allowed.
  - "One spec for both" still needed two parameter sets. At print sizes (6-7 pt type, 89 mm) the interactive chart is unreadable on screen, so the HTML version uses the same spec function with a 1.8× type and cell scale.

  With (i), matplotlib keeps exact sizes, hatching and legend control for the paper figures, and Altair adds hover, group filtering and zoom only where the HTML gains from them. Phase 11d's interactivity scope is the large figures only, so the second implementation stays small.
- **Plotly is out:** its static export (kaleido 1.5.0) raised `RuntimeError` asking for Google Chrome, in an image with no browser.
- **Image for 11d** (if (b) and (i) are chosen): the current report stack plus `typst=0.15.1`, `jinja2=3.1.6`, `altair=6.3.0`, `vl-convert-python=1.9.0.post1`, `font-inter=4.1` and `pyyaml`. It is 1.12 GB as a docker image and 428 MB as a sif, against 895 MB and 338 MB today. Its Seqera docker image and sif twin are built, and both ran every render here.

## Summary

| Question | (a) Quarto 1.10.19 (Typst 0.15.1) | (b) Typst 0.15.1 + Python/Jinja | Recommendation |
| --- | --- | --- | --- |
| PDF out of the box | book: serif body (template ignores `mainfont`), blank pages, tables overprinted at page ends, empty "References" plus a "Bibliography" heading; plain `typst` format: fixed font, same table problem | correct on the first render | (b) |
| PDF after fixes | tables break (raw Typst `breakable` rule); column widths and the wide landscape table still wrong | 12 / 16 pages (4 / 120 samples), 0.5 / 1.7 MB | (b) |
| Render time, whole report (4 / 120 samples) | HTML 39-49 / 188-219 s; PDF 19-36 / 52-79 s | 4-6 / 6-8 s (Typst 0.7 / 1.8 s) | (b) |
| Image (docker / sif) | 1.65 GB / 578 MB (engine alone +645 MB) | 1.12 GB / 428 MB (engine alone +115 MB) | (b) |
| Apptainer (read-only, no network, no home) | works with `QUARTO_*` activation variables and `XDG_*`/`DENO_DIR` in the task directory (2 MB written there) | works unchanged | (b) |
| HTML offline | Google Fonts `@import`; module scripts blocked from `file://` | no external request, works from `file://` | (b) |
| Numbering, cross-references | native (`@fig-`, `@tbl-`, chapter-prefixed in the book) | Typst labels in the PDF; Python counters in the HTML (same numbers) | either |
| Methods and references | citeproc with a shipped `nature.csl` | Typst `bibliography()` (built-in `nature` style); hand-formatted list in the HTML | either |
| Adding a section | one `.qmd`; each table and figure written twice behind `content-visible when-format` | a Typst file and a Jinja template; data and captions once, in Python | (a) on paper, (b) in practice |
| Static figures (P2) | (i) matplotlib: exact size, hatching, legends placed freely | (ii) vl-convert: exact width only after refit, height not controlled, no bar hatching | (i) |
| Interactive figures (P2) | (i) and (ii) use the same Vega-Lite: hover, group filter, zoom on quantitative axes, 0.9 MB library once per site | | Vega-Lite |
| Plotly | static export needs Chrome | | out |

The rest of this document gives the setup and the measurements behind each row.

## Setup

**Data.** The summary tables of the Phase 11a asm reference results were regenerated with dev `c8115c2`'s `results_summary.py` in the pinned report image. That run has all 14 stages on and every ARG branch (K-MARVEL, RGI, ARGs-OAP, contig DeepARG, geNomad). One shared preparation step fed every prototype:

- **"4 samples"**: the real tables with the samples renamed A-D. The asm run has 4 samples, not 3, and none below the 1 M clean-pair cutoff, so the prototype raised the cutoff to 3 M to put two samples (2.9 M and 1.4 M pairs) in low depth and exercise the hatching. This run was started with `--quality_control false`, so its read-fate bars are 100 % clean.
- **"120 samples"**: synthetic, in four groups of 30. Each group is based on one real sample, with log-normal noise and 15 % dropouts. 10 samples are below 1 M clean pairs, and the QC and decontamination fractions are synthetic. The long tables are copies of the real rows: 4,916 method-comparison rows, 9,715 contig ARGs and 4,942 contig-host rows.

The methods section reads the run's own `software_versions.yml` (140 lines, 51 tools and databases). Thirteen references were checked against PubMed: Nextflow, sylph, CARD 2023, DeepARG, geNomad, MEGAHIT, DIAMOND, GTDB-Tk v2, CheckM2, MicrobeCensus, SAMtools, ARGs-OAP v2.0 and Metabuli. The ARGs-OAP 4.0 the pipeline runs has no PubMed record of its own, so the v2.0 paper stood in. Phase 11d needs a curated tool → reference table.

**Both engines got the same inputs.** A shared Python layer (140 lines) selects the columns, formats the numbers, caps every PDF table at 60 rows with a note pointing to the TSV and the HTML, and turns `software_versions.yml` into a tools table with reference keys. Both HTML sites use the same 27-line table script (click a header to sort, type to filter; no library) and the same Vega-Lite specs with the Vega libraries bundled. Both PDFs embed the same matplotlib SVGs.

**Images** (Seqera Containers, `freeze`, conda-forge; checked inside):

| Image | Docker | sif | Contents (versions found inside) |
| --- | --- | --- | --- |
| current report image | 895 MB | 338 MB | Python 3.13.18, pandas 3.1.0, numpy 2.5.4, matplotlib 3.11.2, seaborn 0.13.2, h5py 3.16.0, openpyxl 3.1.5, biopython 1.88 |
| (b) + figures | 1.12 GB | 428 MB | + typst 0.15.1, jinja2 3.1.6, altair 6.3.0, vl-convert-python 1.9.0.post1 (Vega-Lite 6.4), font-inter 4.1, pyyaml 6.0.3 |
| (b) engine only | 1.01 GB | | + typst, jinja2, font-inter, pyyaml |
| figures only | 1.07 GB | | + altair, vl-convert-python, font-inter, pyyaml |
| (a) + figures | 1.65 GB | 578 MB | + quarto 1.10.19 (pandoc 3.10, deno 2.7.14, dart-sass 1.101.0, esbuild 0.28.1, typst 0.15.1 from conda), nbclient 0.11.0, nbformat 5.11.1, ipykernel 7.4.0, jupyter_core 5.9.1, tabulate 0.10.0, altair, vl-convert-python, font-inter, pyyaml |
| (a) engine only | 1.54 GB | | + quarto, the Jupyter packages, tabulate, font-inter, pyyaml |

Both sif twins (Amazon Linux 2023) have every command the renders call (`typst`, `quarto`, `deno`, `pandoc`, `cp`, `mkdir`, `mv`, `rm`, `date`, `fc-list`, and the Python libraries), so either could be pinned for singularity/apptainer under the CONTRIBUTING rule. Quarto 1.10.19 (2026-10-06) was used rather than 1.10.20, which came out two days before the run. Its conda package uses conda-forge's `typst`, the same 0.15.1 as image (b). The vl-convert pin is the 1.9.0 stable line, not the 2.0 release candidates. The renders ran in docker with `--network none`, and again under Apptainer 1.5.3 with `--containall --no-home --cleanenv --net --network none`. PDFs were checked with poppler 26.10.0 (`pdfinfo`, `pdffonts`), and the HTML was screenshotted from `file://` with Playwright 1.63.0 (Chromium) and no network. None of these tools went into the report image.

## P1: report engine

### Output quality

| | (a) book (`orange-book`) | (a2) plain `typst` format | (b) Typst templates |
| --- | --- | --- | --- |
| Pages, 4 / 120 samples | 21 / 29 | 15 / 21 | 12 / 16 |
| PDF size | 480 KB / 1.6 MB | 452 KB / 1.6 MB | 489 KB / 1.7 MB |
| Body font | Libertinus Serif (template); Inter only in tables and figures | Inter | Inter |
| Fonts embedded (`pdffonts`) | yes, all | yes, all | yes, all (static Inter TTFs given with `--font-path`) |
| Long tables | overprinted at the page end until `#show figure: set block(breakable: true)`; then break with a repeated header | same | break with a repeated header (`table.header`), zebra rows |
| Column widths | from pandoc's text widths: short columns wrap ("K-MAR-VEL + RGI" over four lines) | same | automatic; numbers right-aligned |
| Wide tables (12 genera × samples, 21 SARG types × samples) | `.landscape` ignored in the book | landscape pages work; the SARG table still overprints | landscape pages (`page(flipped: true)`); long SARG names wrap at `_` |
| Captions, numbering | chapter-prefixed ("Table 2.2"), caption above tables | "Table 5", above | "Table 5", above; figures "Figure 2" below |
| Cross-references | `@fig-`, `@tbl-` resolve | resolve | Typst labels resolve; the HTML uses the same numbers from Python |
| References | citeproc, Nature CSL (shipped with the project); an empty "References" heading plus a separate "Bibliography" | same | `bibliography(style: "nature")`, numbered in citation order |
| Blank pages | yes (chapters start on odd pages) | no | no |

Engine (b)'s PDF has a contents page, running headers, page numbers ("3 / 12") and section numbering, all from 42 lines of `lib.typ`. The matplotlib SVGs keep their text as text, and Typst renders it with the bundled Inter, so the figures match the body type.

In the HTML, Quarto gives more page furniture for free: a sidebar, search and a per-page table of contents. That furniture needs its own scripts, and those fail from `file://`. Its long tables are not in a scroll box, so the ARG page at 120 samples is 439,028 px tall. Engine (b)'s pages are plain by design: a sticky top navigation and a fact row, with tables in scroll boxes (520 px) with a filter box and sortable headers. Every figure has its interactive version followed by links to its SVG, PDF, 600 dpi PNG and source CSV.

### Render time and size

Wall time in docker (`--network none`) on a shared 80-core host, two or three runs each. The figures were rendered beforehand and are not included.

| | 4 samples | 120 samples |
| --- | --- | --- |
| (b) whole report (Typst PDF + Jinja site) | 3.6-5.8 s (Typst 0.7 s) | 5.6-8.2 s (Typst 1.7-2.5 s) |
| (a) HTML site | 39-56 s | 188-219 s |
| (a) PDF, book | 28-36 s | 58-79 s |
| (a) PDF, plain | 19-24 s | 52-69 s |
| Site size, (a) / (b) | 5.2 / 3.8 MB | 13.0 / 8.7 MB |
| Largest page, (a) / (b) | 174 / 63 KB | 2.3 / 1.1 MB |

Quarto starts a Jupyter kernel per chapter and runs pandoc over 5,000-row markdown tables, so its time grows with the table sizes. Both sites carry the same assets: the figures folder (2.4 MB), the Vega bundle (880 KB, once per site) and three WOFF2 faces of Inter (430 KB). Quarto adds 1.2 MB of `site_libs`.

### Apptainer

Both images ran read-only with no network, no home directory and a clean environment:

- **(b)** needed nothing. The only message was fontconfig's "No writable cache directories", which goes away with `XDG_CACHE_HOME` in the task directory. Typst uses no packages here (`@preview` imports would download at compile time), and its fonts come from `--font-path` with `--ignore-system-fonts`.
- **(a)** needed `XDG_CACHE_HOME`, `XDG_DATA_HOME` and `DENO_DIR` in the task directory, where it wrote 2 MB. Without them it stopped with `Read-only file system (os error 30): mkdir '/home/…/.cache/quarto'`. It also needed the conda activation variables (`QUARTO_DENO`, `QUARTO_PANDOC`, `QUARTO_TYPST`, `QUARTO_SHARE_PATH`, …), because a container started by Nextflow does not source `etc/conda/activate.d/quarto.sh`. That script's `QUARTO_DENO_DOM` points into the package's build directory and has to be overridden with `$CONDA_PREFIX/lib/deno_dom.so`. The `orange-book` Typst package used by the book format is bundled with Quarto, so nothing was downloaded.

Docker and Apptainer gave byte-identical Vega-Lite specs, interactive HTML, Vega-Lite SVGs and source CSVs. Some PNGs differed only by sub-pixel antialiasing (largest channel difference 64 for vl-convert), and matplotlib's PDF and SVG differed by their timestamps until `metadata={"CreationDate": None}` / `{"Date": None}` was set.

### Adding a section, maintenance

| | (a) Quarto | (b) Typst + Jinja |
| --- | --- | --- |
| Per section | one `.qmd` (54 lines each here): prose once, Python cells, and each table and figure emitted twice (`content-visible when-format="html"` / `unless-format="html"`, because the PDF gets the capped table and the static figure) | one Typst file (19-26 lines) and one Jinja template (11-18 lines): prose twice, data and captions once in Python |
| Project code | `_quarto.yml` 51 + pre-render 37 + helpers 38 + render script 25 lines | `build.py` 119 + `lib.typ` 42 + `base.html`/`macros.html` 26 lines |
| Moving parts | Quarto, Deno, pandoc, dart-sass, esbuild, Jupyter kernel, Typst, the `orange-book` template, a CSL file | Typst, Jinja2 |
| Known traps | the five workarounds above; table layout lives in pandoc and in raw Typst injected through `include-in-header` | none found; the HTML references and numbering are hand-written (about 15 lines) |

## P2: figures

### The shared style

The prototype's style module (114 lines) sets:

- **Widths:** 89 mm (single column) up to 14 samples, 183 mm (double) above; heights follow the content, at most 247 mm.
- **Type:** Inter 4.1 (conda-forge `font-inter`, static TTFs) at 7 pt base and 6 pt ticks. Rows of 100 or more samples use 5 pt, and rows above 100 use 4.5 pt (see below).
- **Font embedding:** TrueType (type 42) in PDF. In SVG the text stays text, and the font is added as WOFF2 `@font-face` data URIs, so the file looks the same without the font installed and stays editable. The three faces add about 430 KB per SVG, so 11d could embed only the faces a figure uses.
- **Colour:** Okabe-Ito for categories (phyla, read fates); viridis for quantities, with its top 5 % dropped; light grey for "not detected"; vermillion for low-depth labels and the cutoff line.
- **Files:** PDF, SVG, 600 dpi PNG and a source CSV (the plotted values, long format) for every figure.
- **Samples:** in group then natural order, with group brackets and separators. Low-depth samples are hatched, and their labels are coloured and marked with "†".
- **Large sets:** above 60 samples a heatmap is drawn with samples as rows. With samples as columns, 120 samples at 183 mm leave 1.2 mm per column, less than one line of 5 pt type.

One finding for 11d: the `font-inter` package also ships CFF `.otf` files and a variable font. matplotlib embedded some glyphs from the `.otf` (poppler reported a font-type mismatch) until only the static TTFs were left in its font list, and Typst picked `InterVariable` until it was given a folder with the static TTFs only.

### Static output against paper needs

| Need | (i) matplotlib | (ii) Vega-Lite + vl-convert |
| --- | --- | --- |
| Exact width | exact (89.0 / 183.0 mm; `bbox_inches="tight"` avoided) | 88.9-89.0 / 183.1 mm after a measure-and-refit loop (2-10 renders); a step-sized plot ignores the width until its cells outgrow the legend |
| Exact height, 120 samples | 230-247 mm, capped at the page | 225-260 mm, not controlled (two heatmaps over 247 mm) |
| Grouped sample labels | brackets and separators | facets (heatmaps, and a concat of two row facets for read fate) |
| Hatching for low depth | heatmap columns and read-fate bars (`////`) | heatmaps: one white diagonal per cell, which lines up into stripes; bars: not possible, grey only |
| Legends outside the plot | anywhere; phylum legend beside the plot, colour bar and key below | per sub-chart only; a legend beside the phylum strip widened that sub-chart until it was moved to the far side |
| Fonts embedded | PDF TrueType, SVG text plus `@font-face` | PDF TrueType + CFF, SVG text plus `@font-face`; PNG needed `vl_convert.register_font_directory` |
| Render time, three figures (4 / 120 samples) | 8.7-14.5 / 43-58 s | static only 4-10 / 25-35 s (incl. refit passes) |

The matplotlib figures were also the more finished ones, with a log axis labelled "1 M, 10 M", a phylum colour strip next to the genus labels and a separate key. Phase 11d still has to tune a few things: the 1 M tick disappears when the axis starts above it, and at 120 samples the 4.5 pt row labels are below Nature's 5 pt minimum. That size is the limit of one page: 5 pt labels fit about 110 samples on a 247 mm page, so larger runs need the figure split by group or labelled by group only.

### Interactive output

| | Vega-Lite (both options) |
| --- | --- |
| Hover | tooltips with sample, group, taxon/type, value, phylum |
| Filtering | a group menu (`binding_select`) under the heatmaps; tables filter and sort through the shared script |
| Zoom | drag/scroll zoom on the read-fate log axis (`bind="scales"`); not on heatmaps, whose axes are ordinal |
| Screen size | the print spec (89 mm, 6 pt) is unreadable in a browser, so the HTML version is the same spec function with type and cells × 1.8 |
| Libraries, offline | `vl_convert.javascript_bundle()`: vega-embed + Vega-Lite + Vega in one 880 KB script, once per site |
| Spec size with inline data (4 / 120 samples) | taxonomy heatmap 14 / 877 KB, read fate 5 / 74 KB, SARG heatmap 16 / 369 KB |

A size budget for 11d: at 120 samples the largest inline spec is 877 KB. A budget of about 1 MB per figure would keep the interactive version up to about 140 samples, with the static figure above that.

### Plotly

`plotly` 7.1.0 with `python-kaleido` 1.5.0 (a Seqera image of both, 437 MB): `fig.write_image("x.pdf")` raised `RuntimeError` asking for Google Chrome ("Either download and install Chrome yourself … or install it … by running: `plotly_get_chrome`"). Out, per the phase rule.

## What Phase 11d would carry over

If the recommendations are taken:

- **Image:** the "(b) + figures" image above, re-pinned everywhere the shared report image appears.
- **Fonts:** the static Inter TTFs only, registered explicitly in matplotlib (font list filtered), Typst (`--font-path`, `--ignore-system-fonts`) and vl-convert (`register_font_directory`). `XDG_CACHE_HOME` goes in the task directory.
- **Figures:** deterministic metadata (`CreationDate`/`Date` set to `None`); 89 / 183 mm; samples as rows above 60; a split or group-only labels above about 110.
- **PDF:** the 60-row cap with a pointer to the TSV, landscape pages for wide tables, and a curated tool → reference table for the methods section.
- **HTML:** interactive Vega-Lite for the read-fate and heatmap figures, the static figure above the size budget, and the bundled Vega script.
