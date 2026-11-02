"""The run summary's HTML report: a multi-page offline site (<outdir>/report/) from the Jinja templates in
assets/report/html/ and the report model of bin/summary_report.py, one page per section present. Tables are
sortable and filterable (assets/report/static/tables.js, no library); the large figures are interactive Vega-Lite
(bin/summary_interactive.py) with the Vega libraries bundled once in report/assets/; every figure links its static
SVG / PDF / PNG and its CSV in ../figures/. Nothing is loaded from the network, so the site works from file://.
<outdir>/results_summary.html sends the browser to report/index.html."""

import shutil
import time
from pathlib import Path

import figure_style as fs
from summary_common import die, warn

WOFF2_FACES = ("Regular", "Bold", "Italic", "SemiBold")
REDIRECT = """<!doctype html>
<html lang="en"><head><meta charset="utf-8"><title>{title}</title>
<meta http-equiv="refresh" content="0; url=report/index.html"></head>
<body><p>The report is in <a href="report/index.html">report/index.html</a>.</p></body></html>
"""


def write_site(outdir, model, assets):
    try:
        from jinja2 import Environment, FileSystemLoader, StrictUndefined, select_autoescape
    except ImportError:
        die("the HTML report needs jinja2 (conda-forge jinja2; the pipeline's report image has it)")
    t0 = time.time()
    outdir = Path(outdir)
    site = outdir / "report"
    if site.exists():
        shutil.rmtree(site)
    (site / "assets" / "fonts").mkdir(parents=True)
    static = Path(assets) / "static"
    for f in sorted(static.iterdir()):
        shutil.copy(f, site / "assets" / f.name)
    for f in fs.font_files(WOFF2_FACES, "woff2"):
        shutil.copy(f, site / "assets" / "fonts" / f.name)
    specs, bundle = {}, False  # interactive figures: bin/summary_interactive.py (next commit)
    env = Environment(loader=FileSystemLoader(str(Path(assets) / "html")), autoescape=select_autoescape(["html"]),
                      undefined=StrictUndefined, trim_blocks=True, lstrip_blocks=True, keep_trailing_newline=True)
    nav = [(s["page"], s["title"]) for s in model["sections"]]
    for sec in model["sections"]:
        html = env.get_template(f"{sec['id']}.html").render(run=model["run"], s=sec, nav=nav, specs=specs,
                                                            bundle=bundle, sections=model["sections"])
        (site / sec["page"]).write_text(html)
    (outdir / "results_summary.html").write_text(REDIRECT.format(title=model["run"]["title"]))
    print(f"report/ written: {len(model['sections'])} pages ({time.time() - t0:.1f} s)")
    return site
