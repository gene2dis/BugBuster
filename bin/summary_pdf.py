"""The run summary's PDF: Typst templates (assets/report/typst/) reading the report model of
bin/summary_report.py as data.json, compiled with the static Inter TTFs only (--font-path, --ignore-system-fonts)
and no Typst packages, so it needs no network and writes nothing outside its build directory."""

import json
import os
import shutil
import subprocess
import time
from pathlib import Path

import figure_style as fs
from summary_common import warn

FONT_FACES = ("Regular", "Bold", "Italic", "BoldItalic", "Medium", "SemiBold")


def write_pdf(outdir, model, assets, fig_dir, keep_build=None):
    """Compile <outdir>/results_summary.pdf; returns its path, or None (with a warning) when typst is missing or
    fails, so the rest of the summary is still written."""
    if keep_build is None:  # BB_KEEP_PDF_BUILD=1 keeps .pdf_build/ to debug the templates
        keep_build = os.environ.get("BB_KEEP_PDF_BUILD") == "1"
    typst = shutil.which("typst")
    if typst is None:
        warn("typst is not installed: no results_summary.pdf (the pipeline's report image has it)")
        return None
    build = Path(outdir) / ".pdf_build"
    if build.exists():
        shutil.rmtree(build)
    shutil.copytree(Path(assets) / "typst", build)
    (build / "fonts").mkdir()
    for f in fs.font_files(FONT_FACES):  # static TTFs only: with InterVariable in the folder Typst picks it
        shutil.copy(f, build / "fonts" / f.name)
    (build / "figures").mkdir()
    for sec in model["sections"]:
        for key in sec["figures"]:
            shutil.copy(Path(fig_dir) / f"{key}.svg", build / "figures" / f"{key}.svg")
    methods = next((s for s in model["sections"] if s["id"] == "methods"), None)
    (build / "references.bib").write_text(methods["facts"].get("bibtex", "") if methods else "")
    data = json.loads(json.dumps(model, default=str))
    for sec in data["sections"]:
        for f in sec["figures"].values():
            f.pop("interactive", None)  # the HTML's data, not the PDF's
        for t in sec["tables"].values():
            if t:
                t.pop("rows", None)  # the PDF reads its own capped rows (t["pdf"])
    (build / "data.json").write_text(json.dumps(data))
    out = Path(outdir).resolve() / "results_summary.pdf"
    t0 = time.time()
    env = dict(os.environ)
    if not env.get("SOURCE_DATE_EPOCH", "").strip().isdigit():
        env.pop("SOURCE_DATE_EPOCH", None)  # typst (and fontconfig) reject an empty or malformed value
    res = subprocess.run([typst, "compile", "--ignore-system-fonts", "--font-path", "fonts", "report.typ", str(out)],
                         cwd=build, capture_output=True, text=True, env=env)
    if res.returncode != 0:
        warn(f"typst failed, no results_summary.pdf:\n{res.stderr.strip()}")
        return None
    for line in res.stderr.strip().splitlines():
        if line.strip():
            warn(f"typst: {line}")
    if not keep_build:
        shutil.rmtree(build)
    print(f"results_summary.pdf written ({time.time() - t0:.1f} s)")
    return out
