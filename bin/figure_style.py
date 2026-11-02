"""Shared figure style of BugBuster: the run summary's figures (bin/summary_figures.py) and the stage figures
(report_unify.py, taxonomy_report.py, taxonomy_phyloseq.py, blobplot.py, arg_blobplot.py, bin_quality_report.py,
bin_tax_report.py) all draw with it.

- widths: journal single column 89 mm and double column 183 mm, at the exact size (never bbox_inches="tight",
  which changes it); heights follow the content, at most 247 mm (an A4 page less its margins)
- font: Inter, from the static TTFs of conda-forge's font-inter (in the report image's <prefix>/fonts). The package
  also ships CFF .otf files and a variable font; only the static TTFs are registered, so a PDF embeds TrueType
  (pdf.fonttype 42) and every run picks the same files. SVG text stays text, with the faces it uses embedded as
  WOFF2 @font-face data, so the file renders the same without Inter installed and stays editable
- colour: Okabe-Ito for categories, viridis (its brightest 5 % dropped) for quantities, light grey for 'not
  detected', vermillion for low-depth samples and thresholds
- type: 7 pt base, 6 pt ticks and legends, 8 pt titles; sample labels 6 pt up to 60 samples, 5 pt up to 110
  (Nature's minimum), and above that no sample labels, only group labels (the interactive HTML version names them)
- samples: natural order (V2 before V10) inside their group, group brackets; above 60 samples a heatmap has samples
  as rows (with samples as columns, 120 of them at 183 mm leave 1.2 mm each); low-depth samples hatched and their
  labels marked with a dagger
- files: <key>.pdf, <key>.svg, <key>.png (600 dpi) and <key>.csv (the plotted values) with no timestamp in the
  metadata, so a rerun on the same data writes the same PDF and SVG bytes
"""

import base64
import os
import re
import sys
from pathlib import Path

MM = 1 / 25.4
WIDTH_MM = {"single": 89, "double": 183}
MAX_H_MM = 247
FONT = "Inter"
PNG_DPI = 600
OKABE_ITO = ["#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000"]
# Okabe-Ito without black and yellow first, for bars and strips where a black or pale yellow block reads badly
CATEGORY = ["#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00", "#F0E442", "#000000"]
OTHER = "#BBBBBB"
NOT_DETECTED = "#E8E8E8"
LOW_DEPTH = "#D55E00"
DAGGER = "†"
BASE_PT, SMALL_PT, TITLE_PT = 7, 6, 8
VERTICAL_ABOVE = 60   # samples as rows above this many samples
LABEL_5PT_UP_TO = 110  # 5 pt labels fit about 110 samples on a 247 mm page; above, groups are labelled only
_FACES = {("normal", "normal"): "Regular", ("bold", "normal"): "Bold", ("normal", "italic"): "Italic",
          ("bold", "italic"): "BoldItalic", ("500", "normal"): "Medium", ("600", "normal"): "SemiBold"}
_READY = {"done": False, "inter": False}


def font_dir():
    """Directory with Inter's static TTF and WOFF2 files: $BB_FONT_DIR, else <python prefix>/fonts (conda's
    font-inter installs there; CONDA_PREFIX is not set inside the pipeline's containers)."""
    return Path(os.environ.get("BB_FONT_DIR") or Path(sys.prefix) / "fonts")


def font_files(weights=("Regular", "Bold", "Italic", "BoldItalic", "Medium", "SemiBold"), ext="ttf"):
    d = font_dir()
    return [d / f"{FONT}-{w}.{ext}" for w in weights if (d / f"{FONT}-{w}.{ext}").exists()]


def setup():
    """matplotlib with the Agg backend, Inter registered (static TTFs only) and the shared rcParams. Returns True
    when Inter is available; without it (a standalone run outside the report image) the figures use matplotlib's
    default font, with a warning."""
    import matplotlib
    matplotlib.use("Agg")
    if _READY["done"]:
        return _READY["inter"]
    from matplotlib import font_manager
    files = font_files()
    fm = font_manager.fontManager
    if files:
        for f in files:
            fm.addfont(str(f))
        keep = {str(f) for f in files}
        # drop every other Inter file matplotlib may know (.otf, InterVariable): with them, glyphs can come from a
        # CFF font, and the PDF embeds a font of another type
        fm.ttflist = [e for e in fm.ttflist if not (e.name.startswith(FONT) and e.fname not in keep)]
    else:
        print(f"WARNING: no {FONT} TTF in {font_dir()}: figures use matplotlib's default font", file=sys.stderr)
    matplotlib.rcParams.update({
        "font.family": FONT if files else "DejaVu Sans", "font.size": BASE_PT, "axes.titlesize": TITLE_PT,
        "axes.titleweight": "bold", "axes.labelsize": BASE_PT, "xtick.labelsize": SMALL_PT,
        "ytick.labelsize": SMALL_PT, "legend.fontsize": SMALL_PT, "legend.title_fontsize": SMALL_PT,
        "figure.titlesize": TITLE_PT, "figure.titleweight": "bold", "axes.linewidth": 0.5,
        "xtick.major.width": 0.5, "ytick.major.width": 0.5, "xtick.minor.width": 0.4, "ytick.minor.width": 0.4,
        "xtick.major.size": 2, "ytick.major.size": 2, "xtick.minor.size": 1.2, "ytick.minor.size": 1.2,
        "lines.linewidth": 0.75, "patch.linewidth": 0.5, "hatch.linewidth": 0.4, "pdf.fonttype": 42,
        "svg.fonttype": "none", "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 100,
        "savefig.dpi": PNG_DPI, "savefig.facecolor": "white", "image.cmap": "viridis", "legend.frameon": False,
        "axes.titlelocation": "left", "figure.constrained_layout.h_pad": 0.02,
        "figure.constrained_layout.w_pad": 0.02, "svg.hashsalt": "bugbuster",
    })
    _READY.update(done=True, inter=bool(files))
    return _READY["inter"]


def width_class(n_columns, single_up_to=14):
    """'single' (89 mm) for up to single_up_to sample columns or bars, else 'double' (183 mm)."""
    return "single" if n_columns <= single_up_to else "double"


def figsize(width="double", height_mm=100):
    """(width, height) in inches; width is 'single', 'double' or a number of mm; the height is capped."""
    w = WIDTH_MM[width] if isinstance(width, str) else float(width)
    return w * MM, max(20, min(float(height_mm), MAX_H_MM)) * MM


def sample_layout(n):
    """How n samples are drawn on a sample axis: label font size (pt; 0 = no sample labels, groups only), whether a
    heatmap puts samples on rows, and the size of one sample's row or column in mm (a label's height)."""
    fs = SMALL_PT if n <= VERTICAL_ABOVE else (5 if n <= LABEL_5PT_UP_TO else 0)
    cell = (fs or 5) * 0.3528 * 1.05
    return {"fs": fs, "vertical": n > VERTICAL_ABOVE, "cell_mm": cell, "labels": fs > 0}


def viridis(n=256):
    """viridis without its brightest 5 %, so the top of the scale never looks like the 'not detected' grey."""
    import numpy as np
    from matplotlib import colormaps
    from matplotlib.colors import ListedColormap
    cm = ListedColormap(colormaps["viridis"](np.linspace(0, 0.95, n)))
    cm.set_bad(NOT_DETECTED)
    return cm


def category_colors(names, other=("Other", "other", "no-hit", "unassigned", "unclassified", "multiple classes")):
    """A colour per name in the order given: Okabe-Ito, then repeated in lighter shades; names meaning 'the rest'
    are grey and take no slot."""
    import matplotlib.colors as mc
    out, i = {}, 0
    for n in names:
        if n in other:
            out[n] = {"no-hit": "#E4E4E4", "unclassified": "#DADADA", "unassigned": "#DADADA",
                      "multiple classes": "#5A5A5A"}.get(n, OTHER)
            continue
        base = mc.to_rgb(CATEGORY[i % len(CATEGORY)])
        k = i // len(CATEGORY)
        out[n] = mc.to_hex(tuple(c + (1 - c) * min(0.6, 0.35 * k) for c in base))
        i += 1
    return out


def short_count(n):
    """10000000 -> '10 M', 1500000 -> '1.5 M', 250000 -> '250 k'."""
    for div, suffix in ((1e9, "G"), (1e6, "M"), (1e3, "k")):
        if abs(n) >= div:
            return f"{n / div:.3g} {suffix}"
    return f"{n:.0f}"


def count_axis(axis):
    """Log axis of read counts labelled 1 k / 1 M / 10 M, decades only."""
    from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter
    axis.set_major_locator(LogLocator(base=10))
    axis.set_major_formatter(FuncFormatter(lambda v, _: short_count(v) if v >= 1 else f"{v:g}"))
    axis.set_minor_formatter(NullFormatter())


def clip_label(text, n=40):
    """A long label cut to n characters with an ellipsis (the source CSV and the interactive version keep it)."""
    text = str(text)
    return text if len(text) <= n else text[: n - 1].rstrip() + "\u2026"


def label_text(s, low):
    return f"{s} {DAGGER}" if s in low else str(s)


def sample_ticks(ax, samples, low, axis, fs, rotation=90):
    """Sample tick labels; low-depth samples in vermillion with a dagger. fs 0: no labels."""
    pos = range(len(samples))
    if not fs:
        (ax.set_xticks if axis == "x" else ax.set_yticks)([])
        return
    labels = [label_text(s, low) for s in samples]
    if axis == "x":
        ax.set_xticks(pos, labels, rotation=rotation, fontsize=fs)
        ticks = ax.get_xticklabels()
    else:
        ax.set_yticks(pos, labels, fontsize=fs)
        ticks = ax.get_yticklabels()
    for t, s in zip(ticks, samples):
        if s in low:
            t.set_color(LOW_DEPTH)


def group_spans(samples, groups):
    """[(group, first index, last index)] of consecutive samples with the same group."""
    spans, start = [], 0
    for i in range(1, len(samples) + 1):
        if i == len(samples) or groups[i] != groups[start]:
            spans.append((groups[start], start, i - 1))
            start = i
    return spans


def group_brackets(ax, samples, groups, axis, separator="white", side=None, offset=0.01):
    """Group names beside brackets outside the sample axis, with a separator line between groups. axis 'x': the
    brackets above the plot; 'y': right of it (side 'left' puts them left)."""
    if not any(groups):
        return
    for i, (g, a, b) in enumerate(group_spans(samples, groups)):
        if axis == "x":
            if i and separator:
                ax.axvline(a - 0.5, color=separator, lw=0.8)
            y = 1 + offset
            ax.annotate("", xy=(a - 0.4, y), xytext=(b + 0.4, y), xycoords=("data", "axes fraction"),
                        arrowprops={"arrowstyle": "-", "lw": 0.5, "shrinkA": 0, "shrinkB": 0},
                        annotation_clip=False)
            ax.text((a + b) / 2, y + 0.004, str(g), transform=ax.get_xaxis_transform(), ha="center",
                    va="bottom", fontsize=SMALL_PT, clip_on=False)
        else:
            if i and separator:
                ax.axhline(a - 0.5, color=separator, lw=0.8)
            x = -offset if side == "left" else 1 + offset
            ax.annotate("", xy=(x, a - 0.4), xytext=(x, b + 0.4), xycoords=("axes fraction", "data"),
                        arrowprops={"arrowstyle": "-", "lw": 0.5, "shrinkA": 0, "shrinkB": 0},
                        annotation_clip=False)
            ax.text(x + (-0.008 if side == "left" else 0.008), (a + b) / 2, str(g),
                    transform=ax.get_yaxis_transform(), ha="right" if side == "left" else "left", va="center",
                    rotation=90 if side == "left" else 270, fontsize=SMALL_PT, clip_on=False)


def hatch_samples(ax, samples, hatch_set, n_features, axis):
    """Hatch the cells of the samples in hatch_set ('not detected' there is unknown). axis: the sample axis."""
    from matplotlib.patches import Rectangle
    for j, s in enumerate(samples):
        if s in hatch_set:
            xy, w, h = ((j - 0.5, -0.5), 1, n_features) if axis == "x" else ((-0.5, j - 0.5), n_features, 1)
            ax.add_patch(Rectangle(xy, w, h, fill=False, hatch="////", edgecolor="white", lw=0))


def low_depth_handle(min_pairs=None, label=None):
    from matplotlib.patches import Patch
    text = label or (f"{DAGGER} below {short_count(min_pairs)} clean read pairs" if min_pairs else
                     f"{DAGGER} low depth")
    return Patch(facecolor="#9A9A9A", hatch="////", edgecolor="white", lw=0, label=text)


def faces_used(svg_text):
    """Inter faces the text of an SVG uses, as file stems (Regular, Bold, ...). matplotlib writes each text's font
    in its style attribute: font-size, font-family and, when not regular, font-weight / font-style."""
    used = set()
    for m in re.finditer(r'style="([^"]*font-family:[^"]*)"', svg_text):
        css = dict((k.strip(), v.strip()) for k, _, v in (p.partition(":") for p in m.group(1).split(";")) if k)
        if FONT not in css.get("font-family", ""):
            continue
        style = "italic" if css.get("font-style") in ("italic", "oblique") else "normal"
        weight = {"400": "normal", "700": "bold", "800": "bold", "900": "bold"}.get(css.get("font-weight", "normal"),
                                                                                  css.get("font-weight", "normal"))
        used.add(_FACES.get((weight, style), "Regular" if weight in ("normal", "300", "200", "100") else "Bold"))
    return used or {"Regular"}


def embed_svg_fonts(path):
    """Add the WOFF2 data of the Inter faces an SVG's text uses as @font-face rules."""
    text = Path(path).read_text()
    faces = []
    for stem in sorted(faces_used(text)):
        f = font_dir() / f"{FONT}-{stem}.woff2"
        if not f.exists():
            continue
        style = "italic" if "Italic" in stem else "normal"
        weight = {"Bold": 700, "BoldItalic": 700, "Medium": 500, "SemiBold": 600}.get(stem, 400)
        faces.append(f"@font-face{{font-family:'{FONT}';font-weight:{weight};font-style:{style};"
                     f"src:url(data:font/woff2;base64,{base64.b64encode(f.read_bytes()).decode()}) "
                     "format('woff2')}")
    if faces:
        text = re.sub(r"(<svg\b[^>]*>)", lambda m: m.group(1) + "<style>" + "".join(faces) + "</style>", text,
                      count=1)
        Path(path).write_text(text)


def save(fig, outdir, key, source=None, formats=("pdf", "svg", "png")):
    """<key>.pdf / .svg / .png (600 dpi) at the figure's exact size, plus <key>.csv from source (a DataFrame of the
    plotted values), with no timestamps in the metadata. Returns the paths written."""
    import matplotlib.pyplot as plt
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    meta = {"pdf": {"CreationDate": None, "Producer": None}, "svg": {"Date": None}, "png": {"Software": None}}
    paths = []
    for ext in formats:
        p = outdir / f"{key}.{ext}"
        fig.savefig(p, facecolor="white", metadata=meta.get(ext, {}))
        paths.append(p)
    plt.close(fig)
    if "svg" in formats and _READY["inter"]:
        embed_svg_fonts(outdir / f"{key}.svg")
    if source is not None:
        p = outdir / f"{key}.csv"
        source.to_csv(p, index=False)
        paths.append(p)
    return paths


def placeholder(outdir, key, message, formats=("pdf", "svg", "png")):
    """A single-column figure holding only a message (a stage with nothing to plot), with an empty source CSV."""
    import matplotlib.pyplot as plt
    import pandas as pd
    setup()
    fig = plt.figure(figsize=figsize("single", 30))
    fig.text(0.5, 0.5, message, ha="center", va="center", fontsize=BASE_PT)
    return save(fig, outdir, key, pd.DataFrame({"message": [message]}), formats)
