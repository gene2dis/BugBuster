#!/bin/bash
#
# Refreshes docs/images/bugbuster_overview.svg's embedded fonts and renders
# bugbuster_overview.png from it, in two pinned containers (nothing installed on
# the host; needs docker). Run it after editing the SVG.
#
#  1. The report image (bin/figure_style.py's Inter, fontTools 4.68): subsets
#     Inter Regular and SemiBold to the characters the SVG's text uses and writes
#     them as WOFF @font-face data into its <style id="fonts"> block, so the SVG
#     shows Inter where Inter is not installed (GitHub) and its text stays text.
#  2. librsvg 2.62.3 + font-inter 4.1 (Seqera Containers, conda-forge
#     librsvg=2.62.3 font-inter=4.1 fontconfig; sha256:d5ea0d67...): rsvg-convert
#     at 1.75x (2,240 px wide). librsvg takes Inter from fontconfig, not from
#     the @font-face data.
#
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REPORT_IMAGE="community.wave.seqera.io/library/python_pandas_numpy_matplotlib-base_pruned:6bff5ca407ba1454"
RSVG_IMAGE="community.wave.seqera.io/library/librsvg_font-inter_fontconfig:24f32e0e8f33ce87"
SVG="bugbuster_overview.svg"
PNG="bugbuster_overview.png"

docker run --rm -i -u "$(id -u):$(id -g)" -v "${SCRIPT_DIR}:/work" -w /work "${REPORT_IMAGE}" python -I - "${SVG}" <<'EOF'
import base64, html, io, re, sys
from pathlib import Path
from fontTools import subset
from fontTools.ttLib import TTFont

svg_path = Path(sys.argv[1])
svg = svg_path.read_text()
body = re.sub(r'<style id="fonts">.*?</style>', '', svg, flags=re.S)
chars = set()
for text in re.findall(r"<text\b[^>]*>(.*?)</text>", body, flags=re.S):
    chars.update(html.unescape(re.sub(r"<[^>]+>", "", text)))
chars.add(" ")
fonts = Path(sys.prefix) / "fonts"
rules = []
for face, weight in (("Regular", 400), ("SemiBold", 600)):
    font = TTFont(fonts / f"Inter-{face}.ttf")
    options = subset.Options()
    options.flavor = "woff"
    options.layout_features = ["kern", "liga", "calt"]
    options.name_IDs = ["*"]
    sub = subset.Subsetter(options)
    sub.populate(text="".join(sorted(chars)))
    sub.subset(font)
    out = io.BytesIO()
    font.flavor = "woff"
    font.save(out)
    data = base64.b64encode(out.getvalue()).decode()
    rules.append(f"@font-face{{font-family:Inter;font-weight:{weight};font-style:normal;"
                 f"src:url(data:font/woff;base64,{data}) format('woff');}}")
svg = re.sub(r'<style id="fonts">.*?</style>', lambda m: '<style id="fonts">' + "".join(rules) + "</style>",
             svg, count=1, flags=re.S)
svg_path.write_text(svg)
print(f"{svg_path}: Inter 400 + 600 subset to {len(chars)} characters, {svg_path.stat().st_size // 1024} KB")
EOF

docker run --rm -u "$(id -u):$(id -g)" -v "${SCRIPT_DIR}:/work" -w /work "${RSVG_IMAGE}" \
    rsvg-convert --zoom 1.75 --format png --output "${PNG}" "${SVG}"
echo "${PNG}: $(( $(stat -c %s "${SCRIPT_DIR}/${PNG}") / 1024 )) KB"
