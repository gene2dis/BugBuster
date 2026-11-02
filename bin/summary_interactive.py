"""Interactive Vega-Lite (Altair) versions of the run summary's large figures, for the HTML report: hover for
values, a group menu when the samples have groups, and zoom on quantitative axes. Built from the same data as the
static matplotlib figure (FigureWriter's 'interactive' record), at a screen scale (type and cells 1.8x the print
size). A spec whose JSON is larger than SPEC_BUDGET is not embedded: the page shows the static figure instead.
"""

import json

import numpy as np

import figure_style as fs
from summary_common import warn

SPEC_BUDGET = 1_000_000  # bytes of spec JSON (data inline) per figure
SCALE = 1.8


def available():
    try:
        import altair  # noqa: F401 (availability check)
        import vl_convert  # noqa: F401
    except ImportError:
        return False
    return True


def config():
    small, base = fs.SMALL_PT * SCALE, fs.BASE_PT * SCALE  # px: the print sizes in pt, scaled for the screen
    return {
        "font": fs.FONT, "background": "white", "view": {"stroke": None},
        "axis": {"labelFontSize": small, "titleFontSize": base, "labelFont": fs.FONT, "titleFont": fs.FONT,
                 "domainWidth": 0.8, "grid": False, "tickWidth": 0.8, "tickSize": 3, "titleFontWeight": "normal"},
        "legend": {"labelFontSize": small, "titleFontSize": small, "titleFontWeight": "normal", "labelFont": fs.FONT,
                   "titleFont": fs.FONT, "gradientThickness": 10},
        "header": {"labelFontSize": small, "labelFont": fs.FONT},
        "range": {"category": fs.CATEGORY},
    }


def _low_expr(low):
    return f"indexof({json.dumps(sorted(low))}, datum.value) >= 0"


def _sample_axis(alt, low, angle):
    return alt.Axis(labelExpr=f"{_low_expr(low)} ? datum.value + ' {fs.DAGGER}' : datum.value",
                    labelColor={"condition": {"test": _low_expr(low), "value": fs.LOW_DEPTH}, "value": "black"},
                    ticks=False, domain=False, labelAngle=angle, labelOverlap=False)


def _group_select(alt, groups):
    names = list(dict.fromkeys(g for g in groups if g != ""))
    if len(names) < 2:
        return None
    return alt.selection_point(fields=["group"], name="grp", bind=alt.binding_select(
        options=[None] + names, labels=["all"] + names, name="Group "))


def heatmap(d):
    """Heatmap: samples as columns (rows above 60 samples), grey 'not detected' cells, a white diagonal over the
    cells of hatched samples (Vega-Lite has no fill patterns), tooltips, a group menu."""
    import altair as alt
    df = d["frame"].copy()
    feat, val = d["feature"], d["value"]
    samples, features = d["samples"], d["features"]
    df["log10"] = np.where(df[val] > 0, np.log10(df[val].where(df[val] > 0)), np.nan)
    df[feat] = df[feat].astype(str)
    vertical = len(samples) > fs.VERTICAL_ABOVE
    step = 9 * SCALE  # a feature row: the 6 pt label height at screen scale
    sstep = step if vertical else max(step, min(36, 720 / max(len(samples), 1)))  # wider columns for few samples
    s_enc = {"field": "sample", "type": "ordinal", "sort": samples, "title": None}
    f_enc = {"field": feat, "type": "ordinal", "sort": [str(x) for x in features], "title": None}
    s_axis = _sample_axis(alt, d["low"], 0 if vertical else -90)
    f_axis = alt.Axis(ticks=False, domain=False, labelLimit=320, labelFontStyle="italic" if d.get("italic") else
                      "normal", orient="top" if vertical else "left", labelAngle=-90 if vertical else 0)
    f_axis_main = f_axis
    if vertical:
        enc = {"y": alt.Y(**s_enc, axis=s_axis), "x": alt.X(**f_enc, axis=f_axis_main)}
    else:
        enc = {"x": alt.X(**s_enc, axis=s_axis), "y": alt.Y(**f_enc, axis=f_axis_main)}
    tip = [alt.Tooltip("sample:N"), alt.Tooltip(f"{feat}:N")]
    if any(g != "" for g in d["groups"]):
        tip.insert(1, alt.Tooltip("group:N"))
    if d.get("row_group"):
        tip.append(alt.Tooltip(f"{d['row_group']}:N", title=d.get("row_group_title") or "group"))
    tip.append(alt.Tooltip(f"{val}:Q", format=".3g"))
    base = alt.Chart(df)
    sel = _group_select(alt, d["groups"])
    bg = base.mark_rect(color=fs.NOT_DETECTED).encode(**enc, tooltip=tip)
    cells = base.transform_filter("isValid(datum.log10)").mark_rect().encode(
        **enc, tooltip=tip,
        color=alt.Color("log10:Q", scale=alt.Scale(scheme="viridis", domain=d["domain"], clamp=True),
                        legend=alt.Legend(title=d["legend"], orient="bottom", direction="horizontal",
                                          gradientLength=160)))
    layers = [bg, cells]
    if df["hatched"].any():
        b0, b1 = {"bandPosition": 0}, {"bandPosition": 1}
        if vertical:
            dg = {"x": alt.X(**f_enc, **b0), "x2": alt.X2(feat, **b1), "y": alt.Y(**s_enc, **b1),
                  "y2": alt.Y2("sample", **b0)}
        else:
            dg = {"x": alt.X(**s_enc, **b0), "x2": alt.X2("sample", **b1), "y": alt.Y(**f_enc, **b1),
                  "y2": alt.Y2(feat, **b0)}
        layers.append(base.transform_filter("datum.hatched").mark_rule(color="white", strokeWidth=1, opacity=0.9)
                      .encode(**dg))
    if sel is not None:
        layers = [layers[0].add_params(sel).transform_filter(sel)] + [x.transform_filter(sel) for x in layers[1:]]
    chart = alt.layer(*layers).properties(width=alt.Step(step if vertical else sstep),
                                          height=alt.Step(sstep if vertical else step))
    if d.get("row_group"):
        sdf = df.drop_duplicates(feat)[[feat, d["row_group"]]]
        names = list(dict.fromkeys(sdf[d["row_group"]]))
        pal = fs.category_colors(names)
        strip_enc = alt.X(**f_enc, axis=None) if vertical else alt.Y(**f_enc, axis=None)
        strip = alt.Chart(sdf).mark_rect().encode(
            **{"x" if vertical else "y": strip_enc},
            color=alt.Color(f"{d['row_group']}:N", sort=names, scale=alt.Scale(domain=names,
                                                                                range=[pal[n] for n in names]),
                            legend=alt.Legend(title=d.get("row_group_title") or None, orient="bottom", columns=4,
                                              symbolType="square")),
            tooltip=[alt.Tooltip(f"{feat}:N"), alt.Tooltip(f"{d['row_group']}:N")])
        strip = strip.properties(**({"height": 8, "width": alt.Step(step)} if vertical else
                                    {"width": 8, "height": alt.Step(step)}))
        # the strip on the far side of the labels: its legend widens (or lengthens) only its own column
        chart = (alt.vconcat(chart, strip, spacing=2) if vertical else alt.hconcat(chart, strip, spacing=2)) \
            .resolve_scale(color="independent")
    return chart


def read_fate(d):
    """Read fate: stacked shares of raw read pairs and clean read pairs on a log axis (drag / scroll to zoom)."""
    import altair as alt
    f = d["frame"].copy()
    samples, low = d["samples"], d["low"]
    step = 9 * SCALE
    long = f.melt(id_vars=["sample", "group", "low_depth"], value_vars=["pct_clean", "pct_removed_host",
                                                                        "pct_removed_qc"],
                  var_name="fate", value_name="pct")
    names = {"pct_clean": "clean", "pct_removed_host": "removed: host / contaminant",
             "pct_removed_qc": "removed: read QC"}
    long["fate"] = long["fate"].map(names)
    long["order"] = long["fate"].map({v: i for i, v in enumerate(names.values())})
    y = alt.Y("sample:O", sort=samples, title=None, axis=_sample_axis(alt, low, 0))
    tip_g = [alt.Tooltip("group:N")] if any(g != "" for g in d["groups"]) else []
    sel = _group_select(alt, d["groups"])
    left = alt.Chart(long).mark_bar().encode(
        x=alt.X("pct:Q", stack="zero", title="% of raw read pairs", scale=alt.Scale(domain=[0, 100])),
        y=y, order="order:Q",
        color=alt.Color("fate:N", scale=alt.Scale(domain=list(names.values()),
                                                  range=[fs.CATEGORY[0], fs.CATEGORY[1], "#9A9A9A"]),
                        legend=alt.Legend(title=None, orient="bottom", columns=3, labelLimit=400,
                                          columnPadding=18)),
        tooltip=[alt.Tooltip("sample:N")] + tip_g + [alt.Tooltip("fate:N"), alt.Tooltip("pct:Q", format=".2f")]
    ).properties(width=300, height=alt.Step(step))
    pairs = f["clean_pairs"].astype(float)
    lo = 10 ** np.floor(np.log10(max(1.0, min(pairs[pairs > 0].min() if (pairs > 0).any() else 1,
                                              d["min_pairs"] or np.inf) / 1.5)))
    hi = 10 ** np.ceil(np.log10(max(pairs.max() if pairs.notna().any() else 10, 10)))
    zoom = alt.selection_interval(bind="scales", encodings=["x"])
    bars = alt.Chart(f).mark_bar().encode(
        x=alt.X("clean_pairs:Q", scale=alt.Scale(type="log", domain=[lo, hi]), title="clean read pairs (log scale)",
                axis=alt.Axis(format="~s")),
        x2=alt.X2(datum=lo),
        y=alt.Y("sample:O", sort=samples, axis=None),
        color=alt.condition("datum.low_depth", alt.value("#9A9A9A"), alt.value(fs.CATEGORY[0])),
        tooltip=[alt.Tooltip("sample:N")] + tip_g + [alt.Tooltip("clean_pairs:Q", format=","),
                                                     alt.Tooltip("raw_pairs:Q", format=",")]).add_params(zoom)
    layers = [bars]
    if d.get("min_pairs"):
        layers.append(alt.Chart(f).mark_rule(color=fs.LOW_DEPTH, strokeDash=[4, 3], strokeWidth=1)
                      .encode(x=alt.datum(d["min_pairs"])))
    if sel is not None:
        left = left.add_params(sel).transform_filter(sel)
        layers = [x.transform_filter(sel) for x in layers]
    panels = [left, alt.layer(*layers).properties(width=260, height=alt.Step(step))]
    if d.get("classified"):
        cls = alt.Chart(f).mark_bar(color=fs.CATEGORY[2]).encode(
            x=alt.X("pct_classified:Q", title="% classified (profiler)", scale=alt.Scale(domain=[0, 100])),
            y=alt.Y("sample:O", sort=samples, axis=None),
            tooltip=[alt.Tooltip("sample:N")] + tip_g + [alt.Tooltip("pct_classified:Q", format=".1f")])
        if sel is not None:
            cls = cls.transform_filter(sel)
        panels.append(cls.properties(width=180, height=alt.Step(step)))
    return alt.hconcat(*panels, spacing=10)


BUILDERS = {"heatmap": heatmap, "read_fate": read_fate}


def spec_json(key, data):
    """The figure's Vega-Lite spec as JSON text, or None (no interactive version, or above the size budget)."""
    if not data or data.get("kind") not in BUILDERS:
        return None
    import altair as alt
    alt.data_transformers.disable_max_rows()
    spec = BUILDERS[data["kind"]](data).to_dict()
    spec["config"] = config()
    spec["padding"] = {"left": 12, "top": 4, "right": 4, "bottom": 4}
    # inline in a <script>: no '<', '>' or '&' may appear raw (a label with '</script>' would end the block);
    # \u escapes stay valid JSON and JavaScript
    text = (json.dumps(spec, separators=(",", ":"), default=str)
            .replace("<", "\\u003c").replace(">", "\\u003e").replace("&", "\\u0026"))
    if len(text) > SPEC_BUDGET:
        warn(f"{key}: interactive spec {len(text) / 1e6:.1f} MB is above the {SPEC_BUDGET / 1e6:g} MB budget; the "
             "HTML shows the static figure")
        return None
    return text


def javascript_bundle():
    """vega + vega-lite + vega-embed in one script (vl-convert's bundle), for offline use."""
    import vl_convert as vlc
    return vlc.javascript_bundle()
