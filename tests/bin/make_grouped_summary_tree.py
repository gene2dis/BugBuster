#!/usr/bin/env python3
"""make_grouped_summary_tree.py <co-assembly tree> <out tree> [<old>=<new> ...]

A --coassembly_group_by results tree for tests/bin/test_results_summary.sh, made from the co-assembly stage fixture
(tests/data/summary/stages/coassembly: one co-assembly of sampleC and sampleD): two groups with the co-assembly's
files, g1 = sampleC and g2 = sampleD, each under <stage>/coassembly/<group>/ with its files named <group>_*, the
contig table's rows once per group, pipeline_info/coassembly_groups.tsv; the bin-level ARG clusters are left out.
The optional <old>=<new> pairs rename samples in every path and file (a second batch with other sample names).
"""
import shutil
import sys
from pathlib import Path

src, out = Path(sys.argv[1]), Path(sys.argv[2])
renames = [pair.split("=", 1) for pair in sys.argv[3:]]
groups = {"sampleC": "g1", "sampleD": "g2"}


def put(rel, text):
    for old, new in renames:
        rel, text = rel.replace(old, new), text.replace(old, new)
    p = out / rel
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text)


for p in sorted(src.rglob("*")):
    if p.is_dir():
        continue
    rel = str(p.relative_to(src))
    text = p.read_text()
    if "bin_level/clustering" in rel:
        continue
    if rel.startswith("04_binning/coassembly/coverage/"):
        sample = p.name[:-len("_bin_depth.tsv")]
        put(f"04_binning/coassembly/{groups[sample]}/coverage/{p.name}", text)
    elif "/coassembly/" in rel and "coassembly_" in p.name or rel.startswith("03_assembly/coassembly/"):
        for g in sorted(set(groups.values())):
            put(rel.replace("/coassembly/", f"/coassembly/{g}/").replace("coassembly_", f"{g}_"), text)
    elif rel.endswith("Contig_tax_and_arg_prediction.tsv"):
        head, *rows = text.splitlines()
        put(rel, "\n".join([head] + [r.replace("coassembly\t", f"{g}\t", 1) for g in sorted(set(groups.values()))
                                     for r in rows]) + "\n")
    elif rel.endswith("contig_filtering_summary.txt"):
        line = next(x for x in text.splitlines() if x.startswith("coassembly\t"))
        put(rel, text.replace(line, "\n".join(line.replace("coassembly", g, 1) for g in sorted(set(groups.values())))))
    else:
        put(rel, text)
put("pipeline_info/coassembly_groups.tsv", "sample\tgroup\n" + "".join(f"{s}\t{g}\n" for s, g in groups.items()))
