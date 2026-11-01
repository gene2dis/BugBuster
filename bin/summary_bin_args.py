"""Bin-level ARG clustering tables of bin/results_summary.py, from CLUSTERING's MMseqs2 tables
(05_arg_prediction/bin_level/clustering/Arg_cluster_<identity>_cluster.tsv: representative, member; no header).

The members are ARG_FASTA_FORMATTER's protein headers, <gene id>_from_<assembly>_<bin>_gene_<ARG>_ARG_class_<class>
(bin/extract_sequences.py). ARG and gene names contain underscores themselves, so the header is split on the
literal separators, and <assembly>_<bin> is resolved against the known assemblies.
"""

import pandas as pd

from summary_common import die

IDENTITIES = (90, 95, 99, 100)


def load_bin_clusters(root):
    """{identity: DataFrame(representative, member)} or None when the clustering did not run."""
    cdir = root / "05_arg_prediction" / "bin_level" / "clustering"
    out = {}
    for ident in IDENTITIES:
        p = cdir / f"Arg_cluster_{ident}_cluster.tsv"
        if p.exists():
            df = pd.read_csv(p, sep="\t", header=None, dtype=str, keep_default_na=False)
            if df.shape[1] != 2:
                die(f"{p}: expected two columns (representative, member), found {df.shape[1]}")
            out[ident] = df.set_axis(["representative", "member"], axis=1)
    return out or None


def parse_header(header, assemblies):
    """'<gene>_from_<assembly>_<bin>_gene_<ARG>_ARG_class_<class>' -> (assembly, bin, gene_id, arg_gene,
    arg_class). assemblies: names to resolve <assembly>_<bin> with (the longest matching prefix wins)."""
    gene_id, sep, rest = header.partition("_from_")
    left, sep2, arg_class = rest.rpartition("_ARG_class_")
    asm_bin, sep3, arg_gene = left.partition("_gene_")
    if not (sep and sep2 and sep3):
        die(f"bin ARG cluster member '{header}' is not <gene>_from_<assembly>_<bin>_gene_<ARG>_ARG_class_<class>")
    matches = [a for a in assemblies if asm_bin.startswith(a + "_")]
    if not matches:
        die(f"bin ARG cluster member '{header}': no known assembly ({sorted(assemblies)}) prefixes '{asm_bin}'")
    asm = max(matches, key=len)
    return asm, asm_bin[len(asm) + 1:], gene_id, arg_gene, arg_class


def arg_mag_long_table(parts):
    """arg_mag_long: one row per ARG gene of a MAG and the cluster it belongs to at each identity (the cluster is
    named after its representative gene, <assembly>:<bin>:<gene id>). parts: [(batch label, clusters,
    assemblies)]."""
    frames = []
    for label, clusters, assemblies in parts:
        if not clusters:
            continue
        names = {}

        def short(header):
            if header not in names:
                a, b, g, _, _ = parse_header(header, assemblies)
                names[header] = f"{a}:{b}:{g}"
            return names[header]
        base_ident = max(clusters)
        rows = []
        for m in clusters[base_ident]["member"]:
            a, b, g, gene, cls = parse_header(m, assemblies)
            rows.append({"assembly": a, "batch": label, "bin": b, "gene_id": g, "arg_gene": gene, "arg_class": cls,
                         "_member": m})
        df = pd.DataFrame(rows)
        for ident in IDENTITIES:
            if ident in clusters:
                rep = dict(zip(clusters[ident]["member"], clusters[ident]["representative"]))
                df[f"cluster_{ident}"] = df["_member"].map(lambda m: short(rep[m]) if m in rep else "")
            else:
                df[f"cluster_{ident}"] = ""
        frames.append(df.drop(columns="_member"))
    if not frames:
        return None
    df = pd.concat(frames, ignore_index=True)
    return df.sort_values(["batch", "assembly", "bin", "arg_class", "arg_gene", "gene_id"],
                          kind="stable").reset_index(drop=True)


def arg_cluster_sharing_table(long, identity=90):
    """arg_cluster_sharing: the clusters at <identity> % found in more than one MAG, per batch (clusters are made
    within one run): the ARG names and classes of their members, how many MAGs and assemblies hold them, and the
    MAGs."""
    if long is None or long.empty:
        return None
    col = f"cluster_{identity}"
    df = long[long[col] != ""].assign(mag=lambda d: d["assembly"] + ":" + d["bin"])
    rows = []
    for (batch, cluster), g in df.groupby(["batch", col], sort=True):
        mags = sorted(set(g["mag"]))
        if len(mags) < 2:
            continue
        rows.append({"cluster": cluster, "batch": batch, "arg_gene": ";".join(sorted(set(g["arg_gene"]))),
                     "arg_class": ";".join(sorted(set(g["arg_class"]))), "n_mags": len(mags),
                     "n_assemblies": g["assembly"].nunique(), "mags": ";".join(mags)})
    out = pd.DataFrame(rows, columns=["cluster", "batch", "arg_gene", "arg_class", "n_mags", "n_assemblies", "mags"])
    return out.sort_values(["batch", "n_mags", "cluster"], ascending=[True, False, True],
                           kind="stable").reset_index(drop=True)
