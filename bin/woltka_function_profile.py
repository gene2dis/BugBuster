#!/usr/bin/env python3
"""Turn one sample's Woltka ORF profile into per-function read abundances.

Read-level functional branch, Woltka backend (design doc Section 4.6.2,
task T8a). Runs inside WOLTKA_CLASSIFY, in the pinned woltka container
(stdlib only), right after `woltka classify --coords` produced the ORF-level
profile.

Why the term sets are composed here instead of with `woltka collapse`:
on woltka 0.1.7, collapse sums over every map entry, so a WoLr2 row that
repeats a domain (orf-to-pfam lists e.g. PF03989.16 six times for one ORF)
counts that ORF six times, and chained collapses (ORF -> protein -> enzrxn
-> reaction -> pathway; ORF -> KO -> EC) count an ORF once per path that
reaches the term (verified 2026-10-07). The pipeline's semantics (Section
4.8 step 5, shared with the contig branch) are: a gene contributes its full
abundance ONCE to each distinct term it carries. So each ORF's term set is
built here, de-duplicated, and its count/RPK is added once per term.

Ontologies (WoLr2 maps; Q5 record in Section 4.6.2):
  ko       function/kegg/orf-to-ko.map.xz
  ec       KO-derived: union of function/kegg/ko-to-ec.map over the ORF's KOs
  cog      KO-derived: union of function/kegg/ko-to-cog.map (COG ortholog
           group ids, e.g. COG0604 - NOT the single-letter functional
           categories the contig branch's eggNOG COG_category carries)
  pfam     function/pfam/orf-to-pfam.map.xz (versioned accessions, as emitted)
  metacyc  pathways: orf-to-protein.map.xz -> protein-to-enzrxn.map ->
           enzrxn-to-reaction.map -> reaction-to-pathway.map

Descriptions come from ko_name.txt / pfam_name.txt / pathway_name.txt; ec
and cog have no name file in the release and stay empty.

Abundances per ORF: count = Woltka's (possibly fractional, 1/k-divided)
read count; rpk = count / (length_bp / 1000) with length_bp from
proteins/length.map.xz (Section 6.1's RPK). Profile rows whose ORF has no
length entry are an error (wrong or mismatched database), never skipped.

Outputs:
  <prefix>.woltka_functions.tsv   sample_id ontology accession description count rpk
  <prefix>.woltka_summary.tsv     sample_id ontology reads_assigned
                                  reads_annotated fraction_annotated
                                  ('any' row + one row per ontology;
                                  reads_assigned = reads Woltka assigned to
                                  ORFs, excluding its 'Unassigned' row)
  <prefix>.woltka_unassigned.tsv  sample_id reads_unassigned (Woltka's
                                  'Unassigned' row: reads left unassigned
                                  because of --uniq ambiguity; 0 otherwise)

Usage:
  woltka_function_profile.py --profile orf.tsv --db <wol_db_dir>
      --sample-id S1 --prefix S1
"""

import argparse
import lzma
import sys
from pathlib import Path

ONTOLOGIES = ['ko', 'ec', 'cog', 'pfam', 'metacyc']

FUNCTIONS_COLUMNS = ['sample_id', 'ontology', 'accession', 'description', 'count', 'rpk']
SUMMARY_COLUMNS = ['sample_id', 'ontology', 'reads_assigned', 'reads_annotated',
                   'fraction_annotated']
UNASSIGNED_COLUMNS = ['sample_id', 'reads_unassigned']

# Relative paths inside the WoLr2 layout (Section 4.6.2 Q5 record)
LENGTH_MAP = 'proteins/length.map.xz'
ORF_TO_KO = 'function/kegg/orf-to-ko.map.xz'
KO_TO_EC = 'function/kegg/ko-to-ec.map'
KO_TO_COG = 'function/kegg/ko-to-cog.map'
KO_NAMES = 'function/kegg/ko_name.txt'
ORF_TO_PFAM = 'function/pfam/orf-to-pfam.map.xz'
PFAM_NAMES = 'function/pfam/pfam_name.txt'
ORF_TO_PROTEIN = 'function/metacyc/orf-to-protein.map.xz'
PROTEIN_TO_ENZRXN = 'function/metacyc/protein-to-enzrxn.map'
ENZRXN_TO_REACTION = 'function/metacyc/enzrxn-to-reaction.map'
REACTION_TO_PATHWAY = 'function/metacyc/reaction-to-pathway.map'
PATHWAY_NAMES = 'function/metacyc/pathway_name.txt'

UNASSIGNED = 'Unassigned'


def die(msg):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)


def open_text(path):
    if str(path).endswith('.xz'):
        return lzma.open(path, 'rt')
    return open(path, 'rt')


def require(db, rel):
    path = Path(db) / rel
    if not path.is_file():
        die(f"Woltka database file missing: {path} (expected the WoLr2 layout, "
            f"see --custom_woltka_db in docs/parameters.md)")
    return path


def read_profile(path):
    """Woltka 0.1.7 TSV profile: '#FeatureID<TAB><sample>' header (the
    sample name is empty for a single-file input), then feature<TAB>value.
    Exactly two columns - anything else means a different woltka layout."""
    counts = {}
    unassigned = 0.0
    with open(path) as fh:
        header = fh.readline().rstrip('\n').split('\t')
        if len(header) != 2 or header[0] != '#FeatureID':
            die(f"{path}: unexpected Woltka profile header {header!r}; expected "
                f"['#FeatureID', <sample>] (woltka 0.1.7 TSV layout)")
        for lineno, line in enumerate(fh, start=2):
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 2:
                die(f"{path}:{lineno}: expected 2 columns, got {len(fields)}")
            feature, value = fields
            try:
                value = float(value)
            except ValueError:
                die(f"{path}:{lineno}: non-numeric count {value!r} for {feature}")
            if feature == UNASSIGNED:
                unassigned += value
            elif feature in counts:
                die(f"{path}:{lineno}: duplicate feature {feature}")
            else:
                counts[feature] = value
    return counts, unassigned


def stream_orf_map(path, wanted):
    """{orf: set(targets)} for the ORFs in `wanted` only (the WoLr2 ORF-level
    maps hold tens of millions of rows; a sample touches a small subset)."""
    out = {}
    with open_text(path) as fh:
        for line in fh:
            orf, _, rest = line.rstrip('\n').partition('\t')
            if orf in wanted:
                targets = {t for t in rest.split('\t') if t}
                if targets:
                    out.setdefault(orf, set()).update(targets)
    return out


def read_map(path):
    out = {}
    with open_text(path) as fh:
        for line in fh:
            fields = line.rstrip('\n').split('\t')
            targets = {t for t in fields[1:] if t}
            if fields[0] and targets:
                out.setdefault(fields[0], set()).update(targets)
    return out


def read_names(path):
    out = {}
    with open_text(path) as fh:
        for line in fh:
            key, _, name = line.rstrip('\n').partition('\t')
            if key:
                out[key] = name
    return out


def expand(sources, mapping):
    out = set()
    for s in sources:
        out |= mapping.get(s, set())
    return out


def fmt(value):
    # Full precision in these intermediates: rounding here would leak into
    # the CPGE that AGGREGATE_READ_FUNCTIONS derives (RPK / GE); the final
    # tables are formatted there
    return f"{value:.12g}"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--profile', required=True, help='woltka classify TSV profile (ORF level)')
    ap.add_argument('--db', required=True, help='WoLr2-layout database directory')
    ap.add_argument('--sample-id', required=True)
    ap.add_argument('--prefix', required=True)
    args = ap.parse_args()

    counts, unassigned = read_profile(args.profile)
    orfs = set(counts)

    lengths = {}
    with open_text(require(args.db, LENGTH_MAP)) as fh:
        for line in fh:
            orf, _, length = line.rstrip('\n').partition('\t')
            if orf in orfs:
                lengths[orf] = int(length)
    missing = sorted(orfs - set(lengths))
    if missing:
        die(f"{len(missing)} profile ORF(s) have no entry in {LENGTH_MAP} "
            f"(first: {missing[0]}); the alignment index and the database files "
            f"must come from the same WoL release")
    zero = sorted(o for o in orfs if lengths[o] <= 0)
    if zero:
        die(f"non-positive ORF length in {LENGTH_MAP} for {zero[0]}")

    rpk = {o: counts[o] / (lengths[o] / 1000.0) for o in orfs}

    orf_ko = stream_orf_map(require(args.db, ORF_TO_KO), orfs)
    ko_ec = read_map(require(args.db, KO_TO_EC))
    ko_cog = read_map(require(args.db, KO_TO_COG))
    orf_pfam = stream_orf_map(require(args.db, ORF_TO_PFAM), orfs)
    orf_protein = stream_orf_map(require(args.db, ORF_TO_PROTEIN), orfs)
    protein_enzrxn = read_map(require(args.db, PROTEIN_TO_ENZRXN))
    enzrxn_reaction = read_map(require(args.db, ENZRXN_TO_REACTION))
    reaction_pathway = read_map(require(args.db, REACTION_TO_PATHWAY))

    names = {
        'ko': read_names(require(args.db, KO_NAMES)),
        'pfam': read_names(require(args.db, PFAM_NAMES)),
        'metacyc': read_names(require(args.db, PATHWAY_NAMES)),
    }

    term_sets = {ont: {} for ont in ONTOLOGIES}
    for orf in orfs:
        kos = orf_ko.get(orf, set())
        pathways = expand(expand(expand(orf_protein.get(orf, set()), protein_enzrxn),
                                 enzrxn_reaction), reaction_pathway)
        term_sets['ko'][orf] = kos
        term_sets['ec'][orf] = expand(kos, ko_ec)
        term_sets['cog'][orf] = expand(kos, ko_cog)
        term_sets['pfam'][orf] = orf_pfam.get(orf, set())
        term_sets['metacyc'][orf] = pathways

    total = sum(counts.values())
    function_rows = []
    summary_rows = []
    annotated_any = set()
    for ont in ONTOLOGIES:
        agg = {}
        annotated = 0.0
        for orf in sorted(orfs):
            terms = term_sets[ont][orf]
            if not terms:
                continue
            annotated += counts[orf]
            annotated_any.add(orf)
            for term in terms:
                c, r = agg.get(term, (0.0, 0.0))
                agg[term] = (c + counts[orf], r + rpk[orf])
        for term in sorted(agg):
            c, r = agg[term]
            function_rows.append([args.sample_id, ont, term,
                                  names.get(ont, {}).get(term, ''), fmt(c), fmt(r)])
        summary_rows.append((ont, annotated))
    any_annotated = sum(counts[o] for o in annotated_any)

    def fraction(n):
        return f"{n / total:.6f}" if total > 0 else ''

    with open(f"{args.prefix}.woltka_functions.tsv", 'w') as fh:
        print('\t'.join(FUNCTIONS_COLUMNS), file=fh)
        for row in function_rows:
            print('\t'.join(row), file=fh)

    with open(f"{args.prefix}.woltka_summary.tsv", 'w') as fh:
        print('\t'.join(SUMMARY_COLUMNS), file=fh)
        for ont, annotated in [('any', any_annotated)] + summary_rows:
            print('\t'.join([args.sample_id, ont, fmt(total), fmt(annotated),
                             fraction(annotated)]), file=fh)

    with open(f"{args.prefix}.woltka_unassigned.tsv", 'w') as fh:
        print('\t'.join(UNASSIGNED_COLUMNS), file=fh)
        print('\t'.join([args.sample_id, fmt(unassigned)]), file=fh)


if __name__ == '__main__':
    main()
