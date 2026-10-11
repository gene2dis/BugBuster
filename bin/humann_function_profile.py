#!/usr/bin/env python3
"""Turn one sample's HUMAnN 4.0.0a2 output into per-function abundances.

Read-level functional branch, HUMAnN backend (design doc Section 4.6.1,
task T8c, Q2). Runs inside HUMANN, in the pinned HUMAnN image (stdlib only),
right after `humann --count-normalization RPKs`.

Version-aware parsing (Section 4.6.1, mandatory for the alpha): every table
must carry the exact header of a HUMAnN version listed in KNOWN_HUMANN_VERSIONS,
exactly two columns, non-negative numbers, and only the feature forms verified
on 4.0.0a2 output (2026-10-08). Anything else fails loudly - silent misparsing
is the failure mode this script exists to prevent.

Ontologies (owner decisions at T8c plan review, 2026-10-08):
  metacyc  unstratified MetaCyc pathway abundance (<base>_4_pathabundance.tsv),
           UNMAPPED / UNINTEGRATED excluded; accession = pathway id,
           description = pathway name
  ko       unstratified gene families (<base>_2_genefamilies.tsv) regrouped
           with utility_mapping/map_ko_uniref90.txt.gz (UniRef90 members only:
           UniClust90 families never map to a KO)
  ec       the same, with map_level4ec_uniclust90.txt.gz (UniRef90 and
           UniClust90 members)
The regrouping is done here, not with humann_regroup_table: in 4.0.0a2 that
utility does not protect the READS_UNMAPPED row and folds it into UNGROUPED
(Q2 bug 1). A family carrying several terms of one ontology contributes its
full RPK to each, once (Section 4.8 step 5 semantics).

Abundances are HUMAnN RPKs: count = rpk = the family / pathway RPK (the
aggregator's native unit for this backend is rpk; CPGE = RPK / GE).

Summary rows (decided at T8c plan review):
  any       reads_assigned  = reads given to HUMAnN (--input-reads, counted
                              by the module; mates count separately)
            reads_annotated = reads_assigned - READS_UNMAPPED (READS_UNMAPPED
                              is a raw read count even under RPKs, verified)
  ko/ec/metacyc  reads_* blank; fraction_annotated = share of the
                 unstratified gene-family RPK (excluding READS_UNMAPPED) that
                 carries >= 1 term of the ontology (metacyc: pathway RPK has
                 no family link, so its share is the integrated fraction
                 1 - (UNMAPPED + UNINTEGRATED) / total pathway abundance)

Outputs:
  <prefix>.humann_functions.tsv  sample_id ontology accession description count rpk
  <prefix>.humann_summary.tsv    sample_id ontology reads_assigned reads_annotated
                                 fraction_annotated
Without --genefamilies / --pathabundance (a sample without reads: HUMAnN is
not run) the tables are header-only / zero.

Usage:
  humann_function_profile.py --genefamilies S1_2_genefamilies.tsv
      --pathabundance S1_4_pathabundance.tsv --utility-db <db>/utility_mapping
      --basename S1 --input-reads 12345 --sample-id S1 --prefix S1
"""

import argparse
import gzip
import re
import sys
from pathlib import Path

# HUMAnN versions whose output layout was verified (design doc Section 4.6.1
# Q2 record and T8c implementation record). Extend only after re-verifying the
# layout on real output of the new version.
KNOWN_HUMANN_VERSIONS = {'4.0.0.alpha.2'}

ONTOLOGIES = ['ko', 'ec', 'metacyc']
FUNCTIONS_COLUMNS = ['sample_id', 'ontology', 'accession', 'description', 'count', 'rpk']
SUMMARY_COLUMNS = ['sample_id', 'ontology', 'reads_assigned', 'reads_annotated',
                   'fraction_annotated']

GENEFAMILIES_HEADER = re.compile(r'^# Gene Family HUMAnN v(\S+) RPKs$')
PATHWAY_HEADER = re.compile(r'^# Pathway HUMAnN v(\S+)$')
FAMILY_ID = re.compile(r'^(UniRef90|UniClust90)_\S+$')
PATHWAY_ROW = re.compile(r'^(\S+): (.+)$')
GENEFAMILY_SPECIAL = {'READS_UNMAPPED'}
PATHWAY_SPECIAL = {'UNMAPPED', 'UNINTEGRATED'}

KO_MAP = 'map_ko_uniref90.txt.gz'
EC_MAP = 'map_level4ec_uniclust90.txt.gz'
KO_NAMES = 'map_ko_name.txt.gz'
EC_NAMES = 'map_level4ec_name.txt.gz'


def die(msg):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)


def read_table(path, header_re, expected_column, label):
    """Two-column HUMAnN table: '# <kind> HUMAnN v<version>[ RPKs]<TAB><column>'
    header, then feature<TAB>value rows. Returns the unstratified rows as
    {feature: value} (stratified 'feature|taxon' rows are checked for form,
    not used)."""
    rows = {}
    with open(path) as fh:
        header = fh.readline().rstrip('\n').split('\t')
        if len(header) != 2:
            die(f"{path}: expected a 2-column {label} header, got {header!r}")
        match = header_re.match(header[0])
        if not match:
            die(f"{path}: unexpected {label} header {header[0]!r} - not the HUMAnN "
                f"layout this parser was verified against (expected "
                f"{header_re.pattern!r}; was the run made with --count-normalization RPKs?)")
        version = match.group(1)
        if version not in KNOWN_HUMANN_VERSIONS:
            die(f"{path}: HUMAnN version {version!r} is not among the versions this "
                f"parser was verified against ({sorted(KNOWN_HUMANN_VERSIONS)}). "
                f"Re-verify the output layout on real data of that version and "
                f"extend KNOWN_HUMANN_VERSIONS (design doc Section 4.6.1)")
        if header[1] != expected_column:
            die(f"{path}: {label} sample column {header[1]!r} does not match the "
                f"expected {expected_column!r}")
        for lineno, line in enumerate(fh, start=2):
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 2:
                die(f"{path}:{lineno}: expected 2 columns, got {len(fields)}")
            feature, raw = fields
            try:
                value = float(raw)
            except ValueError:
                die(f"{path}:{lineno}: non-numeric value {raw!r} for {feature!r}")
            if value < 0:
                die(f"{path}:{lineno}: negative value {raw!r} for {feature!r}")
            if '|' in feature:
                continue
            if feature in rows:
                die(f"{path}:{lineno}: duplicate feature {feature!r}")
            rows[feature] = value
    return version, rows


def check_families(path, rows):
    for feature in rows:
        if feature not in GENEFAMILY_SPECIAL and not FAMILY_ID.match(feature):
            die(f"{path}: unexpected gene-family feature {feature!r} (verified forms: "
                f"UniRef90_* / UniClust90_* and {sorted(GENEFAMILY_SPECIAL)})")


def split_pathways(path, rows):
    """{pathway id: (name, value)} plus the special-row values."""
    pathways, special = {}, {}
    for feature, value in rows.items():
        if feature in PATHWAY_SPECIAL:
            special[feature] = value
            continue
        match = PATHWAY_ROW.match(feature)
        if not match:
            die(f"{path}: unexpected pathway feature {feature!r} (verified forms: "
                f"'<id>: <name>' and {sorted(PATHWAY_SPECIAL)})")
        pathway_id, name = match.groups()
        if pathway_id in pathways:
            die(f"{path}: duplicate pathway id {pathway_id!r}")
        pathways[pathway_id] = (name, value)
    return pathways, special


def require(directory, name):
    path = Path(directory) / name
    if not path.is_file() or path.stat().st_size == 0:
        die(f"HUMAnN utility mapping file missing or empty: {path} (expected the "
            f"full_mapping_v4_alpha layout, see --custom_humann_db)")
    return path


def load_map(path, wanted):
    """{family: set(terms)} for the families in `wanted` only (the maps hold
    millions of members). Map lines: term<TAB>member<TAB>member..."""
    out = {}
    with gzip.open(path, 'rt') as fh:
        for line in fh:
            fields = line.rstrip('\n').split('\t')
            term = fields[0]
            for member in fields[1:]:
                if member in wanted:
                    out.setdefault(member, set()).add(term)
    return out


def load_names(path):
    names = {}
    with gzip.open(path, 'rt') as fh:
        for line in fh:
            key, _, name = line.rstrip('\n').partition('\t')
            if key:
                names[key] = name
    return names


def fmt(value):
    # full precision: the aggregator derives CPGE from these values
    return f"{value:.12g}"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--genefamilies', help='<base>_2_genefamilies.tsv (absent: no reads)')
    ap.add_argument('--pathabundance', help='<base>_4_pathabundance.tsv (absent: no reads)')
    ap.add_argument('--utility-db', required=True, help='HUMAnN utility_mapping directory')
    ap.add_argument('--basename', required=True, help='the --output-basename HUMAnN ran with')
    ap.add_argument('--input-reads', required=True, type=int,
                    help='reads given to HUMAnN (mates count separately)')
    ap.add_argument('--sample-id', required=True)
    ap.add_argument('--prefix', required=True)
    args = ap.parse_args()

    if bool(args.genefamilies) != bool(args.pathabundance):
        die("--genefamilies and --pathabundance must be given together")
    if args.input_reads < 0:
        die("--input-reads must be >= 0")

    function_rows = []
    summary = {ont: '' for ont in ONTOLOGIES}
    reads_annotated = 0.0

    if args.genefamilies:
        _, families = read_table(args.genefamilies, GENEFAMILIES_HEADER, args.basename,
                                 'gene-family')
        check_families(args.genefamilies, families)
        _, pathway_rows = read_table(args.pathabundance, PATHWAY_HEADER,
                                     f"{args.basename}_Abundance", 'pathway')
        pathways, pathway_special = split_pathways(args.pathabundance, pathway_rows)

        unmapped_reads = families.pop('READS_UNMAPPED', 0.0)
        if unmapped_reads > args.input_reads + 1e-6:
            die(f"{args.genefamilies}: READS_UNMAPPED ({unmapped_reads:g}) exceeds the "
                f"{args.input_reads} reads given to HUMAnN")
        reads_annotated = args.input_reads - unmapped_reads

        wanted = set(families)
        maps = {
            'ko': load_map(require(args.utility_db, KO_MAP), wanted),
            'ec': load_map(require(args.utility_db, EC_MAP), wanted),
        }
        names = {
            'ko': load_names(require(args.utility_db, KO_NAMES)),
            'ec': load_names(require(args.utility_db, EC_NAMES)),
        }
        family_total = sum(families.values())
        for ont in ('ko', 'ec'):
            agg = {}
            annotated = 0.0
            for family in sorted(families):
                terms = maps[ont].get(family)
                if not terms:
                    continue
                annotated += families[family]
                for term in terms:
                    agg[term] = agg.get(term, 0.0) + families[family]
            for term in sorted(agg):
                function_rows.append([args.sample_id, ont, term, names[ont].get(term, ''),
                                      fmt(agg[term]), fmt(agg[term])])
            summary[ont] = f"{annotated / family_total:.6f}" if family_total > 0 else ''

        for pathway_id in sorted(pathways):
            name, value = pathways[pathway_id]
            function_rows.append([args.sample_id, 'metacyc', pathway_id, name,
                                  fmt(value), fmt(value)])
        integrated = sum(value for _, value in pathways.values())
        pathway_total = integrated + sum(pathway_special.values())
        summary['metacyc'] = f"{integrated / pathway_total:.6f}" if pathway_total > 0 else ''

    with open(f"{args.prefix}.humann_functions.tsv", 'w') as fh:
        print('\t'.join(FUNCTIONS_COLUMNS), file=fh)
        for row in function_rows:
            print('\t'.join(row), file=fh)

    with open(f"{args.prefix}.humann_summary.tsv", 'w') as fh:
        print('\t'.join(SUMMARY_COLUMNS), file=fh)
        any_fraction = (f"{reads_annotated / args.input_reads:.6f}"
                        if args.input_reads > 0 else '')
        print('\t'.join([args.sample_id, 'any', fmt(args.input_reads), fmt(reads_annotated),
                         any_fraction]), file=fh)
        for ont in ONTOLOGIES:
            print('\t'.join([args.sample_id, ont, '', '', summary[ont]]), file=fh)


if __name__ == '__main__':
    main()
