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
  cog      KO-derived: the ORF's COG ortholog ids (function/kegg/ko-to-cog.map)
           mapped to COG functional-category letters with NCBI's COG
           definitions table (--cog-def, cog-24.def.tab) - the same vocabulary
           as the contig branch (design doc Q14, owner option d; Q16)
  pfam     function/pfam/orf-to-pfam.map.xz, reported as Pfam NAMES via
           pfam_name.txt - the contig branch's vocabulary (eggNOG writes
           names); the versioned accession goes in the description (Q14)
  metacyc  pathways: orf-to-protein.map.xz -> protein-to-enzrxn.map ->
           enzrxn-to-reaction.map -> reaction-to-pathway.map

Descriptions come from ko_name.txt / pathway_name.txt; pfam rows carry the
versioned Pfam accession instead; ec and cog stay empty (no name file).

COG letters are de-duplicated per ORF AFTER mapping, so an ORF whose KOs
reach two COGs of category E counts once toward E (Section 4.8 step 5).
Every id in ko-to-cog must be in the --cog-def table or in
WOLR2_UNMAPPABLE_COGS (known defects of the pinned WoLr2 release, skipped
with a warning); anything else fails loudly, whatever the sample contains.

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
      --cog-def cog-24.def.tab --sample-id S1 --prefix S1
"""

import argparse
import lzma
import re
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

# NCBI COG definitions table (same table and validation as the contig
# branch's parse_cog_def in bin/aggregate_functions.py, which cannot be
# imported here: it needs pandas, which the woltka image lacks)
COG_ID_RE = re.compile(r'^COG\d+$')
COG_LETTERS_RE = re.compile(r'^[A-Z]+$')

# WoLr2 ko-to-cog ids that cog-24.def.tab cannot map (all 3,311 distinct ids
# checked 2026-10-07; design doc Q14, owner: enumerated skip-list). Skipped
# with a warning; the KO's other COGs still count. Any OTHER unmappable id
# fails loudly (a different WoL release or COG table needs re-checking).
WOLR2_UNMAPPABLE_COGS = frozenset([
    # malformed upstream (KO in brackets)
    'COG00028',   # K24393 (its other COG, COG4032, maps)
    'COG:1140',   # K24714
    ':COG1216',   # K25205
    'COG:5013',   # K24713
    'OG3395',     # K23247
    # well-formed but absent from NCBI COG2024
    'COG3632',    # K01571
    'COG3699',    # K07280
    'COG3849',    # K06931
    'COG5273',    # K20032
])


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


def parse_cog_def(path):
    """COG id -> category letters from an NCBI COG definitions table
    (cog-24.def.tab layout, no header). Latin-1 tolerant; whitespace around
    the category is stripped (cog-24 has 'O ' for COG6144)."""
    mapping = {}
    with open(path, encoding='latin-1') as fh:
        for lineno, line in enumerate(fh, start=1):
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 2:
                die(f"{path}:{lineno}: expected the NCBI COG definitions layout "
                    f"(COG id <TAB> category letters <TAB> ...)")
            cog_id, letters = fields[0].strip(), fields[1].strip()
            if not COG_ID_RE.match(cog_id) or not COG_LETTERS_RE.match(letters):
                die(f"{path}:{lineno}: malformed COG definition ({cog_id!r}, {letters!r})")
            mapping[cog_id] = letters
    if not mapping:
        die(f"{path}: empty COG definitions table")
    return mapping


def cog_letters_map(ko_cog, cog_def, cog_def_path):
    """{KO: set(category letters)} plus {KO: set(skipped ids)}. Every COG id
    of the map must be in the table or the enumerated skip-list."""
    letters, skipped = {}, {}
    for ko, cogs in ko_cog.items():
        for cog in cogs:
            if cog in cog_def:
                letters.setdefault(ko, set()).update(cog_def[cog])
            elif cog in WOLR2_UNMAPPABLE_COGS:
                skipped.setdefault(ko, set()).add(cog)
            else:
                die(f"COG id {cog!r} ({ko} in {KO_TO_COG}) is not in the COG "
                    f"definitions table {cog_def_path} and is not a known WoLr2 "
                    f"defect; the WoL release and the COG table must be re-checked "
                    f"together (design doc Q14)")
    return letters, skipped


def pfam_name_map(path):
    """Versioned Pfam accession -> name. Names must be unique, otherwise
    reporting by name would merge distinct families."""
    names = read_names(path)
    seen = {}
    for acc, name in names.items():
        if not name:
            die(f"{path}: Pfam accession {acc} has an empty name")
        if name in seen:
            die(f"{path}: Pfam name {name!r} is shared by {seen[name]} and {acc}; "
                f"pfam rows are reported by name (design doc Q14)")
        seen[name] = acc
    return names


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
    ap.add_argument('--cog-def', required=True,
                    help='NCBI COG definitions table (cog-24.def.tab): COG id -> category letters')
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
    ko_cog_letters, ko_cog_skipped = cog_letters_map(ko_cog, parse_cog_def(args.cog_def),
                                                     args.cog_def)
    orf_pfam = stream_orf_map(require(args.db, ORF_TO_PFAM), orfs)
    orf_protein = stream_orf_map(require(args.db, ORF_TO_PROTEIN), orfs)
    protein_enzrxn = read_map(require(args.db, PROTEIN_TO_ENZRXN))
    enzrxn_reaction = read_map(require(args.db, ENZRXN_TO_REACTION))
    reaction_pathway = read_map(require(args.db, REACTION_TO_PATHWAY))

    pfam_names = pfam_name_map(require(args.db, PFAM_NAMES))
    names = {
        'ko': read_names(require(args.db, KO_NAMES)),
        'metacyc': read_names(require(args.db, PATHWAY_NAMES)),
    }
    # pfam rows: accession = name, description = versioned accession (Q14).
    # Names are unique (checked above), so the reverse lookup is exact
    names['pfam'] = {name: acc for acc, name in pfam_names.items()}

    term_sets = {ont: {} for ont in ONTOLOGIES}
    skipped_ids, skipped_orfs = set(), 0
    for orf in orfs:
        kos = orf_ko.get(orf, set())
        accessions = orf_pfam.get(orf, set())
        unnamed = sorted(a for a in accessions if a not in pfam_names)
        if unnamed:
            die(f"Pfam accession {unnamed[0]} (ORF {orf}) has no entry in {PFAM_NAMES}")
        skipped = expand(kos, ko_cog_skipped)
        if skipped:
            skipped_ids |= skipped
            skipped_orfs += 1
        pathways = expand(expand(expand(orf_protein.get(orf, set()), protein_enzrxn),
                                 enzrxn_reaction), reaction_pathway)
        term_sets['ko'][orf] = kos
        term_sets['ec'][orf] = expand(kos, ko_ec)
        term_sets['cog'][orf] = expand(kos, ko_cog_letters)
        term_sets['pfam'][orf] = {pfam_names[a] for a in accessions}
        term_sets['metacyc'][orf] = pathways

    if skipped_ids:
        print(f"WARNING: {args.sample_id}: {skipped_orfs} ORF(s) carry KOs whose "
              f"ko-to-cog entries cannot be mapped to COG categories and were "
              f"skipped (known WoLr2 defects: {', '.join(sorted(skipped_ids))}); "
              f"their other COGs still count", file=sys.stderr)

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
