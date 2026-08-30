#!/usr/bin/env python3
"""Aggregate gene counts, eggNOG annotations and gene coordinates into the
canonical functional-annotation tables (design doc Sections 4.8, 5, 6;
tasks T4 (TPM) and T5 (CPGE) — contig branch).

Inputs per sample (per co-assembly for annotations/GFF under coassembly mode):
  - <id>.featureCounts.txt   gene counts (FEATURECOUNTS_GENES)
  - <id>.emapper.annotations eggNOG-mapper v3 annotations
  - <id>.gff[.gz]            Pyrodigal gene coordinates
  - <id>.ags.tsv             MicrobeCensus AGS table (OPTIONAL, per sample:
                             missing samples fall back to TPM-only with empty
                             cpge fields — MicrobeCensus failure is non-fatal)

Outputs (all TSV with header, deterministic ordering):
  - gene_annotations.tsv          Section 5.1 (long, one row per gene per term)
  - gene_abundance.tsv            Section 5.2 (tpm + cpge)
  - function_abundance.tsv        Section 5.3 (source=contigs)
  - function_wide_<ont>_tpm.tsv   wide TPM matrix per ontology (ko, cog, ec, pfam, cazy)
  - function_wide_<ont>_cpge.tsv  wide CPGE matrix per ontology (blank cells for
                                  samples without AGS)
  - annotated_fraction.tsv        per-sample annotated fraction by count and abundance
  - ags_and_ge.tsv                per-sample AGS/genome-equivalents summary with an
                                  'unavailable' status row for AGS-less samples

Version-aware parsing (design doc 4.2/4.6.1, binding): the eggNOG-mapper
version is read from the EGGNOG_MAPPER_ANNOTATE versions.yml and must be a
version whose output layout this parser was verified against; the annotations
column header must match that layout exactly. Any mismatch is a hard error —
silent misparsing is the failure mode this guards against.

Term splitting (Section 4.8 step 4, corrected 2026-08-29): KEGG_ko, EC, PFAMs
and CAZy are comma-joined; COG_category is an undelimited letter string and
splits per character. A gene carrying several terms of one ontology
contributes its full abundance to each (intentional double-counting).
"""

import argparse
import gzip
import re
import sys
from pathlib import Path

import pandas as pd

# eggNOG-mapper versions whose annotations layout this parser is verified
# against (design doc Section 4.2 implementation record). Extend only after
# re-verifying the column list on real output of the new version.
KNOWN_EMAPPER_VERSIONS = {'3.0.0-beta6'}

# The exact v3.0.0-beta6 column header (22 columns; '#' is glued to 'query').
# --md5 would append a 23rd column; the pipeline never passes it.
EXPECTED_ANNOTATION_COLUMNS = (
    '#query', 'seed_ortholog', 'evalue', 'score', 'eggNOG_OGs', 'tax_ceiling',
    'farthest_donor_lineage', 'COG_category', 'Preferred_name', 'GOs', 'EC',
    'KEGG_ko', 'KEGG_Pathway', 'KEGG_Module', 'KEGG_Reaction', 'KEGG_rclass',
    'BRITE', 'KEGG_TC', 'CAZy', 'BiGG_Reaction', 'PFAMs',
    'annotation_confidence',
)

# (annotations column, ontology name, 5.1 db value, split mode)
ONTOLOGY_FIELDS = [
    ('KEGG_ko', 'ko', 'eggnog_ko', 'comma'),
    ('COG_category', 'cog', 'eggnog_cog', 'chars'),
    ('EC', 'ec', 'eggnog_ec', 'comma'),
    ('PFAMs', 'pfam', 'eggnog_pfam', 'comma'),
    ('CAZy', 'cazy', 'eggnog_cazy', 'comma'),
]

ONTOLOGIES = [ontology for _, ontology, _, _ in ONTOLOGY_FIELDS]

FEATURECOUNTS_HEADER_PREFIX = ['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length']

AGS_COLUMNS = ['sample_id', 'average_genome_size_bp', 'genome_equivalents',
               'total_bases']

TPM_FORMAT = '{:.4f}'
CPGE_FORMAT = '{:.6f}'
FRACTION_FORMAT = '{:.4f}'


def fail(message):
    print(f"Error: {message}", file=sys.stderr)
    sys.exit(1)


def open_text(path):
    """Open plain or gzipped text transparently (Pyrodigal GFFs are .gz)."""
    if str(path).endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path, 'r')


def sample_id_from(path, suffixes):
    name = Path(path).name
    for suffix in suffixes:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    fail(f"cannot derive sample id from '{name}' (expected suffix {suffixes})")


def parse_versions_yml(path):
    """Extract eggnog-mapper and eggnog_db versions; enforce the known-layout set."""
    text = Path(path).read_text()
    tool_match = re.search(r'^\s*eggnog-mapper:\s*(\S+)\s*$', text, re.MULTILINE)
    db_match = re.search(r'^\s*eggnog_db:\s*(\S+)\s*$', text, re.MULTILINE)
    if not tool_match:
        fail(f"no 'eggnog-mapper:' version found in {path} — cannot verify the "
             f"annotations layout (design doc 4.2 requires version-aware parsing)")
    tool_version = tool_match.group(1)
    if tool_version == 'unknown' or tool_version not in KNOWN_EMAPPER_VERSIONS:
        fail(f"eggnog-mapper version '{tool_version}' is not among the layouts this "
             f"parser was verified against ({sorted(KNOWN_EMAPPER_VERSIONS)}). "
             f"Re-verify the annotations column layout on real output of that "
             f"version and extend KNOWN_EMAPPER_VERSIONS before aggregating.")
    db_version = db_match.group(1) if db_match else 'unknown'
    return tool_version, db_version


def parse_gff(path):
    """Pyrodigal/Prodigal GFF3 -> one row per CDS with the FAA-matching gene id.

    Prodigal's ID attribute is <seqnum>_<genenum>; the gene id used everywhere
    downstream is <contig>_<genenum> (same rewrite as bin/gff2saf.sh).
    """
    rows = []
    with open_text(path) as handle:
        for line_no, line in enumerate(handle, 1):
            line = line.rstrip('\n')
            if not line or line.startswith('#'):
                continue
            fields = line.split('\t')
            if len(fields) < 9 or fields[2] != 'CDS':
                continue
            attrs = dict(
                part.split('=', 1)
                for part in fields[8].split(';')
                if '=' in part
            )
            raw_id = attrs.get('ID', '')
            genenum = raw_id.rsplit('_', 1)[-1] if raw_id else ''
            if not genenum.isdigit():
                fail(f"{path}: missing or malformed ID attribute on line {line_no}: "
                     f"'{fields[8]}'")
            rows.append({
                'gene_id': f"{fields[0]}_{genenum}",
                'contig_id': fields[0],
                'start': int(fields[3]),
                'end': int(fields[4]),
                'strand': fields[6],
                'partial': attrs.get('partial', ''),
            })
    genes = pd.DataFrame(rows, columns=['gene_id', 'contig_id', 'start', 'end',
                                        'strand', 'partial'])
    if not genes.empty:
        genes['length_bp'] = genes['end'] - genes['start'] + 1
        duplicates = genes['gene_id'][genes['gene_id'].duplicated()]
        if not duplicates.empty:
            fail(f"{path}: duplicate gene ids after rewrite: "
                 f"{sorted(duplicates.unique())[:5]}")
    else:
        genes['length_bp'] = pd.Series(dtype=int)
    return genes


def parse_ags(path):
    """MICROBECENSUS <id>.ags.tsv -> (sample, values dict).

    Layout guard in the house style: exact header, exactly one data row, the
    embedded sample_id must equal the filename-derived one, and every value
    must be a positive number (the module already sanity-checked them; this
    re-check catches hand-fed or mis-paired files). The dict keeps the raw
    strings for output fidelity plus '<name>_float' parsed values.
    """
    sample = sample_id_from(path, ['.ags.tsv'])
    with open_text(path) as handle:
        lines = [line.rstrip('\n') for line in handle if line.strip()]
    if not lines or lines[0].split('\t') != AGS_COLUMNS:
        fail(f"{path}: unexpected ags.tsv header (expected {AGS_COLUMNS})")
    if len(lines) != 2:
        fail(f"{path}: expected exactly one data row, found {len(lines) - 1}")
    fields = lines[1].split('\t')
    if len(fields) != len(AGS_COLUMNS):
        fail(f"{path}: data row has {len(fields)} fields, expected "
             f"{len(AGS_COLUMNS)}")
    if fields[0] != sample:
        fail(f"{path}: embedded sample_id '{fields[0]}' does not match the "
             f"filename-derived id '{sample}'")
    values = {}
    for name, raw in zip(AGS_COLUMNS[1:], fields[1:]):
        try:
            value = float(raw)
        except ValueError:
            value = None
        if value is None or value <= 0:
            fail(f"{path}: {name} must be a positive number, got '{raw}'")
        values[name] = raw
        values[name + '_float'] = value
    return sample, values


def parse_counts(path):
    """featureCounts table -> one row per gene (header-only tables are valid)."""
    header = None
    rows = []
    with open_text(path) as handle:
        for line in handle:
            line = line.rstrip('\n')
            if not line or line.startswith('#'):
                continue
            fields = line.split('\t')
            if header is None:
                if fields[:6] != FEATURECOUNTS_HEADER_PREFIX or len(fields) < 7:
                    fail(f"{path}: unexpected featureCounts header: {fields[:7]} — "
                         f"expected {FEATURECOUNTS_HEADER_PREFIX} + <bam>")
                header = fields
                continue
            if len(fields) != len(header):
                fail(f"{path}: row with {len(fields)} fields does not match the "
                     f"{len(header)}-column header")
            rows.append({
                'gene_id': fields[0],
                'fc_length': int(fields[5]),
                'count': int(fields[6]),
            })
    if header is None:
        fail(f"{path}: no featureCounts header line found")
    return pd.DataFrame(rows, columns=['gene_id', 'fc_length', 'count'])


def parse_annotations(path, tool_version):
    """eggNOG-mapper annotations -> one row per annotated gene.

    Layout guard: the '#query' header line must match the verified column list
    for the pinned version exactly. '##' comment lines (variable-length leading
    block, trailing summary) are skipped wherever they occur; a '## emapper-'
    line, when present, must agree with versions.yml. All comment lines are
    absent under --no_file_comments, so only the '#query' line is required.
    """
    header_seen = False
    rows = []
    with open_text(path) as handle:
        for line_no, line in enumerate(handle, 1):
            line = line.rstrip('\n')
            if not line:
                continue
            if line.startswith('##'):
                version_match = re.match(r'^## emapper-(\S+)', line)
                if version_match and version_match.group(1) != tool_version:
                    fail(f"{path}: header says emapper-{version_match.group(1)} but "
                         f"versions.yml says {tool_version} — mismatched inputs")
                continue
            if line.startswith('#query'):
                columns = tuple(line.split('\t'))
                if columns != EXPECTED_ANNOTATION_COLUMNS:
                    fail(f"{path}: annotations column header does not match the "
                         f"verified emapper-{tool_version} layout "
                         f"({len(EXPECTED_ANNOTATION_COLUMNS)} columns). Got "
                         f"{len(columns)} columns: {list(columns)[:8]}... — refusing "
                         f"to guess (design doc 4.2: fail loudly on layout drift)")
                header_seen = True
                continue
            if line.startswith('#'):
                fail(f"{path}: unrecognized comment line {line_no}: '{line[:60]}'")
            if not header_seen:
                fail(f"{path}: data row before the '#query' column header "
                     f"(line {line_no}) — not a valid emapper annotations file")
            fields = line.split('\t')
            if len(fields) != len(EXPECTED_ANNOTATION_COLUMNS):
                fail(f"{path}: line {line_no} has {len(fields)} fields, expected "
                     f"{len(EXPECTED_ANNOTATION_COLUMNS)}")
            record = dict(zip(EXPECTED_ANNOTATION_COLUMNS, fields))
            rows.append({
                'gene_id': record['#query'],
                'evalue': record['evalue'],
                'score': record['score'],
                **{column: record[column] for column, _, _, _ in ONTOLOGY_FIELDS},
            })
    if not header_seen:
        fail(f"{path}: no '#query' column header line found — not a valid "
             f"emapper annotations file")
    columns = ['gene_id', 'evalue', 'score'] + [c for c, _, _, _ in ONTOLOGY_FIELDS]
    return pd.DataFrame(rows, columns=columns)


def split_terms(value, mode):
    """Split one annotation field into terms ('-' and empty mean unannotated)."""
    if value in ('', '-'):
        return []
    if mode == 'chars':
        return [ch for ch in value if ch not in ('-', ' ')]
    return [term for term in value.split(',') if term and term != '-']


def explode_annotations(annotations):
    """Long form: one row per gene per (db, ontology, accession)."""
    rows = []
    for record in annotations.itertuples(index=False):
        for column, ontology, db, mode in ONTOLOGY_FIELDS:
            for accession in split_terms(getattr(record, column), mode):
                rows.append({
                    'gene_id': record.gene_id,
                    'db': db,
                    'ontology': ontology,
                    'accession': accession,
                    'evalue': record.evalue,
                    'score': record.score,
                })
    return pd.DataFrame(rows, columns=['gene_id', 'db', 'ontology', 'accession',
                                       'evalue', 'score'])


def compute_rpk(abundance):
    """RPK_i = count_i / (length_i / 1000), the shared numerator of 6.1/6.2."""
    return abundance['count'] / (abundance['length_bp'] / 1000.0)


def compute_tpm(abundance, rpk):
    """Section 6.1: TPM_i = RPK_i / sum(RPK) * 1e6, per sample.

    A sample whose genes attracted zero reads gets tpm 0.0 for every gene
    (never NaN); its TPM column then sums to 0, not 1e6.
    """
    totals = rpk.groupby(abundance['sample_id']).transform('sum')
    tpm = (rpk / totals * 1e6).where(totals > 0, 0.0)
    return tpm.fillna(0.0)


def compute_cpge(abundance, rpk, ge_by_sample):
    """Section 6.2: CPGE_i = RPK_i / genome_equivalents, per sample.

    Samples without a MicrobeCensus genome-equivalents value map to NaN
    (rendered as empty fields downstream — the TPM-only fallback).
    """
    genome_equivalents = abundance['sample_id'].map(ge_by_sample)
    return rpk / genome_equivalents


def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Aggregate gene counts, eggNOG annotations, gene '
                    'coordinates and MicrobeCensus AGS estimates into '
                    'functional abundance tables (TPM and CPGE)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Per-sample assembly mode (matching id sets across the three lists;
  # --ags is optional and may cover only a subset of samples)
  aggregate_functions.py --assembly-mode assembly \\
      --counts s1.featureCounts.txt s2.featureCounts.txt \\
      --annotations s1.emapper.annotations s2.emapper.annotations \\
      --gffs s1.gff.gz s2.gff.gz \\
      --ags s1.ags.tsv s2.ags.tsv \\
      --eggnog-versions-yml versions.yml --output-dir .

  # Co-assembly mode (one shared annotations file and GFF)
  aggregate_functions.py --assembly-mode coassembly \\
      --counts s1.featureCounts.txt s2.featureCounts.txt \\
      --annotations coassembly.emapper.annotations \\
      --gffs coassembly.gff.gz \\
      --eggnog-versions-yml versions.yml --output-dir .
        """
    )
    parser.add_argument('--counts', required=True, nargs='+', type=Path,
                        help='featureCounts tables, one per sample')
    parser.add_argument('--annotations', required=True, nargs='+', type=Path,
                        help='eggNOG-mapper annotations (per sample, or one '
                             'shared file under coassembly)')
    parser.add_argument('--gffs', required=True, nargs='+', type=Path,
                        help='Pyrodigal gene GFFs, plain or gzipped (per '
                             'sample, or one shared file under coassembly)')
    parser.add_argument('--ags', nargs='*', type=Path, default=[],
                        help='MicrobeCensus <id>.ags.tsv tables (optional; a '
                             'subset of samples or none — missing samples fall '
                             'back to TPM-only with empty cpge fields)')
    parser.add_argument('--eggnog-versions-yml', required=True, type=Path,
                        help='versions.yml from EGGNOG_MAPPER_ANNOTATE (source '
                             'of the pinned tool and database versions)')
    parser.add_argument('--assembly-mode', required=True,
                        choices=['assembly', 'coassembly'],
                        help='Drives the join topology and the assembly_mode column')
    parser.add_argument('--output-dir', required=True, type=Path,
                        help='Directory for the output tables')
    return parser.parse_args()


def main():
    args = parse_arguments()

    for path in [*args.counts, *args.annotations, *args.gffs, *args.ags,
                 args.eggnog_versions_yml]:
        if not path.exists():
            fail(f"input file not found: {path}")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    tool_version, db_version = parse_versions_yml(args.eggnog_versions_yml)

    counts_by_sample = {
        sample_id_from(path, ['.featureCounts.txt']): parse_counts(path)
        for path in args.counts
    }
    annotations_by_id = {
        sample_id_from(path, ['.emapper.annotations']):
            parse_annotations(path, tool_version)
        for path in args.annotations
    }
    gffs_by_id = {
        sample_id_from(path, ['.gff.gz', '.gff']): parse_gff(path)
        for path in args.gffs
    }
    ags_by_sample = {}
    for path in args.ags:
        sample, values = parse_ags(path)
        if sample in ags_by_sample:
            fail(f"duplicate --ags table for sample '{sample}'")
        ags_by_sample[sample] = values

    sample_ids = sorted(counts_by_sample)

    unknown_ags = sorted(set(ags_by_sample) - set(counts_by_sample))
    if unknown_ags:
        fail(f"--ags sample ids not present in --counts: {unknown_ags}")
    ge_by_sample = {sample: values['genome_equivalents_float']
                    for sample, values in ags_by_sample.items()}
    ags_missing = sorted(set(sample_ids) - set(ags_by_sample))
    if ags_missing:
        print(f"Warning: no MicrobeCensus AGS for sample(s) "
              f"{', '.join(ags_missing)} — cpge left empty (TPM-only "
              f"fallback, design doc Section 4.7)", file=sys.stderr)

    # Join topology (design doc Section 3.5): per-sample gene sets under
    # 'assembly'; one shared co-assembly gene set under 'coassembly'
    if args.assembly_mode == 'assembly':
        for label, keys in [('annotations', annotations_by_id),
                            ('GFFs', gffs_by_id)]:
            if set(keys) != set(counts_by_sample):
                fail(f"sample ids of --counts {sorted(counts_by_sample)} and "
                     f"--{label.lower()} {sorted(keys)} do not match")
        gff_for = dict(gffs_by_id)
        annotations_for = dict(annotations_by_id)
    else:
        if len(annotations_by_id) != 1 or len(gffs_by_id) != 1:
            fail(f"coassembly mode expects exactly one annotations file and one "
                 f"GFF, got {len(annotations_by_id)} and {len(gffs_by_id)}")
        shared_gff = next(iter(gffs_by_id.values()))
        shared_annotations = next(iter(annotations_by_id.values()))
        gff_for = {sample: shared_gff for sample in sample_ids}
        annotations_for = {sample: shared_annotations for sample in sample_ids}

    # Cross-checks: mis-paired inputs must fail, not silently misjoin
    for sample in sample_ids:
        counted = counts_by_sample[sample]
        genes = gff_for[sample]
        counted_ids = set(counted['gene_id'])
        gff_ids = set(genes['gene_id'])
        if counted_ids != gff_ids:
            missing = sorted(gff_ids - counted_ids)[:5]
            extra = sorted(counted_ids - gff_ids)[:5]
            fail(f"sample '{sample}': counted gene ids do not match the GFF gene "
                 f"set (missing from counts: {missing}; not in GFF: {extra})")
        merged = counted.merge(genes[['gene_id', 'length_bp']], on='gene_id')
        bad_length = merged[merged['fc_length'] != merged['length_bp']]
        if not bad_length.empty:
            fail(f"sample '{sample}': featureCounts Length disagrees with the GFF "
                 f"for genes {sorted(bad_length['gene_id'])[:5]}")
        stray = set(annotations_for[sample]['gene_id']) - gff_ids
        if stray:
            fail(f"sample '{sample}': annotation query ids not present in the GFF: "
                 f"{sorted(stray)[:5]}")

    # --- Section 5.2: gene abundance with TPM and CPGE ---
    abundance_parts = []
    for sample in sample_ids:
        part = counts_by_sample[sample].merge(
            gff_for[sample][['gene_id', 'length_bp']], on='gene_id')
        part.insert(0, 'sample_id', sample)
        abundance_parts.append(part)
    abundance = (pd.concat(abundance_parts, ignore_index=True)
                 if abundance_parts else
                 pd.DataFrame(columns=['sample_id', 'gene_id', 'fc_length',
                                       'count', 'length_bp']))
    if not abundance.empty:
        rpk = compute_rpk(abundance)
        abundance['tpm'] = compute_tpm(abundance, rpk)
        abundance['cpge_num'] = compute_cpge(abundance, rpk, ge_by_sample)
    else:
        abundance['tpm'] = []
        abundance['cpge_num'] = []
    abundance['assembly_mode'] = args.assembly_mode
    # NaN cpge_num = sample without AGS -> empty field (TPM-only fallback)
    abundance['cpge'] = abundance['cpge_num'].map(
        lambda value: '' if pd.isna(value) else CPGE_FORMAT.format(value))
    abundance = abundance.sort_values(['sample_id', 'gene_id'],
                                      kind='mergesort', ignore_index=True)

    gene_abundance = abundance[['sample_id', 'gene_id', 'assembly_mode', 'count',
                                'length_bp', 'tpm', 'cpge']].copy()
    gene_abundance['tpm'] = gene_abundance['tpm'].map(TPM_FORMAT.format)
    gene_abundance.to_csv(args.output_dir / 'gene_abundance.tsv', sep='\t',
                          index=False)

    # --- Section 5.1: gene annotations (long, per annotations unit) ---
    exploded_by_id = {
        unit: explode_annotations(table)
        for unit, table in annotations_by_id.items()
    }
    annotation_parts = []
    for unit in sorted(exploded_by_id):
        exploded = exploded_by_id[unit]
        if exploded.empty:
            continue
        # Under coassembly the single annotations unit pairs with the single
        # shared GFF whatever either file is named; under assembly the ids match
        genes = next(iter(gffs_by_id.values())) \
            if args.assembly_mode == 'coassembly' else gff_for[unit]
        part = exploded.merge(
            genes[['gene_id', 'contig_id', 'start', 'end', 'strand', 'partial',
                   'length_bp']],
            on='gene_id')
        part.insert(0, 'sample_id', unit)
        annotation_parts.append(part)
    if annotation_parts:
        gene_annotations = pd.concat(annotation_parts, ignore_index=True)
    else:
        gene_annotations = pd.DataFrame(columns=['sample_id', 'gene_id', 'db',
                                                 'ontology', 'accession', 'evalue',
                                                 'score', 'contig_id', 'start',
                                                 'end', 'strand', 'partial',
                                                 'length_bp'])
    gene_annotations['description'] = ''  # v3 dropped Description (design doc 4.8 record)
    gene_annotations['tool'] = 'eggnog-mapper'
    gene_annotations['tool_version'] = tool_version
    gene_annotations['db_version'] = db_version
    gene_annotations = gene_annotations.sort_values(
        ['sample_id', 'gene_id', 'db', 'accession'], kind='mergesort',
        ignore_index=True)
    gene_annotations[['sample_id', 'gene_id', 'contig_id', 'start', 'end',
                      'strand', 'length_bp', 'partial', 'db', 'accession',
                      'description', 'evalue', 'score', 'tool', 'tool_version',
                      'db_version']].to_csv(
        args.output_dir / 'gene_annotations.tsv', sep='\t', index=False)

    # --- Section 5.3: function abundance (per-sample term explode x TPM/CPGE) ---
    term_parts = []
    for sample in sample_ids:
        exploded = explode_annotations(annotations_for[sample])
        if exploded.empty:
            continue
        sample_abundance = abundance.loc[abundance['sample_id'] == sample,
                                         ['gene_id', 'tpm', 'cpge_num']]
        part = exploded[['gene_id', 'ontology', 'accession']].merge(
            sample_abundance, on='gene_id')
        part.insert(0, 'sample_id', sample)
        term_parts.append(part)
    if term_parts:
        terms = pd.concat(term_parts, ignore_index=True)
        # Intentional double-counting (Section 4.8 step 5): a gene carrying
        # several terms of one ontology contributes its full TPM/CPGE to each
        function_abundance = (
            terms.groupby(['sample_id', 'ontology', 'accession'],
                          as_index=False)[['tpm', 'cpge_num']].sum()
            .rename(columns={'tpm': 'abundance_tpm'}))
    else:
        function_abundance = pd.DataFrame(columns=['sample_id', 'ontology',
                                                   'accession', 'abundance_tpm',
                                                   'cpge_num'])
    function_abundance['source'] = 'contigs'
    function_abundance['backend'] = 'eggnog-mapper'
    function_abundance['description'] = ''
    # groupby.sum() turns an all-NaN group into 0.0, so AGS availability (a
    # per-sample fact) decides emptiness, not the summed value
    function_abundance['abundance_cpge'] = [
        CPGE_FORMAT.format(value) if sample in ge_by_sample else ''
        for sample, value in zip(function_abundance['sample_id'],
                                 function_abundance['cpge_num'])
    ]
    function_abundance['abundance_native'] = ''  # read branch only
    function_abundance['native_unit'] = ''       # read branch only
    function_abundance = function_abundance.sort_values(
        ['sample_id', 'ontology', 'accession'], kind='mergesort',
        ignore_index=True)
    out_53 = function_abundance[['sample_id', 'source', 'backend', 'ontology',
                                 'accession', 'description', 'abundance_tpm',
                                 'abundance_cpge', 'abundance_native',
                                 'native_unit']].copy()
    out_53['abundance_tpm'] = out_53['abundance_tpm'].map(TPM_FORMAT.format)
    out_53.to_csv(args.output_dir / 'function_abundance.tsv', sep='\t',
                  index=False)

    # --- Wide matrices per ontology (every sample a column, missing -> 0;
    #     CPGE cells are blank for samples without AGS) ---
    ags_samples = [sample for sample in sample_ids if sample in ge_by_sample]
    for ontology in ONTOLOGIES:
        subset = function_abundance[function_abundance['ontology'] == ontology]
        wide = subset.pivot_table(index='accession', columns='sample_id',
                                  values='abundance_tpm', aggfunc='sum',
                                  fill_value=0.0)
        wide = wide.reindex(columns=sample_ids, fill_value=0.0).sort_index()
        wide.index.name = 'accession'
        # float_format instead of per-cell mapping: DataFrame.map needs
        # pandas >= 2.1 and the pinned image ships 2.0
        wide.to_csv(args.output_dir / f'function_wide_{ontology}_tpm.tsv',
                    sep='\t', header=True, float_format='%.4f')

        # CPGE mirror: same accession rows as the TPM matrix; AGS-less sample
        # columns are all-blank (never 0.0 — absence of a value, not a zero)
        cpge_subset = subset[subset['sample_id'].isin(ags_samples)]
        wide_cpge = cpge_subset.pivot_table(index='accession',
                                            columns='sample_id',
                                            values='cpge_num', aggfunc='sum',
                                            fill_value=0.0)
        wide_cpge = wide_cpge.reindex(index=wide.index, columns=sample_ids)
        if ags_samples:
            wide_cpge[ags_samples] = wide_cpge[ags_samples].fillna(0.0)
        wide_cpge.index.name = 'accession'
        # Cell-wise formatting: pandas 2.0 (pinned image) drops float_format
        # when na_rep is also given, and DataFrame.map needs >= 2.1
        for column in wide_cpge.columns:
            wide_cpge[column] = ['' if pd.isna(value) else CPGE_FORMAT.format(value)
                                 for value in wide_cpge[column]]
        wide_cpge.to_csv(args.output_dir / f'function_wide_{ontology}_cpge.tsv',
                         sep='\t', header=True)

    # --- Annotated fraction per sample, by count and by abundance ---
    fraction_rows = []
    for sample in sample_ids:
        sample_genes = abundance[abundance['sample_id'] == sample]
        genes_total = len(sample_genes)
        tpm_total = sample_genes['tpm'].sum()
        exploded = explode_annotations(annotations_for[sample])
        for ontology in ['any'] + ONTOLOGIES:
            if exploded.empty:
                annotated_ids = set()
            elif ontology == 'any':
                annotated_ids = set(exploded['gene_id'])
            else:
                annotated_ids = set(
                    exploded.loc[exploded['ontology'] == ontology, 'gene_id'])
            annotated_ids &= set(sample_genes['gene_id'])
            genes_annotated = len(annotated_ids)
            if genes_total > 0:
                by_count = FRACTION_FORMAT.format(genes_annotated / genes_total)
            else:
                by_count = ''
            if tpm_total > 0:
                annotated_tpm = sample_genes.loc[
                    sample_genes['gene_id'].isin(annotated_ids), 'tpm'].sum()
                by_abundance = FRACTION_FORMAT.format(annotated_tpm / tpm_total)
            else:
                by_abundance = ''
            fraction_rows.append({
                'sample_id': sample,
                'ontology': ontology,
                'genes_total': genes_total,
                'genes_annotated': genes_annotated,
                'fraction_by_count': by_count,
                'fraction_by_abundance': by_abundance,
            })
    pd.DataFrame(fraction_rows, columns=['sample_id', 'ontology', 'genes_total',
                                         'genes_annotated', 'fraction_by_count',
                                         'fraction_by_abundance']).to_csv(
        args.output_dir / 'annotated_fraction.tsv', sep='\t', index=False)

    # --- AGS / genome equivalents summary (Section 9.1 "MicrobeCensus AGS
    #     table"); the 'unavailable' rows are the recorded warning that
    #     Section 4.7's non-fatal-failure fallback asks for ---
    ags_rows = []
    for sample in sample_ids:
        values = ags_by_sample.get(sample)
        ags_rows.append({
            'sample_id': sample,
            'average_genome_size_bp': values['average_genome_size_bp'] if values else '',
            'genome_equivalents': values['genome_equivalents'] if values else '',
            'total_bases': values['total_bases'] if values else '',
            'status': 'ok' if values else 'unavailable',
        })
    pd.DataFrame(ags_rows, columns=AGS_COLUMNS + ['status']).to_csv(
        args.output_dir / 'ags_and_ge.tsv', sep='\t', index=False)

    annotated_any = sum(1 for row in fraction_rows
                        if row['ontology'] == 'any' and row['genes_annotated'])
    print(f"Aggregated {len(sample_ids)} sample(s), "
          f"{len(gene_abundance)} gene abundance rows, "
          f"{len(out_53)} function abundance rows "
          f"({annotated_any}/{len(sample_ids)} samples with annotations; "
          f"{len(ags_samples)}/{len(sample_ids)} samples with CPGE)")


if __name__ == '__main__':
    main()
