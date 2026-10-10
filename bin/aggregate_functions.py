#!/usr/bin/env python3
"""Aggregate gene counts, eggNOG annotations, dbCAN CAZy calls and gene
coordinates into the canonical functional-annotation tables (design doc
Sections 4.8, 5, 6; tasks T4 (TPM), T5 (CPGE) and T6 (dbCAN) — contig branch).

Inputs per sample (per co-assembly for annotations/GFF/dbCAN under coassembly
mode):
  - <id>.featureCounts.txt   gene counts (FEATURECOUNTS_GENES)
  - <id>.emapper.annotations eggNOG-mapper v3 annotations
  - <id>.gff[.gz]            Pyrodigal gene coordinates
  - <id>.overview.tsv        run_dbcan v5 overview (OPTIONAL as a set: absent
                             entirely when --functional_cazy false)
  - <id>.ags.tsv             MicrobeCensus AGS table (OPTIONAL, per sample:
                             missing samples fall back to TPM-only with empty
                             cpge fields — MicrobeCensus failure is non-fatal)

Outputs (all TSV with header, deterministic ordering):
  - gene_annotations.tsv          Section 5.1 (long, one row per gene per term;
                                  db=eggnog_* and db=dbcan rows, per-row tool
                                  provenance)
  - gene_abundance.tsv            Section 5.2 (tpm + cpge)
  - function_abundance.tsv        Section 5.3 (source=contigs; backend
                                  distinguishes eggnog-mapper from run_dbcan)
  - function_wide_<ont>_tpm.tsv   wide TPM matrix per ontology (ko, cog, ec,
                                  pfam, cazy — eggNOG-derived) plus cazy_dbcan
                                  (run_dbcan-derived; the two CAZy backends are
                                  never merged into one matrix)
  - function_wide_<ont>_cpge.tsv  wide CPGE matrix per ontology (blank cells for
                                  samples without AGS)
  - annotated_fraction.tsv        per-sample annotated fraction by count and
                                  abundance ('cazy' counts either backend;
                                  'cazy_dbcan' counts run_dbcan alone)
  - ags_and_ge.tsv                per-sample AGS/genome-equivalents summary with an
                                  'unavailable' status row for AGS-less samples

Version-aware parsing (design doc 4.2/4.6.1, binding): the eggNOG-mapper and
run_dbcan versions are read from their modules' versions.yml files and must be
versions whose output layouts this parser was verified against; the column
headers must match those layouts exactly. Any mismatch is a hard error —
silent misparsing is the failure mode this guards against.

Term splitting (Section 4.8 step 4, corrected 2026-08-29): KEGG_ko, EC, PFAMs
and CAZy are comma-joined; COG_category is an undelimited letter string and
splits per character - and when it holds a COG id instead (eggNOG 7 does for
most genes), the id is mapped to its category letters via --cog-def (NCBI
cog-24.def.tab; design doc Q16). PFAMs values carry '_<start>_<end>' domain coordinates
in v3, stripped to the Pfam name (design doc Q15, 2026-10-07); terms are
de-duplicated per gene. dbCAN overview calls are '+'-joined with optional
'(start-end)' domain ranges (stripped); the --dbcan-consensus policy picks the
'Recommend Results' column ('recommended', calls supported by >= 2 tools) or
the union of the per-tool columns ('any'). dbCAN feeds the cazy ontology only
(eggNOG remains the EC source). A gene carrying several terms of one ontology
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
    ('COG_category', 'cog', 'eggnog_cog', 'cog'),
    ('EC', 'ec', 'eggnog_ec', 'comma'),
    ('PFAMs', 'pfam', 'eggnog_pfam', 'pfam'),
    ('CAZy', 'cazy', 'eggnog_cazy', 'comma'),
]

ONTOLOGIES = [ontology for _, ontology, _, _ in ONTOLOGY_FIELDS]

# run_dbcan versions whose overview layout this parser is verified against
# (design doc Section 4.3 interface record). Extend only after re-verifying
# the column list on real output of the new version.
KNOWN_DBCAN_VERSIONS = {'5.2.9'}

# The exact run_dbcan v5 overview.tsv header, verified on real 5.2.9 output
# (2026-08-30). Note: the tool's OVERVIEW_COLUMNS constant lists only the
# first 7 — the real file appends a 'Substrate' column (dbCAN-sub substrate
# mapping). Substrate is retained in the published file but not parsed here.
DBCAN_OVERVIEW_COLUMNS = ('Gene ID', 'EC#', 'dbCAN_hmm', 'dbCAN_sub',
                          'DIAMOND', '#ofTools', 'Recommend Results',
                          'Substrate')

# Trailing '(start-end)' domain-range suffix on dbCAN calls (e.g. GH5(100-300))
DBCAN_RANGE_RE = re.compile(r'\(\d+-\d+\)$')

# eggNOG-mapper v3 PFAMs values are '<pfam_name>_<start>_<end>' (domain
# coordinates appended; one value per domain hit, so a repeated domain appears
# once per copy). Verified 2026-10-07 on all 73.9 M pfam values of the eggNOG 7
# database (eggnog.db prots.pfam, which emapper passes through verbatim): every
# value matches, start <= end. The coordinates are always the LAST two fields,
# so names that themselves end in digits (AAA_12, Phage_holin_2_4) are safe.
# Design doc Q15: keeping the raw strings fragmented the pfam ontology by
# coordinates.
PFAM_DOMAIN_RE = re.compile(r'^(.+)_(\d+)_(\d+)$')

# eggNOG 7 COG_category values are either functional-category letters (in
# practice only 'S') or a COG ortholog-group id (e.g. COG0450, ~126 k OGs of
# the eggNOG 7 DB; 75-80 % of annotated genes on real data). Ids are mapped to
# their category letters with NCBI's COG definitions table (--cog-def,
# cog-24.def.tab: COG id, category letters in order of importance, ...),
# which covers all 4,911 COG ids eggNOG 7 uses (design doc Q16, owner option
# b, 2026-10-07). Splitting an id per character produced bogus terms before.
COG_ID_RE = re.compile(r'^COG\d+$')
COG_LETTERS_RE = re.compile(r'^[A-Z]+$')
COG_DEF = {}         # COG id -> category letters, filled from --cog-def
COG_DEF_NAME = None  # basename of the --cog-def table (provenance)

# Wide-matrix specs: (file label, backend filter, ontology filter). The two
# CAZy backends get separate matrices, never one merged matrix (a gene called
# by both tools would double-count; owner decision, design doc Section 4.3).
WIDE_MATRIX_SPECS = (
    [(ontology, 'eggnog-mapper', ontology) for ontology in ONTOLOGIES]
    + [('cazy_dbcan', 'run_dbcan', 'cazy')]
)

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


def parse_dbcan_versions_yml(path):
    """Extract run_dbcan and dbcan_db versions; enforce the known-layout set."""
    text = Path(path).read_text()
    tool_match = re.search(r'^\s*run_dbcan:\s*(\S+)\s*$', text, re.MULTILINE)
    db_match = re.search(r'^\s*dbcan_db:\s*(\S+)\s*$', text, re.MULTILINE)
    if not tool_match:
        fail(f"no 'run_dbcan:' version found in {path} — cannot verify the "
             f"overview layout (design doc 4.3 requires version-aware parsing)")
    tool_version = tool_match.group(1)
    if tool_version == 'unknown' or tool_version not in KNOWN_DBCAN_VERSIONS:
        fail(f"run_dbcan version '{tool_version}' is not among the layouts this "
             f"parser was verified against ({sorted(KNOWN_DBCAN_VERSIONS)}). "
             f"Re-verify the overview column layout on real output of that "
             f"version and extend KNOWN_DBCAN_VERSIONS before aggregating.")
    db_version = db_match.group(1) if db_match else 'unknown'
    return tool_version, db_version


def parse_dbcan_overview(path):
    """run_dbcan v5 <id>.overview.tsv -> one row per gene with raw call columns.

    Layout guard in the house style: the header must match the verified
    8-column v5 layout exactly, every data row must have 8 fields, and gene
    ids must be unique (the tool emits one row per gene). Header-only files
    (empty gene sets, RUN_DBCAN short-circuit) are valid.
    """
    header = None
    rows = []
    with open_text(path) as handle:
        for line_no, line in enumerate(handle, 1):
            line = line.rstrip('\n')
            if not line:
                continue
            fields = line.split('\t')
            if header is None:
                if tuple(fields) != DBCAN_OVERVIEW_COLUMNS:
                    fail(f"{path}: overview column header does not match the "
                         f"verified run_dbcan v5 layout "
                         f"({len(DBCAN_OVERVIEW_COLUMNS)} columns). Got "
                         f"{len(fields)} columns: {fields[:8]} — refusing to "
                         f"guess (design doc 4.3: fail loudly on layout drift)")
                header = fields
                continue
            if len(fields) != len(DBCAN_OVERVIEW_COLUMNS):
                fail(f"{path}: line {line_no} has {len(fields)} fields, expected "
                     f"{len(DBCAN_OVERVIEW_COLUMNS)}")
            record = dict(zip(DBCAN_OVERVIEW_COLUMNS, fields))
            rows.append({
                'gene_id': record['Gene ID'],
                'dbCAN_hmm': record['dbCAN_hmm'],
                'dbCAN_sub': record['dbCAN_sub'],
                'DIAMOND': record['DIAMOND'],
                'recommend': record['Recommend Results'],
            })
    if header is None:
        fail(f"{path}: no header line found — not a run_dbcan v5 overview file")
    overview = pd.DataFrame(rows, columns=['gene_id', 'dbCAN_hmm', 'dbCAN_sub',
                                           'DIAMOND', 'recommend'])
    duplicates = overview['gene_id'][overview['gene_id'].duplicated()]
    if not duplicates.empty:
        fail(f"{path}: duplicate gene ids in overview: "
             f"{sorted(duplicates.unique())[:5]}")
    return overview


def split_dbcan_terms(value):
    """One dbCAN overview call field -> CAZy accessions.

    Separators verified on real 5.2.9 output: hmm/sub/DIAMOND columns join
    domains with '+', 'Recommend Results' joins calls with '|', and ';' also
    occurs as a sub-separator — split on all three. hmm/sub calls carry
    '(start-end)' domain ranges, stripped here. '-' means no call.
    Subfamily ids (GH5_4) and dbCAN-sub cluster ids (GH78_e118) are kept as
    emitted.
    """
    if value in ('', '-'):
        return []
    terms = []
    for chunk in re.split(r'[+;|]', value):
        term = DBCAN_RANGE_RE.sub('', chunk.strip())
        if term and term != '-':
            terms.append(term)
    return terms


def explode_dbcan(overview, consensus):
    """Long form: one row per gene per dbCAN CAZy accession, per the consensus
    policy ('recommended' = the tool's >=2-tools column; 'any' = union of the
    per-tool columns). Same column shape as explode_annotations; evalue/score
    stay empty (the overview carries no single per-call value; Section 5.1
    allows empty).
    """
    rows = []
    for record in overview.itertuples(index=False):
        if consensus == 'recommended':
            accessions = split_dbcan_terms(record.recommend)
        else:
            accessions = []
            for value in (record.dbCAN_hmm, record.dbCAN_sub, record.DIAMOND):
                accessions.extend(split_dbcan_terms(value))
        seen = set()
        for accession in accessions:
            if accession in seen:
                continue
            seen.add(accession)
            rows.append({
                'gene_id': record.gene_id,
                'db': 'dbcan',
                'ontology': 'cazy',
                'accession': accession,
                'evalue': '',
                'score': '',
            })
    return pd.DataFrame(rows, columns=['gene_id', 'db', 'ontology', 'accession',
                                       'evalue', 'score'])


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
    """Split one annotation field into distinct terms ('-' and empty mean
    unannotated). Terms are de-duplicated per gene (order kept): a gene
    carrying a term contributes its abundance to it once (Section 4.8 step 5)
    - this matters for 'pfam', where a repeated domain yields one value per
    copy."""
    if value in ('', '-'):
        return []
    if mode == 'cog':
        if COG_ID_RE.match(value):
            if not COG_DEF:
                fail(f"COG_category value '{value}' is a COG id, but no --cog-def "
                     f"table was given to map it to category letters (design doc Q16)")
            if value not in COG_DEF:
                fail(f"COG id '{value}' is not in the --cog-def table "
                     f"({COG_DEF_NAME}) - use the COG release that covers the "
                     f"eggNOG database's ids (design doc Q16)")
            value = COG_DEF[value]
        elif not COG_LETTERS_RE.match(value):
            fail(f"COG_category value '{value}' is neither category letters nor "
                 f"a COG id - not the eggNOG-mapper v3 layout this parser was "
                 f"verified against (design doc Q16)")
        terms = list(value)
    else:
        terms = [term for term in value.split(',') if term and term != '-']
    if mode == 'pfam':
        names = []
        for term in terms:
            match = PFAM_DOMAIN_RE.match(term)
            if not match:
                fail(f"PFAMs value '{term}' is not '<pfam_name>_<start>_<end>', "
                     f"the eggNOG-mapper v3 layout this parser was verified "
                     f"against (design doc Q15) - re-verify the PFAMs format "
                     f"before aggregating")
            names.append(match.group(1))
        terms = names
    return list(dict.fromkeys(terms))


def parse_cog_def(path):
    """NCBI COG definitions table (cog-24.def.tab layout, no header): column 1
    the COG id, column 2 its functional-category letters. Latin-1 tolerant
    (older COG releases are not UTF-8); whitespace around the category is
    stripped (cog-24 has 'O ' for COG6144)."""
    mapping = {}
    with open(path, encoding='latin-1') as handle:
        for lineno, line in enumerate(handle, start=1):
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 2:
                fail(f"{path}:{lineno}: expected the NCBI COG definitions layout "
                     f"(COG id <TAB> category letters <TAB> ...)")
            cog_id, letters = fields[0].strip(), fields[1].strip()
            if not COG_ID_RE.match(cog_id) or not COG_LETTERS_RE.match(letters):
                fail(f"{path}:{lineno}: malformed COG definition "
                     f"({cog_id!r}, {letters!r})")
            mapping[cog_id] = letters
    if not mapping:
        fail(f"{path}: empty COG definitions table")
    return mapping


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
  # --ags is optional and may cover only a subset of samples; --dbcan is
  # optional as a whole set and requires --dbcan-versions-yml when given)
  aggregate_functions.py --assembly-mode assembly \\
      --counts s1.featureCounts.txt s2.featureCounts.txt \\
      --annotations s1.emapper.annotations s2.emapper.annotations \\
      --gffs s1.gff.gz s2.gff.gz \\
      --ags s1.ags.tsv s2.ags.tsv \\
      --dbcan s1.overview.tsv s2.overview.tsv \\
      --dbcan-versions-yml dbcan_versions.yml \\
      --cog-def cog-24.def.tab \\
      --eggnog-versions-yml versions.yml --output-dir .

  # Co-assembly mode (one shared annotations file, GFF and dbCAN overview)
  aggregate_functions.py --assembly-mode coassembly \\
      --counts s1.featureCounts.txt s2.featureCounts.txt \\
      --annotations coassembly.emapper.annotations \\
      --gffs coassembly.gff.gz \\
      --dbcan coassembly.overview.tsv \\
      --dbcan-versions-yml dbcan_versions.yml \\
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
    parser.add_argument('--dbcan', nargs='*', type=Path, default=[],
                        help='run_dbcan <id>.overview.tsv files (optional as a '
                             'set: absent entirely when --functional_cazy '
                             'false; when given, per sample under assembly '
                             'mode or exactly one shared file under '
                             'coassembly)')
    parser.add_argument('--dbcan-versions-yml', type=Path, default=None,
                        help='versions.yml from RUN_DBCAN (required when '
                             '--dbcan is given)')
    parser.add_argument('--dbcan-consensus', choices=['recommended', 'any'],
                        default='recommended',
                        help="Which dbCAN calls feed aggregation: "
                             "'recommended' uses the tool's Recommend Results "
                             "column (>= 2 tools), 'any' the union of the "
                             "per-tool columns")
    parser.add_argument('--cog-def', type=Path, default=None,
                        help='NCBI COG definitions table (cog-24.def.tab) used '
                             'to map COG ids in COG_category to category '
                             'letters (required whenever the annotations carry '
                             'COG ids, i.e. on any real eggNOG 7 output)')
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
                 *args.dbcan, args.eggnog_versions_yml,
                 *([args.cog_def] if args.cog_def else []),
                 *([args.dbcan_versions_yml] if args.dbcan_versions_yml else [])]:
        if not path.exists():
            fail(f"input file not found: {path}")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    global COG_DEF, COG_DEF_NAME
    if args.cog_def:
        COG_DEF = parse_cog_def(args.cog_def)
        COG_DEF_NAME = args.cog_def.name

    tool_version, db_version = parse_versions_yml(args.eggnog_versions_yml)

    if args.dbcan and not args.dbcan_versions_yml:
        fail("--dbcan given without --dbcan-versions-yml — the overview layout "
             "cannot be version-verified (design doc 4.3)")
    if args.dbcan_versions_yml:
        dbcan_tool_version, dbcan_db_version = \
            parse_dbcan_versions_yml(args.dbcan_versions_yml)
    else:
        dbcan_tool_version, dbcan_db_version = None, None

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
    dbcan_by_id = {}
    for path in args.dbcan:
        unit = sample_id_from(path, ['.overview.tsv'])
        if unit in dbcan_by_id:
            fail(f"duplicate --dbcan overview for '{unit}'")
        dbcan_by_id[unit] = parse_dbcan_overview(path)
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
    # 'assembly'; one shared co-assembly gene set under 'coassembly'. dbCAN
    # overviews follow the annotations shape (the tool runs on the same
    # proteins), but the whole set is optional (--functional_cazy false)
    if args.assembly_mode == 'assembly':
        for label, keys in [('annotations', annotations_by_id),
                            ('GFFs', gffs_by_id)]:
            if set(keys) != set(counts_by_sample):
                fail(f"sample ids of --counts {sorted(counts_by_sample)} and "
                     f"--{label.lower()} {sorted(keys)} do not match")
        if dbcan_by_id and set(dbcan_by_id) != set(counts_by_sample):
            fail(f"sample ids of --counts {sorted(counts_by_sample)} and "
                 f"--dbcan {sorted(dbcan_by_id)} do not match")
        gff_for = dict(gffs_by_id)
        annotations_for = dict(annotations_by_id)
    else:
        if len(annotations_by_id) != 1 or len(gffs_by_id) != 1:
            fail(f"coassembly mode expects exactly one annotations file and one "
                 f"GFF, got {len(annotations_by_id)} and {len(gffs_by_id)}")
        if len(dbcan_by_id) > 1:
            fail(f"coassembly mode expects at most one dbCAN overview, got "
                 f"{len(dbcan_by_id)}")
        shared_gff = next(iter(gffs_by_id.values()))
        shared_annotations = next(iter(annotations_by_id.values()))
        gff_for = {sample: shared_gff for sample in sample_ids}
        annotations_for = {sample: shared_annotations for sample in sample_ids}

    # dbCAN exploded frames, per annotations unit and fanned out per sample
    # (empty frame = no dbCAN input for that sample)
    empty_exploded = pd.DataFrame(columns=['gene_id', 'db', 'ontology',
                                           'accession', 'evalue', 'score'])
    dbcan_exploded_by_id = {
        unit: explode_dbcan(table, args.dbcan_consensus)
        for unit, table in dbcan_by_id.items()
    }
    if args.assembly_mode == 'assembly':
        dbcan_exploded_for = {
            sample: dbcan_exploded_by_id.get(sample, empty_exploded)
            for sample in sample_ids
        }
    else:
        shared_dbcan_exploded = (next(iter(dbcan_exploded_by_id.values()))
                                 if dbcan_exploded_by_id else empty_exploded)
        dbcan_exploded_for = {sample: shared_dbcan_exploded
                              for sample in sample_ids}

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
        stray_dbcan = set(dbcan_exploded_for[sample]['gene_id']) - gff_ids
        if stray_dbcan:
            fail(f"sample '{sample}': dbCAN gene ids not present in the GFF: "
                 f"{sorted(stray_dbcan)[:5]}")

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

    # --- Section 5.1: gene annotations (long, per annotations unit; eggNOG
    #     and dbCAN parts carry their own per-row tool provenance) ---
    exploded_by_id = {
        unit: explode_annotations(table)
        for unit, table in annotations_by_id.items()
    }
    annotation_parts = []
    for unit in sorted(set(exploded_by_id) | set(dbcan_exploded_by_id)):
        # Under coassembly the single annotations unit pairs with the single
        # shared GFF whatever either file is named; under assembly the ids match
        genes = next(iter(gffs_by_id.values())) \
            if args.assembly_mode == 'coassembly' else gff_for[unit]
        for exploded, tool, part_tool_version, part_db_version in [
                (exploded_by_id.get(unit), 'eggnog-mapper',
                 tool_version, db_version),
                (dbcan_exploded_by_id.get(unit), 'run_dbcan',
                 dbcan_tool_version, dbcan_db_version)]:
            if exploded is None or exploded.empty:
                continue
            part = exploded.merge(
                genes[['gene_id', 'contig_id', 'start', 'end', 'strand',
                       'partial', 'length_bp']],
                on='gene_id')
            part.insert(0, 'sample_id', unit)
            part['tool'] = tool
            part['tool_version'] = part_tool_version
            part['db_version'] = part_db_version
            annotation_parts.append(part)
    if annotation_parts:
        gene_annotations = pd.concat(annotation_parts, ignore_index=True)
    else:
        gene_annotations = pd.DataFrame(columns=['sample_id', 'gene_id', 'db',
                                                 'ontology', 'accession', 'evalue',
                                                 'score', 'contig_id', 'start',
                                                 'end', 'strand', 'partial',
                                                 'length_bp', 'tool',
                                                 'tool_version', 'db_version'])
    gene_annotations['description'] = ''  # v3 dropped Description (design doc 4.8 record)
    if COG_DEF_NAME:
        # cog categories are eggNOG calls mapped through the COG table (Q16)
        cog_rows = gene_annotations['db'] == 'eggnog_cog'
        gene_annotations.loc[cog_rows, 'db_version'] = \
            gene_annotations.loc[cog_rows, 'db_version'].astype(str) + f'; {COG_DEF_NAME}'
    gene_annotations = gene_annotations.sort_values(
        ['sample_id', 'gene_id', 'db', 'accession'], kind='mergesort',
        ignore_index=True)
    gene_annotations[['sample_id', 'gene_id', 'contig_id', 'start', 'end',
                      'strand', 'length_bp', 'partial', 'db', 'accession',
                      'description', 'evalue', 'score', 'tool', 'tool_version',
                      'db_version']].to_csv(
        args.output_dir / 'gene_annotations.tsv', sep='\t', index=False)

    # --- Section 5.3: function abundance (per-sample term explode x TPM/CPGE;
    #     the backend column separates eggnog-mapper from run_dbcan rows) ---
    term_parts = []
    for sample in sample_ids:
        sample_abundance = abundance.loc[abundance['sample_id'] == sample,
                                         ['gene_id', 'tpm', 'cpge_num']]
        for exploded, backend in [
                (explode_annotations(annotations_for[sample]), 'eggnog-mapper'),
                (dbcan_exploded_for[sample], 'run_dbcan')]:
            if exploded.empty:
                continue
            part = exploded[['gene_id', 'ontology', 'accession']].merge(
                sample_abundance, on='gene_id')
            part.insert(0, 'sample_id', sample)
            part['backend'] = backend
            term_parts.append(part)
    if term_parts:
        terms = pd.concat(term_parts, ignore_index=True)
        # Intentional double-counting (Section 4.8 step 5): a gene carrying
        # several terms of one ontology contributes its full TPM/CPGE to each
        function_abundance = (
            terms.groupby(['sample_id', 'backend', 'ontology', 'accession'],
                          as_index=False)[['tpm', 'cpge_num']].sum()
            .rename(columns={'tpm': 'abundance_tpm'}))
    else:
        function_abundance = pd.DataFrame(columns=['sample_id', 'backend',
                                                   'ontology', 'accession',
                                                   'abundance_tpm', 'cpge_num'])
    function_abundance['source'] = 'contigs'
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
        ['sample_id', 'ontology', 'backend', 'accession'], kind='mergesort',
        ignore_index=True)
    out_53 = function_abundance[['sample_id', 'source', 'backend', 'ontology',
                                 'accession', 'description', 'abundance_tpm',
                                 'abundance_cpge', 'abundance_native',
                                 'native_unit']].copy()
    out_53['abundance_tpm'] = out_53['abundance_tpm'].map(TPM_FORMAT.format)
    out_53.to_csv(args.output_dir / 'function_abundance.tsv', sep='\t',
                  index=False)

    # --- Wide matrices per ontology/backend (every sample a column, missing
    #     -> 0; CPGE cells are blank for samples without AGS). The cazy_dbcan
    #     matrices carry the run_dbcan calls; the plain cazy ones stay
    #     eggNOG-derived — never merged (design doc 4.3). With --dbcan absent
    #     the cazy_dbcan matrices are emitted header-only, keeping the output
    #     file set deterministic ---
    ags_samples = [sample for sample in sample_ids if sample in ge_by_sample]
    for label, backend, ontology in WIDE_MATRIX_SPECS:
        subset = function_abundance[
            (function_abundance['ontology'] == ontology)
            & (function_abundance['backend'] == backend)]
        wide = subset.pivot_table(index='accession', columns='sample_id',
                                  values='abundance_tpm', aggfunc='sum',
                                  fill_value=0.0)
        wide = wide.reindex(columns=sample_ids, fill_value=0.0).sort_index()
        wide.index.name = 'accession'
        # float_format instead of per-cell mapping: DataFrame.map needs
        # pandas >= 2.1 and the pinned image ships 2.0
        wide.to_csv(args.output_dir / f'function_wide_{label}_tpm.tsv',
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
        wide_cpge.to_csv(args.output_dir / f'function_wide_{label}_cpge.tsv',
                         sep='\t', header=True)

    # --- Annotated fraction per sample, by count and by abundance ---
    fraction_rows = []
    for sample in sample_ids:
        sample_genes = abundance[abundance['sample_id'] == sample]
        genes_total = len(sample_genes)
        tpm_total = sample_genes['tpm'].sum()
        exploded_dbcan = dbcan_exploded_for[sample]
        exploded = pd.concat([explode_annotations(annotations_for[sample]),
                              exploded_dbcan], ignore_index=True)
        # 'cazy' counts a CAZy call from either backend; 'cazy_dbcan' counts
        # run_dbcan alone (always present, 0 when --dbcan is absent)
        for ontology in ['any'] + ONTOLOGIES + ['cazy_dbcan']:
            if exploded.empty:
                annotated_ids = set()
            elif ontology == 'any':
                annotated_ids = set(exploded['gene_id'])
            elif ontology == 'cazy_dbcan':
                annotated_ids = set(exploded_dbcan['gene_id'])
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
    dbcan_rows = int((function_abundance['backend'] == 'run_dbcan').sum())
    dbcan_note = (f"{dbcan_rows} dbCAN CAZy rows" if args.dbcan
                  else "dbCAN input absent")
    print(f"Aggregated {len(sample_ids)} sample(s), "
          f"{len(gene_abundance)} gene abundance rows, "
          f"{len(out_53)} function abundance rows "
          f"({annotated_any}/{len(sample_ids)} samples with annotations; "
          f"{dbcan_note}; "
          f"{len(ags_samples)}/{len(sample_ids)} samples with CPGE)")


if __name__ == '__main__':
    main()
