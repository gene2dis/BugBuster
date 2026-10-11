#!/usr/bin/env python3
"""Aggregate read-level functional profiles into the canonical function
abundance schema (design doc Sections 4.6, 5.3, 6.2; tasks T8a — Woltka,
T8b — SUPER-FOCUS, T8c — HUMAnN). One backend per run; BACKEND_SPECS below holds every
backend-specific detail (file suffixes, versions.yml keys, verified versions,
ontologies, whether CPGE is computable).

The read branch is reported SEPARATELY from the contig branch (Section 2:
"Report separately, do not merge"): this script writes its own read_*
tables and never touches function_abundance.tsv.

Inputs per sample (<stem> = woltka | superfocus | humann):
  - <id>.<stem>_functions.tsv   sample_id ontology accession description count rpk
                                (woltka: WOLTKA_CLASSIFY, bin/woltka_function_profile.py;
                                superfocus: SUPERFOCUS, bin/superfocus_function_profile.py,
                                rpk always blank; humann: HUMANN,
                                bin/humann_function_profile.py, count = rpk = HUMAnN RPK)
  - <id>.<stem>_summary.tsv     sample_id ontology reads_assigned reads_annotated
                                fraction_annotated
  - <id>.woltka_unassigned.tsv  sample_id reads_unassigned (woltka only)
  - <id>.ags.tsv                MicrobeCensus AGS table (OPTIONAL per sample:
                                missing samples get blank CPGE, never 0)
plus the backend process's versions.yml (version-aware guard, below).

SUPER-FOCUS (Q17): SEED levels as ontologies seed_level1..3; there is no gene
length behind a SEED hit, so no RPK and no CPGE: abundance_cpge stays blank,
no _cpge wide matrices are written, cpge_status is 'not_applicable' and
reads_unassigned_ambiguous is blank. AGS / genome equivalents are still
reported in read_sample_summary.tsv when MicrobeCensus produced them.

HUMAnN (T8c): ontologies ko / ec / metacyc in RPK (native_unit = rpk, HUMAnN
run with --count-normalization RPKs); CPGE = RPK / GE as for Woltka. Its
summary carries read counts on the 'any' row only; the per-ontology rows have
blank read columns and an RPK-share fraction_annotated computed by the
composer (no per-family read counts exist under RPKs), which is passed
through as-is (fraction_mode 'composer').

Outputs (TSV with header, deterministic ordering):
  - read_function_abundance.tsv      Section 5.3 schema, source=reads,
                                     backend=woltka; abundance_native = Woltka
                                     read count (mates counted separately;
                                     fractional under the default 1/k
                                     multi-hit division), native_unit=reads;
                                     abundance_cpge = RPK / genome equivalents
                                     (Section 6.2, blank without AGS);
                                     abundance_tpm blank (contig branch only)
  - read_function_wide_<ont>_native.tsv  wide read-count matrix per ontology
  - read_function_wide_<ont>_cpge.tsv    wide CPGE matrix (blank columns for
                                         samples without AGS; CPGE backends only)
  - read_annotated_fraction.tsv      per sample: 'any' + per-ontology share of
                                     ORF-assigned reads carrying a term
  - read_sample_summary.tsv          per sample: reads assigned to ORFs, reads
                                     left unassigned by --uniq ambiguity, AGS /
                                     GE and the CPGE status (ok/unavailable —
                                     the recorded MicrobeCensus-fallback warning)

Version-aware parsing (design doc 4.6.1 discipline, applied to every
backend): the backend tool version is read from its process's versions.yml
and must be one whose output this pipeline was verified against; every input
header must match its expected layout exactly. Any mismatch is a hard error.
"""

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

# Tool versions whose output (and the composer's intermediate tables built
# from it) were verified (design doc Section 4.6.2 / 4.6.3 records). Extend
# only after re-verifying on real output of the new version.
KNOWN_WOLTKA_VERSIONS = {'0.1.7'}
KNOWN_SUPERFOCUS_VERSIONS = {'1.8'}
KNOWN_HUMANN_VERSIONS = {'4.0.0.alpha.2'}

BACKEND_SPECS = {
    'woltka': {
        'stem': 'woltka',
        'tool_key': 'woltka',
        'db_key': 'wol_db',
        'known_versions': KNOWN_WOLTKA_VERSIONS,
        'ontologies': ['ko', 'ec', 'cog', 'pfam', 'metacyc'],
        'cpge': True,          # RPK over the WoLr2 ORF length
        'unassigned': True,    # reads left unassigned by --uniq ambiguity
        'native_unit': 'reads',
        'fraction_mode': 'reads',
    },
    'superfocus': {
        'stem': 'superfocus',
        'tool_key': 'superfocus',
        'db_key': 'superfocus_db',
        'known_versions': KNOWN_SUPERFOCUS_VERSIONS,
        'ontologies': ['seed_level1', 'seed_level2', 'seed_level3'],
        'cpge': False,         # no gene length behind a SEED hit (Q17)
        'unassigned': False,
        'native_unit': 'reads',
        'fraction_mode': 'reads',
    },
    'humann': {
        'stem': 'humann',
        'tool_key': 'humann',
        'db_key': 'humann_db',
        'known_versions': KNOWN_HUMANN_VERSIONS,
        'ontologies': ['ko', 'ec', 'metacyc'],
        'cpge': True,          # HUMAnN RPK (--count-normalization RPKs) / GE
        'unassigned': False,
        'native_unit': 'rpk',
        # 'any' row in reads; per-ontology rows: RPK share from the composer
        'fraction_mode': 'composer',
    },
}

FUNCTIONS_COLUMNS = ['sample_id', 'ontology', 'accession', 'description', 'count', 'rpk']
SUMMARY_COLUMNS = ['sample_id', 'ontology', 'reads_assigned', 'reads_annotated',
                   'fraction_annotated']
UNASSIGNED_COLUMNS = ['sample_id', 'reads_unassigned']
AGS_COLUMNS = ['sample_id', 'average_genome_size_bp', 'genome_equivalents',
               'total_bases']
OUT_53_COLUMNS = ['sample_id', 'source', 'backend', 'ontology', 'accession',
                  'description', 'abundance_tpm', 'abundance_cpge',
                  'abundance_native', 'native_unit']

NATIVE_FORMAT = '{:.6f}'
CPGE_FORMAT = '{:.6f}'
FRACTION_FORMAT = '{:.4f}'


def fail(message):
    print(f"Error: {message}", file=sys.stderr)
    sys.exit(1)


def sample_id_from(path, suffix):
    name = Path(path).name
    if not name.endswith(suffix):
        fail(f"cannot derive sample id from '{name}' (expected suffix '{suffix}')")
    return name[: -len(suffix)]


def parse_versions_yml(path, backend, spec):
    text = Path(path).read_text()
    tool_key, db_key = spec['tool_key'], spec['db_key']
    tool_match = re.search(rf'^\s*{tool_key}:\s*(\S+)\s*$', text, re.MULTILINE)
    db_match = re.search(rf'^\s*{db_key}:\s*(.+?)\s*$', text, re.MULTILINE)
    if not tool_match:
        fail(f"no '{tool_key}:' version found in {path} — cannot verify the output "
             f"layout (design doc 4.6 requires version-aware parsing)")
    version = tool_match.group(1)
    if version not in spec['known_versions']:
        constant = f"KNOWN_{backend.upper()}_VERSIONS"
        fail(f"{tool_key} version '{version}' is not among the versions this parser "
             f"was verified against ({sorted(spec['known_versions'])}). Re-verify "
             f"the {backend} output on real data of that version and extend "
             f"{constant} before aggregating.")
    return version, (db_match.group(1) if db_match else 'unknown')


def read_table(path, columns, suffix):
    """Read a per-sample TSV, enforcing the exact header and that every row
    carries the filename-derived sample id."""
    sample = sample_id_from(path, suffix)
    table = pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)
    if list(table.columns) != columns:
        fail(f"{path}: unexpected header {list(table.columns)} (expected {columns})")
    wrong = table.loc[table['sample_id'] != sample, 'sample_id']
    if not wrong.empty:
        fail(f"{path}: row sample_id '{wrong.iloc[0]}' does not match the "
             f"filename-derived id '{sample}'")
    return sample, table


def to_float(table, column, path, allow_empty=False):
    def convert(raw):
        if raw == '' and allow_empty:
            return float('nan')
        try:
            value = float(raw)
        except ValueError:
            fail(f"{path}: non-numeric {column} value '{raw}'")
        if value < 0:
            fail(f"{path}: negative {column} value '{raw}'")
        return value
    return table[column].map(convert)


def parse_ags(path):
    sample = sample_id_from(path, '.ags.tsv')
    with open(path) as handle:
        lines = [line.rstrip('\n') for line in handle if line.strip()]
    if not lines or lines[0].split('\t') != AGS_COLUMNS:
        fail(f"{path}: unexpected ags.tsv header (expected {AGS_COLUMNS})")
    if len(lines) != 2:
        fail(f"{path}: expected exactly one data row, found {len(lines) - 1}")
    fields = lines[1].split('\t')
    if len(fields) != len(AGS_COLUMNS) or fields[0] != sample:
        fail(f"{path}: malformed data row or sample_id not matching '{sample}'")
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


def collect(paths, columns, suffix, label):
    tables = {}
    for path in paths:
        sample, table = read_table(path, columns, suffix)
        if sample in tables:
            fail(f"duplicate {label} input for sample '{sample}'")
        tables[sample] = (path, table)
    return tables


def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Aggregate read-level functional profiles (Woltka, '
                    'SUPER-FOCUS or HUMAnN backend) into the canonical function abundance schema.')
    parser.add_argument('--backend', required=True, choices=sorted(BACKEND_SPECS))
    parser.add_argument('--functions', nargs='+', required=True,
                        help='<id>.<backend>_functions.tsv files')
    parser.add_argument('--summaries', nargs='+', required=True,
                        help='<id>.<backend>_summary.tsv files')
    parser.add_argument('--unassigned', nargs='*', default=[],
                        help='<id>.woltka_unassigned.tsv files (woltka only, required there)')
    parser.add_argument('--ags', nargs='*', default=[],
                        help='MicrobeCensus <id>.ags.tsv files (optional per sample)')
    parser.add_argument('--versions-yml', required=True,
                        help="the backend process's versions.yml (version-aware guard)")
    parser.add_argument('--output-dir', type=Path, default=Path('.'))
    return parser.parse_args()


def main():
    args = parse_arguments()
    spec = BACKEND_SPECS[args.backend]
    stem = spec['stem']
    ontologies = spec['ontologies']
    tool_version, db_version = parse_versions_yml(args.versions_yml, args.backend, spec)

    if spec['unassigned'] and not args.unassigned:
        fail(f"--unassigned <id>.{stem}_unassigned.tsv files are required for the "
             f"{args.backend} backend")
    if not spec['unassigned'] and args.unassigned:
        fail(f"--unassigned is not used by the {args.backend} backend")

    functions = collect(args.functions, FUNCTIONS_COLUMNS, f'.{stem}_functions.tsv', 'functions')
    summaries = collect(args.summaries, SUMMARY_COLUMNS, f'.{stem}_summary.tsv', 'summary')
    unassigned = collect(args.unassigned, UNASSIGNED_COLUMNS, f'.{stem}_unassigned.tsv',
                         'unassigned')

    sample_ids = sorted(functions)
    checks = [('summary', summaries)] + ([('unassigned', unassigned)] if spec['unassigned'] else [])
    for label, tables in checks:
        if sorted(tables) != sample_ids:
            fail(f"the {label} inputs cover samples {sorted(tables)} but the "
                 f"functions inputs cover {sample_ids}; every sample needs "
                 f"{'all three' if spec['unassigned'] else 'both'}")

    ags_by_sample = {}
    for path in args.ags:
        sample, values = parse_ags(path)
        if sample in ags_by_sample:
            fail(f"duplicate ags.tsv input for sample '{sample}'")
        if sample not in sample_ids:
            fail(f"{path}: ags.tsv for unknown sample '{sample}' (read-branch "
                 f"samples: {sample_ids})")
        ags_by_sample[sample] = values
    ge_by_sample = {s: v['genome_equivalents_float'] for s, v in ags_by_sample.items()}
    ags_samples = [s for s in sample_ids if s in ge_by_sample]
    missing_ags = [s for s in sample_ids if s not in ge_by_sample]
    if missing_ags and spec['cpge']:
        print(f"Warning: no MicrobeCensus AGS for sample(s) {', '.join(missing_ags)} "
              f"— their abundance_cpge is left blank (non-fatal fallback, design "
              f"doc Section 4.7)", file=sys.stderr)

    # --- Section 5.3 long table ---
    frames = []
    for sample in sample_ids:
        path, table = functions[sample]
        bad = sorted(set(table['ontology']) - set(ontologies))
        if bad:
            fail(f"{path}: unknown ontology value(s) {bad} (expected {ontologies})")
        dup = table.duplicated(['ontology', 'accession'])
        if dup.any():
            row = table[dup].iloc[0]
            fail(f"{path}: duplicate row for {row['ontology']} {row['accession']}")
        table = table.copy()
        table['count_num'] = to_float(table, 'count', path)
        if spec['cpge']:
            table['rpk_num'] = to_float(table, 'rpk', path)
        else:
            filled = table.loc[table['rpk'] != '', 'rpk']
            if not filled.empty:
                fail(f"{path}: rpk must be blank for the {args.backend} backend "
                     f"(no gene length, CPGE not computable), got '{filled.iloc[0]}'")
            table['rpk_num'] = float('nan')
        frames.append(table)
    long = (pd.concat(frames, ignore_index=True) if frames
            else pd.DataFrame(columns=FUNCTIONS_COLUMNS + ['count_num', 'rpk_num']))

    def cpge_cell(row):
        if not spec['cpge']:
            return ''
        ge = ge_by_sample.get(row['sample_id'])
        return '' if ge is None else CPGE_FORMAT.format(row['rpk_num'] / ge)

    out = pd.DataFrame({
        'sample_id': long['sample_id'],
        'source': 'reads',
        'backend': args.backend,
        'ontology': long['ontology'],
        'accession': long['accession'],
        'description': long['description'],
        'abundance_tpm': '',
        'abundance_cpge': long.apply(cpge_cell, axis=1) if not long.empty else [],
        'abundance_native': [NATIVE_FORMAT.format(v) for v in long['count_num']],
        'native_unit': spec['native_unit'],
    }, columns=OUT_53_COLUMNS)
    out = out.sort_values(['sample_id', 'ontology', 'backend', 'accession'],
                          kind='mergesort')
    out.to_csv(args.output_dir / 'read_function_abundance.tsv', sep='\t', index=False)

    # --- Wide matrices per ontology (every sample a column; missing -> 0 for
    #     read counts, blank CPGE columns for AGS-less samples). Every file of
    #     the backend's set is always written, header-only when an ontology has
    #     no rows (woltka: 10; superfocus: 3, no CPGE) ---
    for ontology in ontologies:
        subset = long[long['ontology'] == ontology]
        native = subset.pivot_table(index='accession', columns='sample_id',
                                    values='count_num', aggfunc='sum', fill_value=0.0)
        # astype: pivot_table downcasts all-integral values to int64, which
        # would bypass float_format
        native = native.reindex(columns=sample_ids, fill_value=0.0).sort_index().astype(float)
        native.index.name = 'accession'
        native.to_csv(args.output_dir / f'read_function_wide_{ontology}_native.tsv',
                      sep='\t', header=True, float_format='%.6f')

        if not spec['cpge']:
            continue
        cpge = subset[subset['sample_id'].isin(ags_samples)].copy()
        cpge['cpge_num'] = [r / ge_by_sample[s] for r, s in
                            zip(cpge['rpk_num'], cpge['sample_id'])]
        wide = cpge.pivot_table(index='accession', columns='sample_id',
                                values='cpge_num', aggfunc='sum', fill_value=0.0)
        wide = wide.reindex(index=native.index, columns=sample_ids)
        if ags_samples:
            wide[ags_samples] = wide[ags_samples].fillna(0.0)
        wide.index.name = 'accession'
        # cell-wise formatting: pandas 2.0 (pinned image) drops float_format
        # when na_rep is also given (aggregate_functions.py precedent)
        for column in wide.columns:
            wide[column] = ['' if pd.isna(v) else CPGE_FORMAT.format(v) for v in wide[column]]
        wide.to_csv(args.output_dir / f'read_function_wide_{ontology}_cpge.tsv',
                    sep='\t', header=True)

    # --- Annotated fraction (share of ORF-assigned reads carrying a term) ---
    fraction_rows = []
    for sample in sample_ids:
        path, table = summaries[sample]
        expected = ['any'] + ontologies
        if list(table['ontology']) != expected:
            fail(f"{path}: expected ontology rows {expected}, got {list(table['ontology'])}")
        if spec['fraction_mode'] == 'composer':
            # 'any' row: reads; ontology rows: blank reads, composer RPK share
            any_row = table.iloc[[0]]
            assigned = to_float(any_row, 'reads_assigned', path).iloc[0]
            annotated = to_float(any_row, 'reads_annotated', path).iloc[0]
            if annotated > assigned + 1e-6:
                fail(f"{path}: reads_annotated exceeds reads_assigned for any")
            fraction_rows.append({
                'sample_id': sample, 'backend': args.backend, 'ontology': 'any',
                'reads_assigned': NATIVE_FORMAT.format(assigned),
                'reads_annotated': NATIVE_FORMAT.format(annotated),
                'fraction_annotated': (FRACTION_FORMAT.format(annotated / assigned)
                                       if assigned > 0 else ''),
            })
            for _, row in table.iloc[1:].iterrows():
                if row['reads_assigned'] != '' or row['reads_annotated'] != '':
                    fail(f"{path}: {row['ontology']} row must have blank read columns "
                         f"for the {args.backend} backend (RPK-share fraction only)")
                raw = row['fraction_annotated']
                if raw != '':
                    try:
                        value = float(raw)
                    except ValueError:
                        fail(f"{path}: non-numeric fraction_annotated '{raw}'")
                    if not 0.0 <= value <= 1.0 + 1e-9:
                        fail(f"{path}: fraction_annotated '{raw}' outside [0, 1]")
                fraction_rows.append({
                    'sample_id': sample, 'backend': args.backend,
                    'ontology': row['ontology'], 'reads_assigned': '',
                    'reads_annotated': '',
                    'fraction_annotated': '' if raw == '' else FRACTION_FORMAT.format(float(raw)),
                })
            continue
        assigned = to_float(table, 'reads_assigned', path)
        annotated = to_float(table, 'reads_annotated', path)
        for ontology, a, n in zip(table['ontology'], assigned, annotated):
            if n > a + 1e-6:
                fail(f"{path}: reads_annotated exceeds reads_assigned for {ontology}")
            fraction_rows.append({
                'sample_id': sample,
                'backend': args.backend,
                'ontology': ontology,
                'reads_assigned': NATIVE_FORMAT.format(a),
                'reads_annotated': NATIVE_FORMAT.format(n),
                'fraction_annotated': FRACTION_FORMAT.format(n / a) if a > 0 else '',
            })
    pd.DataFrame(fraction_rows, columns=['sample_id', 'backend', 'ontology',
                                         'reads_assigned', 'reads_annotated',
                                         'fraction_annotated']).to_csv(
        args.output_dir / 'read_annotated_fraction.tsv', sep='\t', index=False)

    # --- Per-sample summary incl. the CPGE status (recorded fallback warning) ---
    summary_rows = []
    for sample in sample_ids:
        reads_unassigned = None
        if spec['unassigned']:
            upath, utable = unassigned[sample]
            if len(utable) != 1:
                fail(f"{upath}: expected exactly one data row, found {len(utable)}")
            reads_unassigned = to_float(utable, 'reads_unassigned', upath).iloc[0]
        spath, stable = summaries[sample]
        # the 'any' row only: for the humann backend the per-ontology rows
        # have blank read columns
        reads_assigned = to_float(stable.iloc[[0]], 'reads_assigned', spath).iloc[0]
        values = ags_by_sample.get(sample)
        summary_rows.append({
            'sample_id': sample,
            'backend': args.backend,
            'tool_version': tool_version,
            'db_version': db_version,
            'reads_assigned': NATIVE_FORMAT.format(reads_assigned),
            'reads_unassigned_ambiguous': ('' if reads_unassigned is None
                                           else NATIVE_FORMAT.format(reads_unassigned)),
            'average_genome_size_bp': values['average_genome_size_bp'] if values else '',
            'genome_equivalents': values['genome_equivalents'] if values else '',
            'cpge_status': ('not_applicable' if not spec['cpge']
                            else 'ok' if values else 'unavailable'),
        })
    pd.DataFrame(summary_rows, columns=[
        'sample_id', 'backend', 'tool_version', 'db_version', 'reads_assigned',
        'reads_unassigned_ambiguous', 'average_genome_size_bp', 'genome_equivalents',
        'cpge_status']).to_csv(args.output_dir / 'read_sample_summary.tsv',
                               sep='\t', index=False)

    cpge_note = (f"{len(ags_samples)}/{len(sample_ids)} samples with CPGE" if spec['cpge']
                 else "CPGE not applicable for this backend")
    print(f"Aggregated {len(sample_ids)} sample(s) ({args.backend} {tool_version}, "
          f"{db_version}): {len(out)} function abundance rows; {cpge_note}")


if __name__ == '__main__':
    main()
