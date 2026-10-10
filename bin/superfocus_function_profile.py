#!/usr/bin/env python3
"""Turn one sample's SUPER-FOCUS output into per-level SEED read abundances.

Read-level functional branch, SUPER-FOCUS backend (design doc Section 4.6.3,
task T8b, Q17). Runs inside SUPERFOCUS, in the pinned SUPER-FOCUS container
(stdlib only), right after `superfocus` profiled the sample's reads.

Only the raw count column of <prefix>all_levels_and_function.xls is used
(SUPER-FOCUS 1.8, verified in the source and on real output):
  - its '%' columns are computed with numpy.divide(where=...) without an
    out= array, so they are undefined for a zero-total sample;
  - its subsystem_level_N files aggregate by level NAME only, so the level-2
    placeholder '-' (present under 33 level-1 categories) is merged into one
    bogus row.
The levels are therefore summed here from the function-level rows, keyed by
their full path:
  seed_level1  accession 'L1'
  seed_level2  accession 'L1 | L2'
  seed_level3  accession 'L1 | L2 | L3'
with the level's own name as description (owner decision, Q17).

Counting semantics are SUPER-FOCUS's own (-n 1, the default): every read with
an accepted hit contributes exactly 1, divided 1/k across the k distinct
(level 1, level 2, level 3, function) keys of its equal-best hits. A read
therefore adds at most 1 to any level, and the per-level totals equal the
number of reads with an accepted hit. Mates are separate reads.

There is no gene length behind a SEED hit, so no RPK exists and CPGE
(Section 6.2) is not computable for this backend (Q17): the rpk column is
left blank.

Outputs:
  <prefix>.superfocus_functions.tsv  sample_id ontology accession description
                                     count rpk (rpk always blank)
  <prefix>.superfocus_summary.tsv    sample_id ontology reads_assigned
                                     reads_annotated fraction_annotated
                                     ('any' row + one row per level;
                                     reads_assigned = reads given to
                                     SUPER-FOCUS, mates counted separately;
                                     reads_annotated = reads with an
                                     accepted hit)

Usage:
  superfocus_function_profile.py --table output_all_levels_and_function.xls
      --query S1.fastq --input-reads 1000 --aligner diamond --database 90
      --sample-id S1 --prefix S1
  (omit --table for a sample with no reads: header/zero tables only)
"""

import argparse
import sys
from collections import defaultdict

ONTOLOGIES = ['seed_level1', 'seed_level2', 'seed_level3']
LEVEL_COLUMNS = ['Subsystem Level 1', 'Subsystem Level 2', 'Subsystem Level 3', 'Function']
FUNCTIONS_COLUMNS = ['sample_id', 'ontology', 'accession', 'description', 'count', 'rpk']
SUMMARY_COLUMNS = ['sample_id', 'ontology', 'reads_assigned', 'reads_annotated',
                   'fraction_annotated']
PATH_SEPARATOR = ' | '

# --superfocus_aligner value -> the aligner name SUPER-FOCUS writes in its
# 'Aligner used:' line (it calls mmseqs2 by its command name)
ALIGNER_NAMES = {'diamond': 'diamond', 'mmseqs2': 'mmseqs'}


def die(msg):
    print(f"Error: {msg}", file=sys.stderr)
    sys.exit(1)


def fmt(value):
    # full precision: the aggregator formats the published tables
    return f"{value:.12g}"


def read_table(path, query, aligner, database):
    """Parse the all-levels table, enforcing SUPER-FOCUS 1.8's exact layout:
    'Query: ...', 'Database used: <db>', 'Aligner used: <aligner>', a blank
    separator row (written by csv.writer as '""'), then the header (levels, function, one count column named after the
    query file, its '%' column)."""
    with open(path, encoding='utf-8') as handle:
        lines = handle.read().split('\n')
    if lines and lines[-1] == '':
        lines.pop()
    if len(lines) < 5:
        die(f"{path}: expected 4 metadata lines and a header, found {len(lines)} line(s)")
    if not lines[0].startswith('Query: '):
        die(f"{path}: line 1 should start with 'Query: ', got '{lines[0]}'")
    if lines[1] != f"Database used: {database}":
        die(f"{path}: line 2 is '{lines[1]}', expected 'Database used: {database}'")
    if lines[2] != f"Aligner used: {ALIGNER_NAMES[aligner]}":
        die(f"{path}: line 3 is '{lines[2]}', expected "
            f"'Aligner used: {ALIGNER_NAMES[aligner]}'")
    # csv.writer renders the blank separator row as a quoted empty field
    if lines[3] != '""':
        die(f"{path}: line 4 should be the blank separator row '\"\"', got '{lines[3]}'")
    expected = LEVEL_COLUMNS + [query, f"{query} %"]
    header = lines[4].split('\t')
    if header != expected:
        die(f"{path}: unexpected header {header} (expected {expected}) - the SUPER-FOCUS "
            f"output layout changed or the query was not the single file '{query}'")

    rows = []
    for number, line in enumerate(lines[5:], start=6):
        fields = line.split('\t')
        if len(fields) != len(expected):
            die(f"{path}: line {number} has {len(fields)} fields, expected {len(expected)}")
        levels = fields[:3]
        if any(level == '' for level in levels):
            die(f"{path}: line {number} has an empty subsystem level")
        try:
            count = float(fields[4])
        except ValueError:
            die(f"{path}: line {number} has a non-numeric count '{fields[4]}'")
        if count < 0:
            die(f"{path}: line {number} has a negative count '{fields[4]}'")
        rows.append((levels, count))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--table', help='<prefix>all_levels_and_function.xls')
    parser.add_argument('--query', required=True,
                        help='the single query file name given to superfocus -q')
    parser.add_argument('--input-reads', required=True, type=int,
                        help='number of reads in the query file')
    parser.add_argument('--aligner', required=True, choices=sorted(ALIGNER_NAMES))
    parser.add_argument('--database', required=True, help="cluster level, e.g. '90'")
    parser.add_argument('--sample-id', required=True)
    parser.add_argument('--prefix', required=True)
    args = parser.parse_args()

    if args.input_reads < 0:
        die(f"--input-reads must be >= 0, got {args.input_reads}")
    if args.table is None and args.input_reads != 0:
        die("--table is required unless the sample has no reads")

    rows = read_table(args.table, args.query, args.aligner, args.database) if args.table else []

    abundance = {ontology: defaultdict(float) for ontology in ONTOLOGIES}
    description = {}
    total = 0.0
    for levels, count in rows:
        total += count
        for depth, ontology in enumerate(ONTOLOGIES, start=1):
            accession = PATH_SEPARATOR.join(levels[:depth])
            abundance[ontology][accession] += count
            description[(ontology, accession)] = levels[depth - 1]

    # every read contributes at most 1 (Q17 counting semantics)
    if total > args.input_reads + 1e-6:
        die(f"{args.table}: total count {total} exceeds the {args.input_reads} input reads - "
            f"was the table produced with -n 0, or for a different query?")

    with open(f"{args.prefix}.superfocus_functions.tsv", 'w') as out:
        out.write('\t'.join(FUNCTIONS_COLUMNS) + '\n')
        for ontology in ONTOLOGIES:
            for accession in sorted(abundance[ontology]):
                value = abundance[ontology][accession]
                if value <= 0:
                    continue
                out.write('\t'.join([args.sample_id, ontology, accession,
                                     description[(ontology, accession)], fmt(value), '']) + '\n')

    fraction = f"{total / args.input_reads:.6f}" if args.input_reads > 0 else ''
    with open(f"{args.prefix}.superfocus_summary.tsv", 'w') as out:
        out.write('\t'.join(SUMMARY_COLUMNS) + '\n')
        for ontology in ['any'] + ONTOLOGIES:
            out.write('\t'.join([args.sample_id, ontology, str(args.input_reads), fmt(total),
                                 fraction]) + '\n')

    print(f"{args.sample_id}: {fmt(total)} of {args.input_reads} reads with an accepted SEED hit; "
          + ', '.join(f"{o} {len(abundance[o])} terms" for o in ONTOLOGIES))


if __name__ == '__main__':
    main()
