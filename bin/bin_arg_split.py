#!/usr/bin/env python3
"""Bin-level ARGs from the assembly-wide DeepARG run (performance review M8).

usage: bin_arg_split.py <sample> <contig_bins.tsv> <proteins.faa[.gz]> <deeparg .mapping.ARG>

DeepARG scales its features over the whole input file, so a per-bin run gives
probabilities that depend on the bin (a lone candidate is scaled against
itself). The bin-level ARGs are therefore DEEPARG_CONTIGS' hits on the
assembly, split by the bin of each gene's contig (BIN_CONTIG_LABELS), into the
layout ARG_FASTA_FORMATTER reads:
  <sample>_bins_proteins/<sample>_<bin>_proteins.faa            the bin's Pyrodigal proteins
  <sample>_deeparg_results/<sample>_<bin>_deep_arg.out.mapping.ARG  DeepARG's header and the bin's hits
Every bin gets both files (a bin without hits: header only). A gene id is
<contig>_<gene number>.
"""
import gzip
import os
import sys


def contig_of(gene):
    return gene.rsplit('_', 1)[0]


def main():
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    sample, bins_tsv, faa, mapping = sys.argv[1:]
    contig_bin = {}
    with open(bins_tsv) as fh:
        header = fh.readline().rstrip('\n').split('\t')
        if header[:2] != ['contig', 'bin']:
            sys.exit(f"ERROR: {bins_tsv} is not a contig to bin table (contig, bin, ...)")
        for line in fh:
            fields = line.rstrip('\n').split('\t')
            if len(fields) >= 2 and fields[0]:
                contig_bin[fields[0]] = fields[1]
    bins = sorted(set(contig_bin.values()))

    prot_dir, arg_dir = f'{sample}_bins_proteins', f'{sample}_deeparg_results'
    os.makedirs(prot_dir)
    os.makedirs(arg_dir)

    opener = gzip.open if faa.endswith('.gz') else open
    outs = {b: open(os.path.join(prot_dir, f'{sample}_{b}_proteins.faa'), 'w') for b in bins}
    n_proteins = 0
    with opener(faa, 'rt') as fh:
        out = None
        for line in fh:
            if line.startswith('>'):
                out = outs.get(contig_bin.get(contig_of(line[1:].split()[0])))
                n_proteins += out is not None
            if out is not None:
                out.write(line)
    for out in outs.values():
        out.close()

    with open(mapping) as fh:
        arg_header = fh.readline()
        cols = arg_header.rstrip('\n').split('\t')
        if 'read_id' not in cols:
            sys.exit(f"ERROR: {mapping} has no read_id column (not a DeepARG .mapping.ARG file)")
        gene_col = cols.index('read_id')
        hits = {b: [] for b in bins}
        n_hits = 0
        for line in fh:
            b = contig_bin.get(contig_of(line.rstrip('\n').split('\t')[gene_col]))
            if b is not None:
                hits[b].append(line)
                n_hits += 1
    for b in bins:
        with open(os.path.join(arg_dir, f'{sample}_{b}_deep_arg.out.mapping.ARG'), 'w') as out:
            out.write(arg_header)
            out.writelines(hits[b])
    print(f"{sample}: {len(bins)} bins, {n_proteins} proteins, {n_hits} DeepARG hits in bins")


if __name__ == '__main__':
    main()
