# Reporting roadmap

Planned improvements to how BugBuster presents its results. The first round (October 2026) was driven by a real
32-sample run with QC, sylph and K-MARVEL. It produced the run-level summary (`00_summary/`,
`bin/results_summary.py`) and fixed the QC, taxonomy and classification plots (CHANGELOG, Unreleased).

The modules below have **no example results yet**. Every item comes from reading the code, not from looking at real
output, and is marked **unverified** until example data confirm it. When a module gets example data:

1. Confirm or drop each item.
2. Fix what is confirmed.
3. Add the module to `results_summary.py`.
4. Update this file.

## Lessons from the first round (apply to every module)

These problems were found on real data in the Sylph and K-MARVEL outputs. They are the checklist for each module below.

- **Depth normalisation.** Clean depth varied ~550x between samples, so raw counts are not comparable. Normalise to clean read pairs (RPM, plus RPKM for genes), or per genome equivalent when MicrobeCensus ran.
- **State the units.** `Reads_report.csv` counts single reads (both mates), while K-MARVEL and sylph count read pairs.
- **Low-depth flag.** "Not detected" in a shallow sample is unknown, not absent (`NA_low_depth` in the summary).
- **Names, not IDs.** Name the rows of a table by taxon or gene, not by accession or taxid.
- **"Other" in bar plots.** Picking the global top N by mean abundance hides taxa that dominate a few samples. Use the union of each sample's top-k, and use heatmaps for detail.
- **Figure layout.** Use natural sample order (V2 before V10), height that scales with the sample count, and the legend outside the plot.
- **Mutation-dependent ARG hits.** Detecting the gene is not detecting resistance. Check the mutations against CARD `snps.txt` (or the database's equivalent).
- **Alleles share reads.** Close alleles of one gene take the same reads, so represent the gene by its best-supported allele and never sum them.
- **Cross-batch merge.** Large studies run in batches: a summary must outer-join taxa and genes across results directories.
- **Join keys.** Modules spell the sample key differently (`Id`, `sample`, `Sample`, `sample_id`). The summary should accept all of them.

## Taxonomy: Kraken2/Bracken and sourmash (unverified)

These paths share `bin/taxonomy_phyloseq.py` and `bin/taxonomy_report.py` with sylph. They already get the
plot fixes (top-k union, natural order, scaled labels), and `results_summary.py` reads their `*_otu_table.tsv` /
`*_tax_table.tsv` without changes. The remaining issues are:

- **Kraken2 OTU values are Bracken `clade_reads` for `S` rows only** (`taxonomy_phyloseq.py`, kraken parser).
  - `S` is hardcoded even when `--bracken_tax_level` is not `S`. Check what the table holds then.
  - Reads assigned above species are outside the denominator.
- **Sourmash OTU values are abundance-weighted hashes.**
  - The value is `unique_intersect_bp / scaled * average_abund`, not reads.
  - Samples with zero matches are dropped, not zero-filled.
- **`Classified` units differ by profiler.**
  - Kraken2 is 0-100; sourmash and sylph are 0-1. Sourmash's is the summed `f_unique_weighted`, a k-mer fraction, not reads.
  - The plot detects the scale with `max() > 1`, which misreads a kraken2 batch where every sample is at most 1% classified.
  - Make the column a fraction for every profiler.
- **Different denominators.** Kraken2's `Classified` comes from the pre-Bracken report, while the OTU table is post-Bracken.
- **Duplicate accessions.** The pivot uses the default `aggfunc='mean'`, so a duplicate accession within a sample would be averaged, not summed.
- **Summary work.** Pass the profiler's unclassified fraction into `samples_qc`. The summary must not assume sylph's `relative_abundance` semantics: say what "relative abundance" means per profiler in the summary report's methods section.

Example data needed: a kraken2 + Bracken run and a sourmash run on the same samples as a sylph run.

## Read-level ARGs: KARGA / KARGVA / ARGs-OAP (verified on Pudahuel, 3 samples, 2026-10-17)

`bin/arg_norm_report.py` publishes `05_arg_prediction/read_level/summary/KARGA_norm.csv` and `KARGVA_norm.csv`;
`results_summary.py` reads them (`bin/summary_megares.py`).

- **Units: k-mer depth, labelled.** `CPM_ARGs` and `copies_per_cell` divide `AverageKMerDepth`, a k-mer depth. The
  summary reports copies per cell as the primary unit and `cpm_kmer` as k-mer depth per million single reads, apart
  from the read-pair RPKM of K-MARVEL and RGI. Normalising read counts instead remains open.
- **Cells are not counted twice.** ARGs-OAP stage one writes one row per sample (R1 and R2 merged), so the
  `groupby().sum()` sees one row. Confirmed.
- **Alleles share k-mers.** Confirmed: the AAC(6') alleles of one sample had depths 46.8 / 51.1 / 47.6 / 51.7. The
  summary represents each MEGARes group by its deepest allele (`group_representative`). `Gen` is the MEGARes group.
- **No SNP confirmation for KARGA.** The summary carries MEGARes' `RequiresSNPConfirmation` flag as
  `snp_confirmation_required`; no Pudahuel hit had it. KARGVA has its own SNP-hit filter.
- **Absent samples.** A sample with no passing KARGVA hit has one row with blank ARG fields: the summary reads it as no
  hits. `nCell = 0` leaves a blank, kept blank in the matrices. A sample missing ARGs-OAP output is still dropped by
  `arg_norm_report.py` (unverified).
- **Not compared with CARD methods.** MEGARes names differ from CARD's; the owner chose a separate section without
  tiers (2026-10-17). A MEGARes-to-tier map remains possible later.
- **Open: ARGs-OAP `extracted.fa`** (~90 MB per sample, the stage-one candidate reads) is published but nothing reads
  it. Only stage one runs, so there is no SARG abundance.
- **Open: KARGVA `versions.yml`.** `modules/local/kargva/main.nf` indents its heredoc with tabs, so `END_VERSIONS` is
  left in the file and `software_versions.yml` is not valid YAML. The summary reads the versions line by line. Fixing
  the module reruns KARGVA (~4.5 h per sample) on `-resume`.

## Read functional (Woltka) and MicrobeCensus (verified on Pudahuel, 2026-10-17)

`results_summary.py` reads AGGREGATE_READ_FUNCTIONS' tables and MicrobeCensus (`bin/summary_functional.py`); the
tables match the pipeline's own CPGE matrices exactly.

- **Open: MetaCyc mixes reactions with pathways.** WoLr2's `function/metacyc/reaction-to-pathway.map` has 159 targets
  that are reaction IDs (`1.1.1.39-RXN`, ...), not pathways, and absent from `pathway_name.txt`.
  `bin/woltka_function_profile.py` follows the map, so 135 such IDs appear as MetaCyc "pathways" in a Pudahuel sample.
  Fix: keep only targets listed in `pathway_name.txt` (changes WOLTKA_CLASSIFY's output; ask first).
- The per-sample ORF table's sample column is `<sample>.wol` (the SAM name), unlike every other table.
- ARGs-OAP's cell estimate (`nCell`) and MicrobeCensus genome equivalents differ by 16-22% on the same samples; the
  summary reports both and labels them as different estimates.

## Software versions (found 2026-10-17)

- **Open: RGI_BWT / RGI_KMER** write multi-line junk (`rgi:` with no version, `kma --version` is not a KMA option).
  Fixing the modules reruns RGI_BWT on `-resume` (hours per sample).

## Read-level ARGs: RGI bwt / kmer (verified on batch1, 32 samples, 2026-10-16)

`RGI_REPORT` (inline Python in `modules/local/rgi_report/main.nf`) publishes `RGI_summary_report.csv` and four PNGs;
`RESULTS_SUMMARY` reads the RGI allele and k-mer gene tables.

Done:

- **RGI in the summary.** Genes count at `--rgi_summary_min_coverage` (80) % coverage of their best reference
  sequence (`rgi bwt` lists every gene with any mapped read: 38-548 per batch1 sample, 0-140 at 80%). Rows go into
  `arg_long.tsv` with a `method` column; RGI is the primary method when K-MARVEL did not run.
- **Pathogen of origin is parsed** from `<sample>_rgi_kmer_61mer_analysis.gene.txt` (`CARD*kmer Prediction`):
  `rgi_kmer_hosts` and `rgi_host_check` (top host a watchlist pathogen or not). On batch1 it separated real
  priority genes (KPC, NDM -> *K. pneumoniae*; vanB -> *E. faecium*) from WildCARD sequences of gut commensals
  (poxtA -> *Blautia producta*, *Enterocloster*; vanD in one sample -> *Enterocloster clostridioformis*).
- **K-MARVEL vs RGI.** Gene-level comparison (alleles collapsed): batch1 212 genes found by both, 47 K-MARVEL only
  (24 of them seen by RGI below the cutoff, 16 mutation-dependent genes RGI bwt does not test), 43 RGI only.
- **`RGI_REPORT` counts at the cutoff.** `*_cov<N>` columns beside the unfiltered ones; the plots use them.

Still open:

- **Fallback table headers (confirmed).** Zero-hit samples are reported with 0 genes, but the headers of
  `bin/rgi_bwt_lowcov_fallback.sh` are those of `rgi main` (16 and 27 columns), not RGI's `rgi bwt` tables
  (26 and 24). RGI_REPORT and the summary treat a table without hits as "no hits" whatever its header; anything else
  reading them would not.
- **vanG is tier 1 and RGI finds it in 26 of 32 batch1 samples**, mostly from WildCARD sequences of Clostridia
  (*Ruthenibacterium*, *Faecalibacterium*); K-MARVEL calls none. The host check flags them; decide whether vanG
  belongs in tier 1 at all for gut metagenomes.
- **Tier 4 for RGI.** `rgi bwt` reports homolog matches; mutation-dependent genes get `not_checked_rgi`.
- **RGI_REPORT plots pool all samples.** `value_counts()` counts sample x gene rows, and the per-sample bar is in
  alphabetical order. The summary's figures cover the per-sample view; consider dropping these plots.
- **No MAPQ filter.** Reads mapping equally well to several WildCARD sequences are counted for one of them by KMA.
- **KARGA and KARGVA are not in the comparison.** They use MEGARes; the summary reports them in their own section
  (see the KARGA section).

## Contig-level ARG and taxonomy: DeepARG + BLASTN/BlobTools (unverified)

`bin/arg_contig_level_report.py` writes `Contig_tax_and_arg_prediction.tsv`. `bin/arg_blobplot.py` and `bin/blobplot.py` draw the blob plots.

- **Fixed:** the table is now tab-separated, and `ARG_blob_plot.png` shows only the ARG-carrying contigs (contigs
  without an ARG were labelled "Other"). Both were confirmed on the test fixtures, not yet on a real run.
- **A contig with several ARGs is plotted once per ARG.** The table has one row per contig and ARG hit, and the blob
  plot draws every row.
- **Blob plots pool all samples.**
  - `blobplot.py` shows only 8 phyla (`<=` threshold), and lets `no-hit` take a colour slot.
  - Percentages are of contigs or rows, not of length or coverage.
- **Coverage is raw.** It is BlobTools `bam0`, not normalised to depth.
- **Fragile contig names.** `contig_name` is taken from the first two `_` tokens of the ORF id, which assumes MEGAHIT `k141_N` naming.
- **Summary work.** This is real ARG→taxon linkage: the ARG and the taxon are on the same contig. Present it as linkage, separate from the read-level co-occurrence hosts. Normalise coverage to clean read pairs.

Example data needed: an `--assembly_mode assembly --contig_tax_and_arg` run.

## Binning (unverified)

The binning reports are `bin/bin_quality_report.py`, `bin/bin_tax_report.py` and, for co-assembly, `bin/bin_summary.py`.

- **Quality classes are not MIMAG.**
  - Contamination above 5% is classed Low.
  - The plot legend says "Completeness <50%, Contamination >5%", but the code applies OR.
  - Use MIMAG tiers (high / medium / low, with tRNA/rRNA where available) and say so.
- **Missing dimensions.** `Mag_quality_summary.csv` has no sample dimension. `Bin_summary.csv` has no quality class and joins on `Name` alone.
- **Taxonomy plot drops bins.** Bins without a phylum are dropped, and every phylum is plotted (no top-N).
- **Coverage is an unweighted mean of per-contig mean depths, not normalised.** On the bowtie2 path, cross-bin reads also inflate it.
- **Summary work.** Add a MAG table (sample, bin, quality tier, completeness, contamination, GTDB lineage, coverage per clean read pair) and a MAG × sample abundance heatmap.

Example data needed: an assembly + binning run (≥2 binners for Binette refinement), and a co-assembly run.

## Bin-level ARG clustering (unverified)

- **Little is published.** Only the mmseqs `Arg_cluster_{90,95,99,100}_cluster.tsv` tables are published (representative→member, no header).
- **Unpublished inputs.** The ARG FASTAs from `ARG_FASTA_FORMATTER` are not published, and its `extract_sequences.py` lives in the container, not the repo, so the headers (sample and bin provenance) cannot be checked.
- **Possible ID collisions.** Prodigal ORF ids from different samples' MEGAHIT contigs may collide unless the headers carry the sample.
- **Summary work.**
  - Publish the FASTAs.
  - Document the headers.
  - Add a cluster × sample/bin presence table, showing which ARG clusters are shared across bins, samples or groups.

Example data needed: an `--arg_bin_clustering` run.

## Functional annotation, contig and read branches (unverified)

The tables in `07_functional_annotation/summary/` follow the design doc. There are no figures.

- **Contig-branch descriptions are blank.** The wide matrices are ID-only (KO, COG letter, EC, Pfam). Add names from the ontologies; the read branch already has them.
- **TPM is relative to assembled genes only** (unmapped reads are outside the denominator). CPGE is reads per kb per genome equivalent. Say so in a summary, next to the annotated fraction.
- **Mates may be counted separately.** featureCounts runs with `-p` but without `--countReadPairs`. Check whether this matches the design doc's intent.
- **Summary work.**
  - Add a functional section with the top functions per ontology as a heatmap (CPGE when available, TPM otherwise) and the annotated fraction per sample.
  - Keep read and contig branches separate (design doc rule).
  - When MicrobeCensus ran, give read-level ARGs a per-genome-equivalent normalisation (RPKG) too.

Example data needed: a `--contig_level_functional` run and a `--read_level_functional` run on the same samples.

## Out of scope

MultiQC is deliberately not part of the pipeline (see CLAUDE.md). The run-level summary does not replace it, and
should not grow per-tool QC panels that MultiQC would provide.
