# Reporting roadmap

Planned improvements to how BugBuster presents its results. The first round (October 2026) was driven by a real
32-sample run with QC, sylph and K-MARVEL. It produced the run-level summary (`00_summary/`,
`bin/results_summary.py`) and fixed the QC, taxonomy and classification plots (CHANGELOG, Unreleased).

Items that come from reading the code, not from real output, are marked **unverified** until example data confirm
them. In November 2026 (roadmap Phases 11a and 11b) the taxonomy, contig, binning, bin ARG and functional sections
were checked against reference runs of every branch (a per-sample assembly and a co-assembly of the same subsampled
real reads, and Kraken2 and sourmash runs on a 3-sample, 2 M-pair fixture); each item now says verified, fixed,
deferred or dropped, with the evidence. When a module gets example data:

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

## Taxonomy: Kraken2/Bracken and sourmash (verified 2026-11-01)

These paths share `bin/taxonomy_phyloseq.py` and `bin/taxonomy_report.py` with sylph. They already get the
plot fixes (top-k union, natural order, scaled labels), and `results_summary.py` reads their `*_otu_table.tsv` /
`*_tax_table.tsv` without changes.

- **Kraken2 OTU values are Bracken species estimates: verified.** Each sample's OTU column equals the sum of
  Bracken's `new_est_reads` (the `S` rows; Bracken writes no other level by default).
  - Open (not testable on the default level): `S` is hardcoded even when `--bracken_tax_level` is not `S`.
  - Reads assigned above species are outside the denominator: verified, small (the OTU table holds 97.4-97.7 % of
    the classified reads). Stated in the summary's methods notes.
- **Sourmash OTU values are abundance-weighted hashes: verified** (`unique_intersect_bp / scaled * average_abund`
  for every row). Open (not testable: every sample matched): samples with zero matches are dropped, not zero-filled.
- **`Classified` units differ by profiler: verified in the stage file** (Kraken2 0-100, sourmash and sylph 0-1).
  Handled downstream: `samples_qc.classified_fraction` is a fraction for every profiler, and `taxonomy_report.py`
  chooses the scale by profiler; the `max() > 1` item is **dropped** (fixed earlier). Making the stage column a
  fraction stays optional.
- **Different denominators: verified.** Kraken2's `Classified` is pre-Bracken, the OTU table post-Bracken.
- **Duplicate accessions: dropped.** Neither profiler writes one (Bracken one row per taxid, gather each genome once).
- **Summary work: fixed.** The unclassified fraction is in `samples_qc`, and the README methods notes and the units
  table say what "relative abundance" means per profiler and that the classified fractions are not comparable between
  profilers (30 % vs 9 % on the same reads).
- **Fixed (roadmap Phase 11b): a custom database was labelled as the registry default.** With `--custom_kraken_db`
  or `--custom_sourmash_db` the labels are now the database's own name.

## Read-level ARGs: ARGs-OAP (verified on a 3-sample, 2 M-pair fixture, 2026-10-27)

ARGS_OAP runs ARGs-OAP 4.0.0 stages one and two; `bin/args_oap_report.py` (ARGS_OAP_REPORT) publishes
`05_arg_prediction/read_level/args_oap_summary/args_oap_long.tsv` and `args_oap_cells.tsv`, and
`results_summary.py` reads them (`bin/summary_sarg.py`). KARGA and KARGVA were removed the same day (roadmap Phase
8), with `bin/arg_norm_report.py` and the summary's MEGARes tables.

- **Units: SARG copies.** A copy is a read weighted by aligned length / reference length (and 1/2 or 1/3 on a two-
  or three-component reference); copies per cell divides by ARGs-OAP's nCell. The summary reports copies per cell as
  the primary unit, plus copies per 16S and ppm (single reads) in `arg_sarg_long`, apart from the read-pair RPKM of
  K-MARVEL and RGI.
- **Singletons are added to the pairs.** ARGs-OAP treats `<id>_tmp_reads_S` as a separate sample; the report sums a
  sample's rows and recomputes stage two's formulas. On a sample without singletons the report equals stage two's
  own `normalized_cell`, `normalized_16S` and `ppm` tables at every level (checked).
- **Nothing to normalise is not a failure.** ARGS_OAP skips stage two when no metadata row has reads, 16S and cells,
  and accepts stage two's exit 2 when no read survives the SARG search; it writes the outcome to
  `stage_two_status.txt` (`ok`, `no_sarg_hit`, `not_normalised`), which the report carries, and the summary leaves a
  not-normalised sample blank rather than 0. A file stage one drops after an error (logged only as a warning) fails
  the task (retried with more memory when a process was killed for memory); an empty singleton file is skipped.
- **Not compared with CARD methods.** SARG names differ from CARD's; ARGs-OAP has its own section without tiers, as
  KARGA had. A SARG-to-tier map remains possible later.
- **Fixture numbers.** 14-17 SARG types and 64-74 subtypes per sample, 0.21-0.31 copies per cell in total; cell and
  16S estimates identical to the stage-one-only runs before Phase 8.
- **Open: `extracted.fa` and `blastout.txt`** (the stage-one candidate reads and the raw BLASTX hits, ~90 MB per
  full-depth sample) are published but only the filtered files are read.

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

## Contig-level ARG and taxonomy: DeepARG + Metabuli/BlobTools (verified 2026-11-01)

`bin/arg_contig_level_report.py` writes `Contig_tax_and_arg_prediction.tsv` (Metabuli replaced BLASTN in roadmap
Phase 9b). `bin/arg_blobplot.py` and `bin/blobplot.py` draw the blob plots.

- **Fixed and verified on real runs:** the table is tab-separated, and `ARG_blob_plot.png` shows only the
  ARG-carrying contigs.
- **A contig with several ARGs is plotted once per ARG: verified** (co-assembly: 301 hits on 228 contigs, 44
  contigs with more than one). Deferred to roadmap Phase 11d (figure redesign).
- **Blob plots pool all samples: verified** (the phylum plot draws every assembly in one panel and never uses its
  `Sample` column). Deferred to Phase 11d, with:
  - `blobplot.py` shows 8 phyla (`<=` threshold) and `no-hit` takes a colour slot: verified.
  - `arg_blobplot.py` keeps the top 9 ARG classes by summed contig length, not by count, and writes `MLS` as `Mls`
    (found 2026-11-01).
  - Percentages are of contigs or rows, not of length or coverage: verified.
- **Coverage is raw: verified** (BlobTools `bam0`; under co-assembly pooled over all samples). The summary's
  `arg_contig_long` says so, and `arg_contig_abundance` gives each sample's reads on the gene (featureCounts).
- **Fragile contig names: dropped** (fixed earlier: everything before the last `_` of the gene id).
- **Summary work: fixed (Phase 11b).** `arg_contig_long` (each hit with its contig's taxon at the lowest rank
  Metabuli reached, bin, MIMAG tier and geNomad class), `arg_contig_hosts` (genes per host genus) and
  `arg_contig_abundance`, presented as linkage, apart from the read-level candidate hosts. Found on the way: the
  `counts` column is always 1, `taxonomy_rank` takes 14 values (the summary collapses them), and 65 % of DeepARG hits
  are below 80 % identity (stated in the caveats).

## Binning (verified 2026-11-01)

The binning reports are `bin/bin_quality_report.py`, `bin/bin_tax_report.py` and, for co-assembly, `bin/bin_summary.py`.

- **Quality classes are not MIMAG: verified.** Of the 27 MAGs called High in the per-sample assembly reference
  run (45 MAGs), 3 are MIMAG high (most lack a complete rRNA set); its 4 "Low" MAGs are MIMAG medium
  (contamination 5-10 %).
  - **Fixed (Phase 11b):** the summary's `mags` table gives MIMAG tiers, with rRNA and tRNA genes from Bakta
    (`tier_basis`; without Bakta at most medium). The stage classes stay as they are (owner decision).
  - **Fixed:** the plot legend said "Completeness <50%, Contamination >5%" for Low while the code applies OR; it now
    says "or".
- **Missing dimensions: verified, superseded.** `Mag_quality_summary.csv` has no sample and `Bin_summary.csv` no
  quality class; the summary's `mags` table has both, in both assembly modes.
- **Taxonomy plot drops bins: verified in code, not exercised** (every reference MAG had a GTDB phylum). Open.
- **Coverage is an unweighted mean of per-contig mean depths: verified**, and only under co-assembly. The summary
  keeps it as published, says so, and adds `pct_sample_bases`; under per-sample assembly `mag_abundance` uses the
  length-weighted contig coverage of the contig table. The bowtie2 cross-bin item is not testable on the bwa-mem3
  reference. **Fixed:** the depth file's fields were padded with spaces.

## Bin-level ARG clustering (verified 2026-11-01)

- **Little is published: dropped as worded.** The per-bin DeepARG hits and proteins are published, and the ARG FASTAs
  were published too, but to a stray top-level `arg/` folder. **Fixed (Phase 11b):** they are in
  `05_arg_prediction/bin_level/arg_fasta/<sample>/`.
- **Headers: verified and documented** (`<gene id>_from_<sample>_<bin name>_gene_<ARG>_ARG_class_<class>`, e.g.
  `k141_92578_35_from_S1_binette_bin10_gene_ROSA_ARG_class_fosmidomycin`; ARG names contain `_`, so the summary
  splits on `_from_`, `_gene_` and `_ARG_class_`).
- **ID collisions: verified, harmless.** Two per-sample assemblies shared 85,711 contig names, but the headers carry
  the sample and bin, so every clustered sequence stays distinct.
- **Summary work: fixed (Phase 11b).** `arg_mag_long` (each ARG gene of a MAG and its cluster at 90-100 %) and
  `arg_cluster_sharing` (clusters found in more than one MAG, at 90 %). Within one co-assembly the identity levels
  add little (193 genes, 189 clusters at 90 %); between per-sample assemblies 21-28 clusters span samples.

## Functional annotation, contig and read branches (verified 2026-11-01)

The tables in `07_functional_annotation/summary/` follow the design doc. There are no figures.

- **Contig-branch descriptions are blank: verified.** **Fixed in the summary:** COG categories and eggNOG CAZy
  families are named (`AA10|Auxiliary Activities 10.` is split into family and name); KO, EC and Pfam stay without
  names (no name list ships with the pipeline; owner decision). dbCAN's calls are a separate `cazy_dbcan` table.
- **TPM is relative to assembled genes only: verified, and large** (14-84 % of a sample's reads were on predicted
  genes). **Fixed:** the summary reports `pct_reads_on_genes` (`assembly_read_use`, `samples_qc`) and states the
  TPM and CPGE definitions.
- **Mates may be counted separately: dropped (intended).** featureCounts counts reads (design doc 5.2), and the
  genome equivalents count both mates too.
- **Summary work: fixed (Phase 11b)** as tables: `function_contig_long` and `function_contig_<ontology>` (CPGE when
  MicrobeCensus ran, else TPM), the annotated fraction in `samples_qc`, never merged with the read branch; and
  read-level ARGs per genome equivalent (`arg_*_rpkg`). The heatmap is left to Phase 11d.

## Out of scope

MultiQC is deliberately not part of the pipeline (see CLAUDE.md). The run-level summary does not replace it, and
should not grow per-tool QC panels that MultiQC would provide.
