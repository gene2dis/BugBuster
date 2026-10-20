# Changelog

All notable changes to BugBuster will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

The pipeline self-reports this state as `2.0.0dev` (manifest version) until the
next release is tagged.

### Added

- **Run-level results summary (`00_summary/`, `--results_summary`, on by default)**
  - `RESULTS_SUMMARY` (`bin/results_summary.py` and its `summary_*.py`
    modules) consolidates the QC, taxonomy and read-level ARG (K-MARVEL, RGI)
    results into tidy and wide TSV tables, an Excel workbook
    (README sheet with column definitions), PDF/PNG/SVG figures and a
    self-contained HTML report (inline SVG, sortable tables, no network)
  - Counts in read pairs (`Reads_report.csv` counts single reads); a
    low-depth flag (`--summary_min_clean_pairs`, default 1,000,000 clean
    pairs) turns "not detected" into `NA_low_depth` instead of absent
  - K-MARVEL abundance as RPM and RPKM of clean read pairs (the
    `kmarvel_summary/` matrices are raw counts); each gene is represented by
    its best-supported allele, as alleles of one gene share their reads
  - ARG reporting tiers from an editable rule table
    (`assets/arg_priority_tiers.tsv`, `--arg_priority_table`): 1 clinical
    priority, 2 other acquired, 3 intrinsic/efflux/regulator, 4 mutation
    dependent; tier 4 hits are checked against the K-MARVEL database's CARD
    `snps.txt` (`resistance_snp_status`), since K-MARVEL lists every
    difference from the reference, not resistance mutations
  - Taxonomy tables named by taxon per rank, a long table with sylph ANI,
    ANI interval and coverage (`species_call_uncertain`), species prevalence,
    and a watchlist of priority pathogens (`assets/summary_watchlist.tsv`,
    `--summary_watchlist`)
  - Optional `--sample_metadata` with `--summary_group_by` /
    `--summary_order_by`: grouped and ordered figures, per-group ARG
    detection patterns
  - RGI in the summary: `rgi bwt` genes count at `--rgi_summary_min_coverage`
    (default 80) % coverage of their best reference sequence (RGI lists
    every gene with a mapped read), with RGI's k-mer pathogen of origin
    (`rgi_kmer_hosts`) and whether the top host is a watchlist pathogen
    (`rgi_host_check`: WildCARD holds resistance-gene sequences of gut
    commensals, e.g. poxtA from *Blautia*); RGI is the primary ARG method
    when K-MARVEL did not run
  - K-MARVEL vs RGI comparison when both ran: gene-level overlap tables
    (`arg_method_comparison`, `arg_method_overlap`, with RGI calls below the
    cutoff marked), an UpSet plot and agreement by tier and category
  - QC figure marks 10M and 20M clean read pairs (`--summary_depth_marks`)
    next to the low-depth line, and writes the % classified on each bar; the
    report header counts the samples at each mark (`depth_class` column)
  - Taxonomy heatmaps show each sample's 5 most abundant taxa
    (`--summary_taxa_per_sample`, no row cap; `--summary_top_n` now caps the
    ARG heatmap only), grouped by phylum and alphabetical within each; group
    names sit in their own strip, so they no longer overlap the title or each
    other
  - An "Attention" list at the top of the HTML report (tier 1 genes found by
    every method, by one only or with a commensal RGI host, and low-depth
    samples), and a first "Reading this report" section explaining the units,
    the ARG tiers (from the tier rule file in use) and the ARG methods
  - The watchlist notes name each taxon's source (ESKAPE, WHO Bacterial
    Priority Pathogens List 2024, CDC AR Threats Report 2019), and the
    report's Methods section shows the list
  - Optional `--card_prevalence` (WildCARD): candidate ARG hosts among the
    species of the same sample (co-occurrence, not linkage) and CARD plasmid
    association
  - `batch` column: the output directory's name, or `--summary_batch_label`
    (script: `--batch-labels`) when batches share an output directory name
  - The script also runs on finished results directories, merging several
    batches without rerunning any tool (manual, Example 13); pinned Seqera
    image with pandas, matplotlib and openpyxl
  - `docs/REPORTING_ROADMAP.md`: reporting improvements planned for the
    modules without example results yet
  - KARGA / KARGVA in the summary (ARG_NORM_REPORT's tables): long tables,
    copies per cell per MEGARes group and class (each group's deepest allele,
    as alleles share k-mers), the KARGVA variants, a class heatmap, and the
    ARGs-OAP 16S and cell estimates in `samples_qc`; reported apart from the
    CARD methods (MEGARes names differ) and in their own unit
  - Read functional profile (Woltka, SUPER-FOCUS or HUMAnN) and MicrobeCensus
    in the summary: `function_read_long`, one wide table per ontology in CPGE
    (`NA_no_ags` without a MicrobeCensus estimate), a COG-category heatmap,
    and reads assigned, annotated fraction, average genome size and genome
    equivalents in `samples_qc`
  - `README.md` and `README.docx` (standard-library writer, no new
    dependency) and `LLM_GUIDE.md` in `00_summary/`: the stages that ran per
    results directory (from `software_versions.yml`, read line by line as it
    is not always valid YAML), every table and figure, units, tiers, caveats,
    versions; the guide covers the stages that ran, methods, every table
    column and figure, pitfalls and recipes. Both describe the summary
    folder only (the part that is shared): paths are relative to it, batches
    appear by label and input files by name, never as local paths
  - One generic input (`stage_reports`) lays out other stages' run-level
    tables by file name, so a stage is added with one case in the module and
    one `mix()` in `main.nf`
  - Every emitted column is documented (`bin/summary_docs.py`); the test
    suite fails on an undocumented one, and the Excel README sheet lists only
    the tables written
  - Table audit on a real run: RbpA and the vanO cluster (chromosomal in
    Actinomycetota, 74-82% of acquired RPKM in one batch) move to tier 3;
    a gene gets one tier across ARG methods (RGI tiers also use CARD
    `snps.txt`); an allele name that collapses onto another CARD gene keeps
    its own name (CRP-1, a CRP beta-lactamase, is not CRP, an efflux
    regulator); `rgi_gene_read_pairs` is written, so RGI gene RPKM can be
    rebuilt from the tables; taxonomy wide and watchlist tables mark taxa
    not found in a low-depth sample `NA_low_depth`, and species prevalence
    counts the samples with taxonomy output and leaves low-depth absences
    out of its denominator; gene family and drug class
    lists use one separator; the report says where tier 1 genes found only
    by another method are; wide tables are filled in one step, not one
    column per sample (pandas PerformanceWarning on 100+ samples)

- **Sylph read-level taxonomic profiler (`--taxonomic_profiler sylph`)**
  - `SYLPH_PROFILE` (sylph 1.0.0) profiles the paired host-removed reads by
    containment ANI against a sylph database with `-u` (unknown-fraction
    estimation); `SYLPH_TAXPROF` (sylph-tax 1.9.2, `--no-config` with the
    staged metadata file) turns each genome profile into a lineage profile
    (`.sylphmpa`), writing a header-only profile for a sample without genomes
  - New `--sylph_db` registry (`gtdb_r232` default: GTDB r232 c200 two-stage
    `.syl2db`, 19 GB; `gtdb_r226`, `gtdb_r220`), each downloaded with its
    sylph-tax metadata by `FORMAT_SYLPH_DB` / `FORMAT_SYLPH_TAXONOMY`
    (`bin/sylph_db_reformat.sh`) into `<databases_dir>/sylph/`;
    `--custom_sylph_db` and `--custom_sylph_taxonomy` (required together,
    checked at launch, so a database is never paired with another release's
    taxonomy) skip both downloads; outputs are then labelled with the custom
    database's file name (e.g. `gtdb-r232-c200-dbv2`), not the `--sylph_db` key
  - `--sylph_read_seq_id` passes `--read-seq-id` to sylph: its automatic
    read-accuracy estimate, which drives the unknown fraction (Classified and
    the estimated read pairs), is skewed by reads from genomes absent from the
    database (simulated 1/3 unknown: 100% classified automatically, 66.1%
    with the identity set, truth 66.7%)
  - `SYLPH_PROFILE` requests 12 GB with a two-stage `.syl2db` (GTDB r232:
    ~4.3 GB peak measured, it loads only the sparse first-stage index) and
    36 GB with a single-stage `.syldb`, doubling on an OOM retry; the read-pair
    count runs concurrently with sylph
  - Outputs: per-sample profiles in `02_taxonomy/sylph/`; the phyloseq tables
    hold sylph's relative (taxonomic) abundance per genome, and a new
    `sylph_<db>_read_counts.tsv` holds the estimated read pairs per genome
    (`sequence_abundance` × the read pairs given to sylph); `Reads_report.csv`
    gains the `Sylph DB` / classified fraction columns. `TAXONOMY_PHYLOSEQ`
    gains an optional `read_counts` output
  - Tests: real-container nf-tests for both modules on a committed two-genome
    database and simulated reads (`tests/bin/make_sylph_fixtures.sh`), sylph
    cases in the report / phyloseq / DB reformat script harnesses and two
    pipeline stub runs (registry and custom databases)

- **Read-level ARG detection with K-MARVEL (`--kmarvel_prediction`)**
  - New optional read-level step, independent of assembly and of
    `--read_arg_prediction` / `--rgi_prediction`: `KMARVEL` runs
    [K-MARVEL](https://bitbucket.org/amr-avengers/k-marvel) 1.0.0 (commit
    7108d2d) on each sample's host-removed reads (R1, R2 and singletons
    pooled, short-read mode) against the CARD protein models: assembled ARGs
    are chimera-filtered, variant models need a CARD resistance SNP, and
    every ARG is confirmed with blastp
  - `KMARVEL_REPORT` (`bin/kmarvel_report.py`) merges the per-sample tables
    into a long table, ARG (CARD ARO) x sample and AMR gene family x sample
    matrices of supporting read pairs (both mates count once), and a drug
    class x sample ARG-count matrix. The gene-family matrix is the one to
    compare samples with: the allele K-MARVEL reports among close variants
    of one gene is not stable (`05_arg_prediction/read_level/kmarvel/`
    and `kmarvel_summary/`); samples without ARGs stay as zero columns
  - `FORMAT_KMARVEL_DB` (`bin/kmarvel_db_reformat.sh`) builds the database
    from the md5-verified CARD 4.0.1 data tarball (~4.6 MB, `--kmarvel_db
    card_v4.0.1`); `--custom_kmarvel_db` takes an existing one
  - New parameters `--kmarvel_min_coverage`, `--kmarvel_min_identity` (90 %
    each) and `--kmarvel_min_read_length` (100 bp)
  - K-MARVEL has no bioconda package: the pipeline builds its own image
    (`containers/kmarvel/Dockerfile`: K-MARVEL built with cargo 1.93.1 against
    a refreshed `Cargo.lock`, since upstream's is stale; BLAST+ 2.17.0),
    published to `ghcr.io/gene2dis/bugbuster-kmarvel` by the manual
    `.github/workflows/kmarvel-image.yml`
  - Works around three K-MARVEL 1.0.0 behaviours: it reads only the first
    member of a multi-member gzip (the reads, gzipped or plain, are
    streamed through a named pipe in a local temp dir instead of
    concatenating `.gz` files; no uncompressed copy on disk, and a task
    whose reads K-MARVEL never consumed fails instead of hanging), it exits 0 when its blastp step
    fails and writes a header-less table for zero ARGs
    (`bin/kmarvel_check_output.sh` fails the task, resp. writes the header),
    and it cuts output names at the first `.` (it runs on fixed file names)

- **bwa-mem3 short-read aligner (`--read_aligner`, now the default; see Changed)**
  - Covers the pipeline's own alignment steps: host/phiX decontamination (`BWAMEM3_BUILD_COMBINED`, `BWAMEM3_DECONTAMINATE`), reads-to-contigs mapping for binning, contig coverage and gene counting (`BWAMEM3_INDEX`, `BWAMEM3_SAMTOOLS`, `BWAMEM3_SAMTOOLS_PER_SAMPLE`), and co-assembly bin depth (`BINS_REFERENCE`, `BWAMEM3_BINS_INDEX`, `BWAMEM3_BIN_DEPTH`, `BEDTOOLS_BIN_DEPTH`). Woltka and HUMAnN/MetaPhlAn keep bowtie2 with their own databases; RGI keeps its own aligner (KMA by default, `--rgi_aligner`)
  - bwa-mem3 0.14.0 and samtools 1.22.1 run from one pinned Seqera Containers image. Output is independent of the thread count (`-K` fixed)
  - Decontamination is strict: a pair is clean only when both mates are unmapped. bowtie2 `--un-conc` kept pairs with a single host mate. Clean-read counts are therefore not directly comparable with bowtie2 runs, and the FASTQ comment field is not carried over. Output file names are set explicitly, so sample ids containing dots no longer break decontamination. Mates are aligned without pairing (`-P`) so that `-T` holds for each mate: in paired mode bwa mem skips it for a chosen pair, and random ~20 bp host hits on overlapping mates removed 17.5% of a bacterial mock community (0.1% with `-P`). New `--decontam_min_score` (bwa `-T`, default 30) and `--bwamem3_dense_sa` (denser suffix array: faster mapping for ~12 GB more RAM)
  - The decontamination index is published to `<databases_dir>/bwamem3_index/`. A human host needs ~75 GB RAM to build it (minutes rather than bowtie2-build's hours). The request is sized from the genome FASTA files (~24 bytes per base, a gzipped file counted at 4x its size; ~90 GB for the default human genome, 8 GB minimum), capped by `--max_memory`. On bwa-mem3's out-of-budget exit status 3 this build and the contig and bin index builds retry with more memory only while the request is below `--max_memory`; at the cap the task fails (the decontamination build with a message naming the three ways out). `--custom_decontamination_index` accepts a bwa-mem3 or bwa-mem2 index under `bwamem3` and rejects a Bowtie2 one at launch
  - Contig alignments drop secondary and supplementary records, keeping one record per read as with bowtie2. Short or divergent local hits become unmapped reads: an alignment needs >= 95% identity (`--bwamem3_min_identity`) and >= 50 aligned bases (`--bwamem3_min_aligned_length`), or 80% of a read shorter than that (2x50 libraries, reads trimmed below 50 bp). Without the filter, reads from genomes missing from the assembly inflated gene counts (Spearman vs bowtie2 0.74 in a poorly covered sample) and the depth of near-absent bins (~1.9x). Each library of a pooled co-assembly is mapped separately and the results merged. The co-assembly index is built once instead of once per alignment task
  - Co-assembly bin depth maps each sample once against all refined bins concatenated, instead of once per bin (bwa mem aligns locally, so per-bin mapping would let other bins' reads inflate coverage). It writes the same `*_bin_depth.tsv` table; a refined bin that adds no contig to the reference (no sequences, or only contigs of an earlier bin) still gets a 0-coverage row
  - The bioconda bwa-mem3 build needs an AVX2-capable CPU (x86-64-v3)
  - `--read_aligner` is listed in `--help`

- **Read-level functional profiling (`--read_level_functional woltka`)**
  - New optional `READ_FUNCTIONAL` subworkflow, independent of assembly (runs
    with `--assembly_mode none`) and of the contig/MAG branches; one backend
    per run, validated at launch
  - Woltka 0.1.7 backend: Bowtie2 alignment of the host-removed reads against
    the Web of Life release 2 genomes (WoLr2) with the SHOGUN multi-hit
    settings (`WOLTKA_ALIGN`, trimmed SAM), Woltka ORF classification
    (`WOLTKA_CLASSIFY`; mates counted separately, multi-hit reads divided 1/k,
    or left unassigned with `--woltka_uniq`), and per-ORF composition to KO,
    EC, COG-category, Pfam and MetaCyc pathway read counts — each ORF counted once per
    distinct term (Woltka's own `collapse` was found to inflate counts through
    duplicate map entries and chained maps, so it is not used)
  - Study-level `AGGREGATE_READ_FUNCTIONS` writes `read_function_abundance.tsv`
    (same schema as the contig branch, `source = reads`, native unit = reads,
    copies per genome equivalent from MicrobeCensus), wide
    `read_function_wide_<ontology>_{native,cpge}.tsv` matrices,
    `read_annotated_fraction.tsv` and `read_sample_summary.tsv` to
    `07_functional_annotation/summary/`; per-sample tables to
    `07_functional_annotation/reads/woltka/<sample>/`. Reported separately
    from, never merged with, the assembly-based tables. Inputs are collected
    in a fixed (sorted) order, so resumed runs are fully cached
  - Same `cog` / `pfam` vocabularies as the contig branch: WoLr2's COG ortholog
    ids are mapped to COG functional categories with NCBI's `cog-24.def.tab`
    (downloaded for `--read_level_functional woltka` too; letters counted once
    per ORF), and Pfam families are reported by name, with the versioned
    accession in `description`. Five malformed COG strings in the WoLr2
    `ko-to-cog` map (e.g. `COG:1140`) and four ids absent from COG 2024 are
    skipped with a warning; any other unmapped id fails loudly
  - WoLr2 subset (~94 GB: Bowtie2 index + ORF coordinates + KEGG / MetaCyc /
    Pfam maps) auto-downloaded from the official host with md5 verification
    (`--woltka_db wolr2`), or a local mirror via `--custom_woltka_db`; the
    alignment needs ≥ 68 GB RAM
  - MicrobeCensus now also runs for the read branch (once per sample when both
    branches are on)
- **SUPER-FOCUS read-level backend (`--read_level_functional superfocus`)**
  - SUPER-FOCUS 1.8 (`SUPERFOCUS`, pinned Seqera Containers image with
    DIAMOND 2.2.1 and MMseqs2 18.8cc5c; docker works) against the DB_90 SEED
    subsystem cluster database; search backend via `--superfocus_aligner
    diamond|mmseqs2` (default `diamond`); R1, R2 and singleton reads
    concatenated into one query (mates counted separately), each read with a
    hit contributing 1, divided 1/k across its best-hit SEED assignments
  - SEED subsystem levels 1-3 summed by the pipeline from the function-level
    counts (ontologies `seed_level1`, `seed_level2`, `seed_level3`;
    path-qualified accessions `L1` / `L1 | L2` / `L1 | L2 | L3`) — SUPER-FOCUS's
    own per-level files merge the level-2 `-` placeholder across level-1
    categories and are kept only as raw outputs. SEED is never mapped to KO/EC
  - No copies per genome equivalent for this backend (a SEED hit has no gene
    length): `abundance_cpge` blank, `cpge_status = not_applicable`, only the
    `read_function_wide_seed_level{1,2,3}_native.tsv` matrices; AGS/GE still
    reported. Per-sample tables to `07_functional_annotation/reads/superfocus/<sample>/`
  - DB_90 archive for the selected aligner only (~0.74 GB DIAMOND / ~0.9 GB
    MMseqs2, figshare CC0, md5-verified, plus `database_PKs.txt` from the
    SUPER-FOCUS v1.8 tag; `FORMAT_SUPERFOCUS_DB`, `--superfocus_db db90`), or
    a local database root via `--custom_superfocus_db`
  - `bin/aggregate_read_functions.py` generalized to per-backend settings
    (file suffixes, versions keys, verified versions, ontologies, CPGE
    applicability); Woltka outputs unchanged (byte-identical on the fixtures)
- **HUMAnN 4 read-level backend (`--read_level_functional humann`)**
  - HUMAnN 4.0.0a2 (an alpha, pinned exactly) with MetaPhlAn 4.1.2 and the
    MetaPhlAn database `mpa_vOct22_CHOCOPhlAnSGB_202403` — the only
    combination this HUMAnN build accepts — in a pinned Seqera Containers
    image (`HUMANN`; no HUMAnN 4 conda package or biocontainer exists; docker
    works). R1, R2 and singleton reads concatenated into one input; run with
    `--count-normalization RPKs`
  - Ontologies `metacyc` (pathway abundance), `ko` and `ec` (gene families
    regrouped by the pipeline with the HUMAnN 4 mapping files; UniClust90
    families reach EC but never KO; KO coverage is limited by HUMAnN 4.0.0a2's
    legacy KO map, ~17-18 % of gene-family abundance on the test data vs
    ~39-61 % for EC). `abundance_native` = RPK
    (`native_unit = rpk`), CPGE = RPK / genome equivalents; six wide
    matrices. Per-sample raw HUMAnN tables (MetaPhlAn profile, gene families,
    reactions, pathway abundance, log) to `07_functional_annotation/reads/humann/<sample>/`
  - Version-aware parsing (`bin/humann_function_profile.py`): the exact
    4.0.0a2 table layout is required, anything else fails with an explanatory
    error. Works around two 4.0.0a2 defects: its own regroup / renorm
    utilities mishandle the `READS_UNMAPPED` row (the pipeline regroups
    itself), and sample names containing `s__` or `t__` crash its
    MetaPhlAn-profile parsing (rejected at launch)
  - Databases (`FORMAT_HUMANN_DB`, `--humann_db v4_alpha-full` (default,
    ChocoPhlAn 44.8 GB) or `v4_alpha-ec_filtered` (6.9 GB), recorded in
    provenance): ChocoPhlAn, UniRef90 EC-filtered (0.94 GB), the v4 utility
    mapping (2.8 GB) and the MetaPhlAn database with its prebuilt Bowtie2
    index (~22.5 GB, md5-verified) — ~71 GB / ~33 GB in total; or a local
    database root via `--custom_humann_db`
  - `bin/aggregate_read_functions.py`: per-backend native unit and
    annotated-fraction mode; Woltka and SUPER-FOCUS outputs unchanged
    (byte-identical on the fixtures)
- **MAG-level functional annotation (`--mag_level_functional`)**
  - Bakta 1.12.1 annotation of every MetaWRAP-refined bin, one task per bin
    (nf-core `bakta/bakta` module, patched to emit `versions.yml`; per-bin
    fan-out of the refined-bins directory), published to
    `07_functional_annotation/mags/<sample>/` as `<sample>_<bin>.{gff3,gbff,faa,fna,tsv,txt,hypotheticals.tsv,hypotheticals.faa}`
  - Requires `--include_binning` and at least two `--binners` (validated at
    launch): the MetaWRAP completeness/contamination filter only runs with
    ≥2 binners, and only quality-filtered bins are annotated. Bins are
    expected to be bacterial; the pipeline does not exclude
    archaeal/eukaryotic/viral bins (documented)
  - Bakta database v6.0 (schema 6) auto-downloaded from the pinned Zenodo
    release: `--bakta_db v6.0-full` (default, 31.9 GB download) or
    `v6.0-light` (1.3 GB; explicit choice, recorded in provenance via the
    `bakta_db` versions entry) / `--custom_bakta_db`. Provisioning refreshes
    the bundled AMRFinderPlus database with `amrfinder_update` (~300 MB; the
    official tarball's copy is too old for the AMRFinderPlus in the pinned
    bakta container and would fail every annotation at the AMR expert step —
    custom databases need the same one-time refresh, see troubleshooting)
  - Independent of `--contig_level_functional`; runs on any container engine
- **Contig-level functional annotation (`--contig_level_functional`)**
  - Shared Pyrodigal gene calling on contigs (one pass feeds both DeepARG and
    functional annotation), published to `07_functional_annotation/gene_calling/`
  - Two-stage eggNOG-mapper v3 annotation (DIAMOND search + orthology-transfer
    annotation) against the eggNOG 7 database (~44 GB, auto-downloaded;
    `--eggnog_db` / `--custom_eggnog_db`), published to
    `07_functional_annotation/eggnog/`
  - Per-gene abundance quantification with featureCounts over the Pyrodigal
    gene coordinates, with gene ids matching the annotated protein ids;
    per-sample counts in both assembly modes (under co-assembly, dedicated
    per-sample alignments against the co-assembly are added since the pooled
    binning BAM cannot yield per-sample counts), published to
    `07_functional_annotation/gene_abundance/`; multi-mapping policy via
    `--featurecounts_multimap` (`primary`/`all`/`none`)
  - Study-level aggregation (`bin/aggregate_functions.py`) into canonical
    tables published to `07_functional_annotation/summary/`: per-gene
    abundance with TPM and copies per genome equivalent (CPGE), long-format
    gene annotations, per-ontology function abundance (KO, COG, EC, Pfam,
    CAZy; intentional double-counting of multi-term genes), wide TPM and
    CPGE matrices per ontology, a per-sample annotated-fraction report and
    an AGS summary (`ags_and_ge.tsv`); the eggNOG annotations parser is
    version-aware and fails loudly on layout drift. Inputs are collected in a
    fixed (sorted) order, so resumed runs are fully cached
  - COG and Pfam terms: eggNOG 7 writes a COG ortholog id (e.g. `COG1629`)
    into `COG_category` for most genes (75-80 % on real data); these are
    mapped to COG functional-category letters with NCBI's COG definitions
    table (`cog-24.def.tab`, ~410 KB, auto-downloaded with md5 verification;
    `--cog_db` / `--custom_cog_db`). eggNOG-mapper v3 writes `PFAMs` values as
    `<pfam_name>_<start>_<end>`; the domain coordinates are stripped, so each
    Pfam family is reported by name and counted once per gene (format
    verified on all 73.9 M Pfam values of the eggNOG 7 database). An unknown
    COG id or an unexpected `PFAMs` form fails loudly
  - eggNOG-mapper 3.0.0-beta6 truncates OG names whose family part contains
    `|` (upstream issue #620; ~2 % of eggNOG 7 OGs), which cost those genes
    their COG category or gave them another OG's (1-4 % of annotated genes
    on real data; KO, EC, Pfam and CAZy unaffected). The annotate step runs
    through `bin/emapper_ogs_fix.py`, which applies the fix proposed on the
    issue at run time and refuses to run on any eggNOG-mapper version other
    than the beta6 it targets
  - MicrobeCensus average genome size estimation on the host-removed reads
    (`--microbecensus`, on by default with the branch), published to
    `07_functional_annotation/microbecensus/`, enabling the CPGE
    normalization. Failure is non-fatal by design: affected samples fall
    back to TPM-only with empty `cpge` fields and `status = unavailable` in
    `ags_and_ge.tsv`. Both published biocontainers of the unmaintained
    upstream tool are broken (Python-3 build: unreleased upstream str/bytes
    fix; Python-2 build: missing libstdc++ for the bundled RAPsearch2); the
    pipeline pins the Python-3 image and routes the call through the
    `bin/run_microbe_census_py3fix.py` shim that patches the one broken
    function
  - run_dbcan v5 CAZy annotation of the predicted proteins
    (`--functional_cazy`, on by default with the branch): protein-mode
    `CAZyme_annotation` (DIAMOND vs CAZy + pyHMMER vs dbCAN/dbCAN-sub HMMs)
    against the pinned dbCAN database release (~7.4 GB, auto-downloaded;
    `--dbcan_db` / `--custom_dbcan_db`), published to
    `07_functional_annotation/dbcan/` with the per-tool overview columns and
    the dbCAN-sub substrate predictions retained. The calls feed the summary
    tables as `db = dbcan` / `backend = run_dbcan` rows alongside the
    eggNOG-derived CAZy calls — reported separately, never merged — plus
    dedicated `function_wide_cazy_dbcan_{tpm,cpge}.tsv` matrices and
    `cazy_dbcan` annotated-fraction rows; the consensus policy is a
    documented parameter (`--dbcan_consensus`: `recommended` = calls
    supported by >= 2 tools, `any` = per-tool union), and the overview
    parser is version-aware like the eggNOG one
  - Runs on every container engine. eggNOG-mapper 3.0.0-beta6 is published
    only as an Apptainer image, which singularity/apptainer use; docker/podman
    and the cloud profiles use `ghcr.io/gene2dis/bugbuster-eggnog-mapper`, an
    image of the same version built from the authors' recipe with the same
    tool versions (`containers/eggnog-mapper/Dockerfile`, published by the
    `eggnog-image` workflow; same annotations as the official image). It is
    replaced by the bioconda image once v3.0.0 final is released
  - `docs/WORKFLOW_DIAGRAM.md` refreshed to cover the functional annotation
    branch (shared Pyrodigal gene calling, eggNOG/dbCAN database preparation,
    per-sample co-assembly counting alignments, the FUNCTIONAL_ANNOTATION
    subworkflow and its parameters)

### Changed

- **`software_versions.yml` also collects the `versions` topic**: current
  nf-core modules report `[process, tool, version]` tuples on the `versions`
  topic instead of a `versions.yml` file. A new `SOFTWARE_VERSIONS`
  subworkflow writes them in the same format as the local modules' files
  (one block per process), so the file and `RESULTS_SUMMARY` are unchanged.
  No module emits the topic yet; the nf-core updates (issue #35) will

- **MetaWRAP requests 80 GB on the first attempt** (was 72 GB from
  `process_high`): its peak RSS was 67-69 GB in every binning run, within 5 %
  of the old request. Still doubled on retry and capped by `--max_memory`

- **CPU requests sized to measured use** (`config/modules.config`):
  `SYLPH_PROFILE` 8 → 4, `MICROBECENSUS` 8 → 2, `WOLTKA_CLASSIFY` 8 → 2,
  `TAXONOMY_PHYLOSEQ` 4 → 2 (each used well under its request in real runs),
  and `CHECKM2_BATCH` 8 → 16 (one task over every bin at the end of a binning
  run, 585-673 % CPU of 8). Outputs are unchanged; more tasks fit at once

- **The execution trace records requested resources**:
  `pipeline_info/execution_trace_*.txt` adds the `start`, `complete`, `%mem`,
  `cpus`, `memory`, `time`, `attempt` and `workdir` columns, so CPU and memory
  requests can be compared with what each task used

- **`SOURMASH` no longer caps itself at 24 parallel tasks** (`maxForks 24`
  removed): the executor's CPU and memory pool already limits concurrency

- **`RGI_REPORT` adds coverage-filtered counts**: `rgi bwt` lists every gene
  with at least one mapped read (38-548 per sample in a 32-sample run), so
  `RGI_summary_report.csv` gains `AMR_Genes_cov<N>`, `Mapped_Reads_cov<N>`,
  `Unique_Drug_Classes_cov<N>` and `Unique_AMR_Families_cov<N>`, counting the
  genes whose best-covered reference sequence reads cover at least
  `--rgi_summary_min_coverage` % (default 80). The existing columns are
  unchanged; the four plots now use the filtered genes and say so in their
  titles

- **Breaking (default aligner)**: `--read_aligner` defaults to `bwamem3` instead of
  `bowtie2`, after a benchmark on five real metagenomes and a spiked mock community
  (`docs/benchmarks/bwamem3_vs_bowtie2.md`): host removal 4.1x faster, reads-to-contigs
  mapping 3.9-6.2x faster with the alignment filter, co-assembly bin depth ~34x faster,
  no human pairs leaked from the spike-in (bowtie2 0.157%) and 0.11% of the mock
  removed (bowtie2 3.0%), the same bins and GTDB species. `--read_aligner bowtie2`
  stays as the fallback, unchanged. What changes for an existing setup:
  - building the human decontamination index needs ~75 GB RAM (bowtie2: ~6 GB) and
    mapping peaked at ~16 GB per task (~4 GB). The pipeline warns at launch when
    `--max_memory` is lower (for the built-in `--host_db`; not for `--custom_host_fasta`)
  - the index is published to `<databases_dir>/bwamem3_index/`; a Bowtie2
    `--custom_decontamination_index` is rejected at launch, with a message naming
    `--read_aligner bowtie2`
  - clean-read counts are lower than in bowtie2 runs (a pair with one host mate is
    removed), and `-resume` of an earlier project re-runs from decontamination
  - the CPU needs AVX2 (x86-64-v3); without it the `BWAMEM3_*` tasks stop with
    "Illegal instruction" (exit 132)
  - the `--bowtie_*` options apply only under `--read_aligner bowtie2`

- **Breaking (Nextflow version)**: minimum Nextflow raised again, from
  24.04.0 (v1.1) to **25.10.0**, and nf-schema from 2.4.2 to 2.7.3. The code now
  follows Nextflow's strict syntax, whose parser needs nf-schema >= 2.7.2,
  which in turn needs Nextflow >= 25.10. nf-schema 2.8.0 needs Nextflow >= 26.04,
  so the plugin stays pinned at 2.7.3
- **Breaking (GTDB-Tk databases)**: GTDB-Tk re-pinned 2.5.2 → 2.7.2 and the
  reference-data registry moved from GTDB R220 to **R232** (~61 GB download,
  pinned release URL; default `--gtdbtk_db release_232`). GTDB-Tk 2.7.x
  accepts only R232 data, so a `databases/gtdbtk/` directory cached by earlier
  pipeline versions (R220) no longer works — delete it to re-download, or pass
  an R232 directory with `--custom_gtdbtk_db`. The upstream-removed
  `--skip_ani_screen` flag was dropped from the `classify_wf` invocation (the
  ANI pre-screen now always runs, against the skani DB bundled in the R232
  package), and the module's synthetic empty-report headers were updated to
  the 2.7 column schema (`fastani_*` → `closest_genome_*`;
  `bin_tax_report.py` already handled both namings)
- **Breaking (output layout)**: contig gene calling switched from the nf-core
  Prodigal module to a shared Pyrodigal step; its outputs moved from
  `05_arg_prediction/contig_level/prodigal/` to
  `07_functional_annotation/gene_calling/` (same predictions, gzipped, now
  produced once for both DeepARG and functional annotation)
- **Breaking (output layout)**: under `--assembly_mode coassembly`, the CheckM2
  and GTDB-Tk reports moved from `04_binning/per_sample/coassembly/quality/`
  and `.../taxonomy/` to `04_binning/coassembly/quality/checkm2/` and
  `04_binning/coassembly/taxonomy/gtdbtk/`, next to the co-assembly bins; a
  co-assembly run no longer creates `04_binning/per_sample/`. Per-sample mode
  is unchanged
- **K-MARVEL database provenance**: a `--custom_kmarvel_db` without a
  `DB_VERSION` file was recorded as `card: "custom"` in
  `software_versions.yml`; it is now recorded as
  `custom (CARD_protein_all.fasta md5 <md5>)`, so the CARD reference used can be
  traced (`bin/kmarvel_db_version.sh`, tested by
  `tests/bin/test_kmarvel_db_version.sh`). A `DB_VERSION` file, as written by
  `FORMAT_KMARVEL_DB`, is still recorded as is (first line)

### Removed

- **Breaking**: the MetaCerberus contig annotation path —
  `--contig_level_metacerberus`, `--metacerberus_hmm`,
  `--metacerberus_minscore`, `--metacerberus_evalue` and its
  `07_functional_annotation/contigs/` output. Superseded by the contig
  functional branch (`--contig_level_functional`: eggNOG-mapper v3 +
  run_dbcan, TPM / copies-per-genome-equivalent tables), which also runs under
  `--assembly_mode coassembly`; one contig annotation backend is kept so the
  results have a single threshold model and ontology mapping. Passing any of
  these params aborts the run at startup

### Fixed

- **ARGS_OAP time limit**: stage one ran with `process_low` (4 CPUs, 8 h) and
  took 6-8 h on 55-70 M read pairs; one sample hit the limit, and on the
  local executor the kill surfaced as `process hasn't exited`, which stops
  the run instead of retrying. ARGS_OAP now uses `process_medium` (8 CPUs,
  36 GB, 24 h); `docs/troubleshooting.md` explains the message
- **Two `Reads_report.csv` written to one path**: with a taxonomic profiler,
  READS_REPORT and TAXONOMY_REPORT both published `Reads_report.csv` to
  `01_quality_control/summary/`, and the one published last was kept, so a
  run could keep the copy without the profiler's DB, Unclassified and
  Classified columns. TAXONOMY_REPORT's copy is now
  `02_taxonomy/tables/Reads_report_taxonomy.csv`; the QC path always holds
  READS_REPORT's, also when the taxonomy step fails. The summary reads the
  taxonomy copy first and falls back to the QC path (older results)
- **Results summary review fixes**:
  - a gene called without an RPKM (no `Reads_report.csv`, or no reference
    length) is blank in `arg_gene_rpkm`, `arg_family_rpkm`,
    `arg_drug_class_rpkm` and `arg_priority_panel`, and 1 in
    `arg_presence`, instead of 0 (not detected)
  - a tier rule outside 1-4 stops with a clear error instead of failing
    while the report is written; the tier table keeps its categories when
    the rule file has no `note` column
  - `n_species` and `shannon` leave out the pooled `unassigned species`
    row (taxa without a species name); Shannon is renormalised over the
    named species
  - the default batch label is the output directory's name, or its
    parent's when it is called `results` (as `bin/results_summary.py`
    labels a directory), reduced to the `--summary_batch_label` pattern so
    it is safe on the command line
  - the `--summary_group_by` / `--summary_order_by` check runs only with
    `--results_summary true`
  - a `--custom_rgi_card_db` is checked at launch for both k-mer files
    RGI_KMER needs (`<k>mer_database.json` and `amr_<k>mer.txt`)
  - tiers are assigned once per distinct ARO annotation, not once per row;
    the status and metadata-merge code is shared; the unused
    `rgi_gene_detected` column is gone
  - the summary no longer fails when KARGA or KARGVA has no passing hit:
    ARG_NORM_REPORT now always writes `CPM_ARGs` and `copies_per_cell`
    (they were left out of a file without hits), and the summary reads a
    file without hits, with or without those columns, as no MEGARes gene
  - `--summary_min_clean_pairs 0` no longer crashes the QC figure (the
    threshold line at 0 on a log axis)
  - when merged batches ran different ARG methods, the report header's
    tier 1 list also covers the samples outside the method comparison
  - a K-MARVEL row raised to tier 4 to match the reference method gets
    its CARD SNP check (`resistance_snp_status`) instead of a blank
  - a 0 setting (e.g. `--rgi-min-coverage 0`) is listed in the README's
    inputs instead of being dropped as unset
  - the RGI placeholder table of a sample without hits is recognised from
    one column or row, not by stripping the whole table
  - the `--results_summary` help and schema text list every result the
    summary covers

- **MEGAHIT out of memory on a large host**: MEGAHIT was run without `-m`, so
  it sized its graph build to 90% of the host's RAM instead of the task's
  memory, and the container killed it at the task limit (a 404M-read
  co-assembly needed ~86 GB at k=21 against the 72 GB first attempt). Its
  exit status, 247, was not in the retry list, so the run stopped instead of
  retrying with more memory. MEGAHIT now runs through `bin/megahit_run.sh`:
  `-m` is 90% of the task memory (a `-m` in `ext.args` still wins), and
  running out of memory exits 137, so the task is retried with twice the
  memory

- **kraken2 classified-reads plot at low classification**: `taxonomy_report.py`
  guessed percent vs fraction from the values, so a kraken2 run where every
  sample was at most 1% classified was drawn as up to 100%; it now goes by
  the profiler (kraken2 reports percent)

- **Report plots readable on large runs**
  - `Box_plot_reads.png`: the y axis had a tick every 5M reads (57
    overlapping labels at 285M); ticks are now chosen automatically
  - Taxonomy bar plots: `--taxonomy_top_n_taxa` taxa are chosen so that each
    sample's most abundant taxa are shown (the global top-N by mean abundance
    put whole samples in "Other"); samples in natural order (`V2` before
    `V10`), figure height and label size follow the sample count
  - Classified-reads plot: samples sorted, legend outside the bars

- **RGI summary counted a zero-hit sample as one gene**: when `rgi bwt` maps too
  few reads, `RGI_BWT` writes placeholder tables with one all-NA row, and
  `RGI_REPORT` counted that row (`Total_AMR_Genes` 1). The row is now dropped:
  the sample is reported with 0 genes and 0 reads, and shown as 0 in
  `RGI_sample_amr_counts.png`. `Total_Mapped_Reads` is written as an integer
  (a mixed run turned it into `17.0`), and samples are sorted

- **ARG blob plot showed every contig without an ARG as "Other"**
  (`ARG_blob_plot.png`): their class is empty, and the top-nine relabel turned
  it into "Other" (82.5% of the contigs in the test data). The plot now shows
  the ARG-carrying contigs only, as the original R script did, so "Other"
  holds only the ARG classes beyond the top nine; the bar panel gives
  percentages of the ARG-carrying contigs

- **`Contig_tax_and_arg_prediction.tsv` is tab-separated**: it was written
  comma-separated despite its name. Scripts that read it as CSV must switch
  to tab-separated

- **RGI pathogen of origin was always empty**: neither database step
  (`RGI_LOAD`, `RGI_LOAD_WILDCARD`) loaded the WildCARD k-mer database, and
  `rgi kmer_query` exits 0 without output when it is missing, so `RGI_KMER`
  wrote `{}` for every sample. Both steps now load it
  (`bin/rgi_load_kmers.sh`, for `--rgi_kmer_size`). A missing k-mer database
  is never fatal, so `rgi bwt` and its reports always run: the load step warns
  when the WildCARD files lack it, `RGI_KMER` warns and writes an empty result
  when the database has none, and the step is skipped with a warning at launch
  for a `--custom_rgi_card_db` without it and for the downloaded CARD-only
  database (`--rgi_include_wildcard false`). Databases prepared by earlier
  runs lack the k-mers: load them as in `docs/RGI_WILDCARD_USAGE.md` (Step 4)
  or prepare again
- **RGI CARD version read as empty on the download path**: `RGI_LOAD` looked for
  `"version"` in `card.json`, which stores it as `"_version"`, so
  `rgi wildcard_annotation` got no version and `software_versions.yml` a blank
  `card:`
- **`RGI_BWT` took hours per sample and kept tens of GB of temporary files**:
  RGI 6.0.3 writes every read KMA reports, unplaced ones included (over 99%
  of a metagenome), to the SAM, two BAM copies and a read table, and then
  scans that table once per reference sequence with reads (81M rows times
  7,859 references for an 81M-read sample). It also removes its temporary
  files only with `--clean`, which the pipeline did not pass: a task kept
  about 390 bytes per clean read (1.2 TB for a 32-sample batch). `RGI_BWT` now
  runs rgi through `bin/rgi_bwt_kma_fix.py`, which runs RGI's own KMA command
  but keeps only the SAM records placed on a reference, and stops when KMA
  fails (RGI ignored it, so a failed alignment looked like a sample without
  hits), and passes `--clean`. The placed reads, and so the gene and allele
  tables, are unchanged; `overall_mapping_stats.txt` and
  `artifacts_mapping_stats.txt` now count the placed reads only. The log no
  longer holds one samtools warning per unplaced KMA read (3.5 GB for an
  81M-read sample). Tests: `tests/bin/test_rgi_bwt_kma_fix.sh`

- **Co-assembly read pairing for sample ids such as `S_1.2`**: MEGAHIT, BOWTIE2_SAMTOOLS
  and BWAMEM3_SAMTOOLS sorted the pooled reads into R1/R2 lists by a substring match on
  `_R1_` or `_1.`, so `S_1.2_R2_clean.fastq.gz` went to R1 and the R1/R2 lists no longer
  matched. The role now comes from the last explicit R1/R2/singleton token of the file
  name, falling back to a bare 1/2 only when there is none, so a chunk number
  (`S1_R1_2.fastq.gz`) does not override it. A file whose name has no role (with
  `--quality_control false`, e.g. `sample.fastq.gz`) was silently left out of the
  co-assembly; it is now an error, and `singletons` counts as a singleton token. All
  three modules use `bin/read_role.sh`, tested by `tests/bin/test_read_role.sh`
- **Bin reports can no longer be skipped silently**: a CheckM2 or GTDB-Tk report whose
  file name matched no sample was dropped with a warning; if that happened for every
  sample, BIN_QUALITY_REPORT, BIN_TAX_REPORT and BIN_SUMMARY never ran. It is now an error
- **Index / aligner mismatch under `--read_aligner bowtie2`** is caught at launch (a
  bwa-mem3 index, e.g. the one the pipeline publishes by default, was failing inside
  BOWTIE2_DECONTAMINATE with an unclear bowtie2 error), and BOWTIE2_DECONTAMINATE itself
  stops with a clear message when its index has no `*.bt2` files. Both index checks now run
  with the other launch checks, before any channel starts: raised later, the error could
  be reported as "Cannot find samplesheet file" instead
- `nextflow_schema.json`: the `deeparg_arg_alignment_evalue` default is the string
  `"1e-10"`, matching its type and `nextflow.config`
- **A run in which no sample yields a bin now finishes**: BIN_TAX_REPORT, BIN_QUALITY_REPORT
  and BIN_SUMMARY warn when every input report is empty and write header-only tables (and
  placeholder plots), instead of failing (audit #13 made them fail; with SemiBin and MetaWRAP
  now finishing on binless samples, that case is a valid outcome). A missing input file
  still fails, and so does BIN_SUMMARY when CheckM2 or GTDB-Tk reports bins the other does
  not. In co-assembly mode the bin depth step is skipped for an empty refined bin set (both
  aligners failed on an empty reference)
- **SemiBin and MetaWRAP stopped the run for a sample without usable bins**:
  - SemiBin2 fails when it finds no bins (its reclustering gets an empty ORF file). SEMIBIN now writes an empty bin directory with a note; any other failure still fails.
  - MetaWRAP stopped on any empty bin set: a binner with no bins, a hybrid set when two binners' bins share too little sequence, or a final set when no bin passes the thresholds. It also counted notes such as `SKIPPED.txt` as bins.
  - METAWRAP now passes only binners with FASTA bins. When an empty hybrid set stops refinement, it refines the binners' own sets (`--skip-refinement`). When no bin reaches `--metawrap_completeness` / `--metawrap_contamination`, it writes an empty refined set with a note.
  - MetaWRAP 1.2 cannot refine a single bin set. When only one binner found bins, METAWRAP keeps the ones that pass MetaWRAP's own CheckM thresholds.
  - MetaWRAP exits 0 when it stops on an error, so the outcome is now read from its output directories.
- **eggNOG-mapper search failed on real-size samples instead of retrying with more
  memory**: diamond peaks at ~80 GB against the eggNOG 7 database, above `process_high`'s
  72 GB, and emapper reports the out-of-memory kill as exit 1, which the memory retry
  (exit 130-145) skips. EGGNOG_MAPPER_SEARCH now exits 137 when diamond was killed or ran
  out of memory, and asks for 96 GB on the first attempt (192 GB on the retry, clamped by
  `--max_memory`). Once `--max_memory` caps the request, an out-of-memory failure is not
  retried: the retry would repeat a multi-hour search with the same memory
- **`software_versions.yml` was invalid YAML for every assembly run**: Nextflow strips a
  script's common indentation, and heredoc bodies at column 0 (BBMAP, CONTIG_FILTER_SUMMARY
  and ten module stubs), `\n` in `printf` strings that Groovy turned into real newlines
  (BOWTIE2_SAMTOOLS*) and a three-tab python heredoc (GTDB_TK_BATCH) left `END_VERSIONS`
  indented. `<<-` strips only tabs, so the terminator was written into `versions.yml`
  and the aggregate failed to parse. Heredoc bodies now sit at the script's base
  indentation. BBMAP's and METAWRAP's versions were also empty, and CALCULATE_DEPTH wrote the
  usage text of `jgi_summarize_bam_contig_depths` (it has no `--version`).
  The stub pipeline tests now parse `software_versions.yml`
- **bowtie2 decontamination failed for sample ids containing a dot** (#21): without a
  `%` in the `--un-conc-gz` path, bowtie2 inserts the mate number before the last
  extension, which for `S1.2` lands inside the id, so the clean reads were not found
  (exit 1). The path now carries `%`. New nf-test `tests/modules/bowtie2_decontaminate.nf.test`
- **Documentation checked against the code**:
  - databases published to `--databases_dir` are not read back by later runs (only through
    the `--custom_*` options or `-resume`), so the index is rebuilt on a fresh run;
  - the launch memory warning covers only the built-in `--host_db`;
  - Bracken and RGI output names, the `kargva/` output, the DeepARG download trigger
    (`arg_bin_clustering` too), the Kraken2/CheckM2 manual-download commands;
  - examples that failed or did nothing: a Bowtie2 index under the default aligner, a
    128 GB limit for a three-genome bwa-mem3 index, an end-to-end `--bowtie_score_min`
    under `--local`, MetaWRAP thresholds with a single binner;
  - the WildCARD files given with `--custom_rgi_wildcard` must be decompressed, and
    `--rgi_include_wildcard` also applies to `rgi bwt`;
  - the workflow diagram's default mode and edges, the cloud retry note, and
    CONTRIBUTING's module and lint rules
- **`--custom_host_fasta` lists**: the docs and examples (including
  `examples/params_multiple_genomes.yaml`) give several host genomes as a
  comma-separated list, but the pipeline read the whole string as one path, which did
  not exist. The list is now split on commas (spaces around them are ignored) and all
  the FASTAs go into the one combined index. New pipeline nf-test with a two-file list
- **bowtie2 decontamination failed with a large (`.bt2l`) index**: BOWTIE2_DECONTAMINATE
  read the basename of `contaminants.1.bt2l` as `contaminantsl`, although the launch check
  accepts `.bt2l`. bowtie2-build writes a large index on its own for references over
  ~4 Gbp, such as several host genomes. BOWTIE2_BUILD_COMBINED now also passes `ext.args`
  (the new test builds with `--large-index`)
- **`low_disk` was described as freeing disk during the run**: the `--help` text, the
  profile description, the launch warning, the README, the manual and
  `docs/DISK_OPTIMIZATION.md` called its work-dir cleanup "progressive" and promised a
  60-70 % lower peak. Nextflow's `cleanup = true` deletes the work directory only once the
  run succeeds, and the in-task cleanup of temporary files runs in every profile, so peak
  disk use is the same with or without `low_disk`. Behaviour is unchanged; the text now
  says what happens: work dir deleted after a successful run (that run cannot be
  resumed; an interrupted one can), task directories reused from the cache or left by
  the interrupted run survive a resumed run's cleanup (delete `work/` once resuming is no
  longer needed), and running out of disk mid-run calls for a
  larger `-work-dir` filesystem or fewer concurrent tasks
- **`max_*` set in a `-c` config file were ignored**: `process.resourceLimits` copied
  `max_cpus`/`max_memory`/`max_time` while `nextflow.config` was parsed, so values from a
  later `-c` file (the documented way in manual §10.1, and `-c conf/my_institution.config`)
  never reached it. After `-profile test`, tasks stayed at 2 CPUs / 6 GB / 6 h. The limits are
  now read at task time, so `-c`, `-params-file` and the command line all work. The local
  executor's pool still follows only the command line and `-params-file` (executor settings
  are fixed at parse time): a `-c` file that raises the limits must also set
  `executor { cpus; memory }`. Otherwise the run fails with "Process requirement exceeds
  available CPUs" and stops: a local task that never gets an exit status because it can
  never be scheduled now ends the run (`terminate`) instead of the default `finish`, which waited forever for
  queued tasks that could not run either (the same applies to the download labels, `RGI_BWT` and
  `MICROBECENSUS`). Documented in the manual, parameters and troubleshooting; regression test
  `tests/bin/test_resource_limits.sh`
- **SUPER-FOCUS, HUMAnN and every database download under singularity/apptainer**: the
  separately built singularity images of the SUPER-FOCUS and HUMAnN containers have no
  `gzip`, so `SUPERFOCUS`/`HUMANN` failed at once (exit 127) under apptainer and singularity.
  Under `-profile singularity`, the singularity image of the shared download container has
  no `tar`, so the kraken, GTDB-Tk, CheckM2, BLAST nt, taxdump and HUMAnN downloads would
  have failed. These processes now use the docker image for every engine (singularity and
  apptainer convert it at pull time). Found by the end-to-end apptainer runs; CI is
  docker-only
- **`BLOBTOOLS` failed under singularity/apptainer**: `blobtools create` wrote its parsed
  taxonomy into its own install directory, which is read-only in a singularity image. It
  now writes `nodesDB.txt` in the task directory (`--db`); regression test
  `tests/bin/test_blobtools_readonly.sh`
- **`BLOBPLOT` / `ARG_BLOBPLOT` always failed**: saving the figure object (`.pkl`) failed on
  its lambda tick formatters after the PNG had been written, so `--contig_tax_and_arg` runs
  stopped there. The formatters are now named functions. PNGs are unchanged
- **Steps that combine many samples re-ran on `-resume`**: the co-assembly (`MEGAHIT` and
  the co-assembly read alignment), `METAWRAP`, `CHECKM2_BATCH`, `GTDB_TK_BATCH`,
  `BIN_SUMMARY` and the report steps (reads, contig filtering, taxonomy, ARG, RGI, BlobTools,
  ARG clustering, bin quality/taxonomy) received their inputs in task-completion order, so
  their cache key could change between otherwise identical runs (observed: 2 of 8 resumes of
  one configuration re-ran the CheckM2 or GTDB-Tk batch). Inputs are now ordered by sample
  (pooled co-assembly reads keep each sample's R1/R2 together). An existing run may re-run
  each of these steps once on its next `-resume`. MEGAHIT multi-threaded output is not
  byte-reproducible regardless of input order, so a re-run co-assembly can differ slightly
- **`--<option> false` turned options on with Nextflow 26.04**: Nextflow 26.04 (verified
  26.04.4 and 26.04.6) passes every command-line parameter as text, so `--include_binning
  false` was the truthy string `"false"`: the branch ran, or a launch check rejected a valid
  run (e.g. `--assembly_mode none --include_binning false`), and a numeric
  `--min_read_sample` crashed QC with `Cannot compare java.lang.Integer ... with
  java.lang.String`. Every on/off option is now read through one helper
  (`subworkflows/local/utils_params.nf`), and `min_read_sample` is converted before the
  comparison. Nextflow 25.10 converted these values itself and was not affected. Also fixed
  on every version: `--azure_delete_pools false` was turned back on by a `?: true` default
- **Co-assembly bin coverage counted one contig too many**: `BEDTOOLS` treated the
  whole-genome summary rows of `genomeCoverageBed` (chrom `genome`) as one more contig, so
  `Total_contigs` in `04_binning/coassembly/coverage/*_bin_depth.tsv` and
  `04_binning/coassembly/summary/Bin_summary.csv` was one too high for every bin, and the
  genome-wide mean depth was averaged in with the per-contig means in `average_cov`. On a
  real co-assembly the averages moved by up to 1.6 % (e.g. 39.04 → 39.67). Regression case
  in `tests/modules/bedtools.nf.test`
- **`RGI_REPORT` never found the k-mer results**: it is given `RGI_KMER`'s
  `*_<k>mer_analysis.json` files but looked for `*mer_analysis.txt`, so its log always
  said "Found 0 k-mer files". It now counts the JSON files. The k-mer results are still
  not part of the summary or plots (they never were); the report output is unchanged.
  New test `tests/modules/rgi_report.nf.test`
- **`RGI_KMER` failed with `--rgi_aligner bowtie2` or `bwa`**: the module required a
  `*_<k>mer_analysis.txt` output that `rgi kmer_query --bwt` (RGI 6.0.3) never writes. It
  writes `_<k>mer_analysis.json`, `.allele.txt` and `.gene.txt`. The required output is
  removed (nothing used it), and the manual's output tree lists the real files. Stub runs
  never showed this, because the stub created the missing file. (With the default `kma`
  the task hung before reaching this check, a separate bug)
- **`RGI_KMER` hung forever with the default `--rgi_aligner kma`, or when a worker ran
  out of memory**: RGI 6.0.3's `kmer_query` read KMA BAMs from the wrong SAM columns (the
  reference name as the flag, base qualities as the sequence), so every worker stopped
  with "failed to parse BAM file". Separately, a worker killed for exceeding the task
  memory (the default 8 workers needed more than 32 GB on the test data) never reported
  back. In both cases the parent process waited forever: the task sat idle until its 24 h
  limit, then was retried at 48 h and 72 h. `RGI_KMER` now runs `kmer_query` through
  `bin/rgi_kmer_kma_fix.py`, which reads every BAM with the standard columns, fails as soon
  as a worker dies (exit 137 for an out-of-memory kill, so Nextflow retries with more
  memory), and stops a run that takes longer than 2 h (`ext.time_limit`, exit 124). The
  patch is version-guarded to RGI 6.0.3. `RGI_KMER` now uses 2 CPUs (2 workers): on the
  test data it used 21 GB, within its 36 GB. Tests: `tests/bin/test_rgi_kmer_kma_fix.sh`
- **`RGI_BWT` / `RGI_KMER` could not run `bin/` scripts**: the RGI image's
  `/usr/local/env-execute`, which both modules source, resets `PATH`, dropping the
  pipeline's `bin/`. `RGI_BWT`'s low-coverage fallback (`rgi_bwt_lowcov_fallback.sh`)
  would have failed with "command not found" exactly when it was needed. The pipeline's
  `PATH` entries are now kept
- Test harness: the blast bad-md5 fixture could leave the checksum unchanged
  (1-in-16 CI flake); it now always corrupts it

## [1.1] - 2026-08-25

### Added

- **Provenance tracking (audit #27)**
  - `versions.yml` emitted by every live module (previously ~20 in-use modules
    emitted none) and, for the first time, aggregated: each run now writes a
    deduplicated `pipeline_info/software_versions.yml` covering every executed
    process plus the pipeline and Nextflow versions
  - DeepARG database downloads record the tool version and download date; the
    CARD version captured by RGI_LOAD now actually reaches the aggregate report

- **RGI AMR Prediction**
  - New `--rgi_prediction` parameter to enable AMR gene prediction with pathogen-of-origin analysis
  - RGI_LOAD module for automatic CARD database download and preparation
  - RGI_LOAD_WILDCARD module for combining separate CARD and WildCARD databases
  - RGI_BWT module for read-level AMR gene alignment using KMA aligner
  - RGI_KMER module for k-mer based pathogen-of-origin prediction
  - RGI_REPORT module for aggregated results and visualizations
  - Support for WildCARD variants for extended allelic diversity
  - Custom database support via `--custom_rgi_card_db` and `--custom_rgi_wildcard` parameters
  - Flexible database options: automatic download, pre-prepared complete database, or separate CARD + WildCARD
  - Comprehensive documentation in `docs/RGI_WILDCARD_USAGE.md` (originally also `docs/RGI_IMPLEMENTATION_{PLAN,SUMMARY}.md`, since consolidated)
  - RGI-specific parameters: `rgi_card_version`, `rgi_include_wildcard`, `rgi_aligner`, `rgi_kmer_size`, `rgi_min_kmer_coverage`
  - Integration with PREPARE_DATABASES subworkflow for automatic database management
  - Output includes per-sample results, pathogen predictions, and multi-sample summary reports with plots

### Changed

- Updated README.md with RGI feature description and usage examples
- Updated docs/manual.md with RGI parameters, output structure, and usage examples
- Updated docs/parameters.md with complete RGI parameter reference
- Enhanced ARG prediction capabilities with complementary tool (RGI alongside KARGA/KARGVA)
- Documentation corrected to match real behavior (audit #26): sample read-count
  filter default, database storage location, output trees, profile lists;
  MultiQC removed from docs, dead wiring, and vendored modules (audit #28;
  re-enabling requires re-vendoring via `nf-core modules install multiqc`)
- Minimum Nextflow version raised from 23.04.0 to **24.04.0** with parameter
  validation now enforced by the nf-schema 2.4.2 plugin (audit #25): unknown
  `--params` are a hard startup error
- Dead documented knobs fixed (audit #14): `fastp_qualified_quality_phred`
  wired, METABAT2 selectors collapsed (pTNF/minCV/minCVSum now delivered),
  `bbmap_lenght` doc typo corrected, unused `mmseqs_*` params removed

### Removed

**Breaking**: parameter validation is now strict (nf-schema
`failUnrecognisedParams`), so passing any removed parameter aborts the run at
startup instead of being silently ignored. Update existing command lines and
`-params-file` YAMLs accordingly.

- `--custom_phiX_index` and `--custom_bowtie_host_index` — replaced by
  `--custom_decontamination_index` (pre-built combined index),
  `--custom_phiX_fasta`, and `--custom_host_fasta` (the pipeline builds the
  combined Bowtie2 index from FASTAs); the old params had been silent no-ops
- `--enable_work_cleanup` — was never consumed; work-dir cleanup is Nextflow's
  `cleanup = true`, enabled by `-profile low_disk`
- `--kraken_db_used`, `--sourmash_db_name` — report database names are derived
  from the selected database
- `--store_filtered_contigs`, `--store_refined_bins` — were no-ops; filtered
  contigs and refined bins are always published
- `--tracedir` — trace/report/timeline/DAG always go to
  `<output>/pipeline_info/`
- `--validationShowHiddenParams`, `--validationSchemaIgnoreParams` —
  nf-validation 1.x options superseded by the nf-schema `validation {}` scope
- All eleven `--mmseqs_*` parameters — unused; clustering settings are fixed
  tiers in `modules/local/clustering`
- Renamed: `--bbmap_lenght` → `--bbmap_length`

### Fixed

- Extensive audit-fix series on branch `fix-pending-issues` (2026-08): host
  decontamination DB and `--local` scoring, database download containers and
  script hardening, kraken2/bracken report parsers, contig tax/ARG arm wiring,
  singleton read handling through QC, storeDir misuse, report failure masking,
  silent no-op feature combos, edge-data crashes, publishing gaps (the
  itemized 30-point list and fixing commits are recorded in the per-item
  `fix:`/`feat:` commit messages on that branch)

## [1.0.0] - 2024-01-10

### Added

- **Pipeline Infrastructure**
  - Manifest block with pipeline metadata and Nextflow version requirement (>=23.04.0)
  - Parameter validation schema (`nextflow_schema.json`)
  - Samplesheet validation schema (`assets/schema_input.json`)
  - Input validation subworkflow with file existence checks
  - Dynamic resource allocation with `check_max()` function
  - Execution reports: timeline, trace, and DAG

- **Modular Architecture**
  - Subworkflows for logical pipeline components:
    - `INPUT_CHECK` - Samplesheet validation
    - `PREPARE_DATABASES` - Database download and formatting
    - `QC` - Quality control and host filtering
    - `TAXONOMY` - Taxonomic profiling
    - `ASSEMBLY` - Genome assembly
    - `BINNING` - Metagenomic binning and refinement

- **Cloud Support**
  - AWS Batch configuration with S3 and Fusion filesystem
  - Google Cloud Batch/Life Sciences configuration
  - Azure Batch configuration with auto-scaling pools
  - SLURM HPC configuration with Singularity

- **Container Profiles**
  - Docker profile
  - Singularity profile
  - Podman profile
  - Apptainer profile

- **Module Improvements**
  - Fixed bash null checks using Groovy conditionals
  - Added `versions.yml` output to some modules (full coverage and aggregation
    landed later, in [1.1])
  - Added `stub` blocks for some modules (full coverage landed later)
  - Added `meta.yml` descriptors for key modules

- **Documentation**
  - Comprehensive deployment guide (`docs/deployment.md`)
  - Updated README with profiles and cloud examples
  - Contributing guidelines (`CONTRIBUTING.md`)
  - Example samplesheet

### Changed

- Refactored main.nf with improved help/version handling
- Consolidated params block in nextflow.config
- Updated process labels with memory and time definitions
- Improved error handling and retry strategies

### Fixed

- Bash null string comparison issues in FASTP, MEGAHIT, BOWTIE2 modules
- Resource allocation respects max limits on retries

## [0.1.0] - Initial Release

### Added

- Initial pipeline implementation
- Basic QC, assembly, binning workflow
- Kraken2 and Sourmash taxonomic profiling
- ARG prediction with DeepARG, KARGA, KARGVA
- Binning with MetaBAT2, SemiBin, COMEBin
- Bin refinement with MetaWRAP
- Quality assessment with CheckM2
- Taxonomic classification with GTDB-TK

[Unreleased]: https://github.com/gene2dis/BugBuster/compare/v1.1...HEAD
[1.1]: https://github.com/gene2dis/BugBuster/releases/tag/v1.1
[1.0.0]: https://github.com/gene2dis/BugBuster/releases/tag/v1.0.0
[0.1.0]: https://github.com/gene2dis/BugBuster/tree/v0.1.0
