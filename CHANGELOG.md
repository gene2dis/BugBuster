# Changelog

All notable changes to BugBuster will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

The pipeline self-reports this state as `2.0.0dev` (manifest version) until the
next release is tagged.

### Added

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
  - Covers the pipeline's own alignment steps: host/phiX decontamination (`BWAMEM3_BUILD_COMBINED`, `BWAMEM3_DECONTAMINATE`), reads-to-contigs mapping for binning, contig coverage and gene counting (`BWAMEM3_INDEX`, `BWAMEM3_SAMTOOLS`, `BWAMEM3_SAMTOOLS_PER_SAMPLE`), and co-assembly bin depth (`BINS_REFERENCE`, `BWAMEM3_BINS_INDEX`, `BWAMEM3_BIN_DEPTH`, `BEDTOOLS_BIN_DEPTH`). Woltka, HUMAnN/MetaPhlAn and RGI keep bowtie2 with their own databases
  - bwa-mem3 0.14.0 and samtools 1.22.1 run from one pinned Seqera Containers image. Output is independent of the thread count (`-K` fixed)
  - Decontamination is strict: a pair is clean only when both mates are unmapped. bowtie2 `--un-conc` kept pairs with a single host mate. Clean-read counts are therefore not directly comparable with bowtie2 runs, and the FASTQ comment field is not carried over. Output file names are set explicitly, so sample ids containing dots no longer break decontamination. Mates are aligned without pairing (`-P`) so that `-T` holds for each mate: in paired mode bwa mem skips it for a chosen pair, and random ~20 bp host hits on overlapping mates removed 17.5% of a bacterial mock community (0.1% with `-P`). New `--decontam_min_score` (bwa `-T`, default 30) and `--bwamem3_dense_sa` (denser suffix array: faster mapping for ~12 GB more RAM)
  - The decontamination index is published to `<databases_dir>/bwamem3_index/`. A human host needs ~75 GB RAM to build it (minutes rather than bowtie2-build's hours), and the build retries on bwa-mem3's out-of-budget exit status 3. `--custom_decontamination_index` accepts a bwa-mem3 or bwa-mem2 index under `bwamem3` and rejects a Bowtie2 one at launch
  - Contig alignments drop secondary and supplementary records, keeping one record per read as with bowtie2. Short or divergent local hits become unmapped reads: an alignment needs >= 95% identity (`--bwamem3_min_identity`) and >= 50 aligned bases (`--bwamem3_min_aligned_length`). Without the filter, reads from genomes missing from the assembly inflated gene counts (Spearman vs bowtie2 0.74 in a poorly covered sample) and the depth of near-absent bins (~1.9x). Each library of a pooled co-assembly is mapped separately and the results merged. The co-assembly index is built once instead of once per alignment task
  - Co-assembly bin depth maps each sample once against all refined bins concatenated, instead of once per bin (bwa mem aligns locally, so per-bin mapping would let other bins' reads inflate coverage). It writes the same `*_bin_depth.tsv` table
  - The bioconda bwa-mem3 build needs an AVX2-capable CPU (x86-64-v3)

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

- **Breaking (default aligner)**: `--read_aligner` defaults to `bwamem3` instead of
  `bowtie2`, after a benchmark on five real metagenomes and a spiked mock community
  (`docs/benchmarks/bwamem3_vs_bowtie2.md`): host removal 4.1x faster, reads-to-contigs
  mapping 3.9-6.2x faster with the alignment filter, co-assembly bin depth ~34x faster,
  no human pairs leaked from the spike-in (bowtie2 0.157%) and 0.11% of the mock
  removed (bowtie2 3.0%), the same bins and GTDB species. `--read_aligner bowtie2`
  stays as the fallback, unchanged. What changes for an existing setup:
  - building the human decontamination index needs ~75 GB RAM (bowtie2: ~6 GB) and
    mapping ~14 GB per task (~4 GB). The pipeline warns at launch when `--max_memory`
    is lower, and the build no longer repeats an exit-3 retry that `--max_memory`
    caps at the same memory
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

- **Co-assembly read pairing for sample ids such as `S_1.2`**: MEGAHIT, BOWTIE2_SAMTOOLS
  and BWAMEM3_SAMTOOLS sorted the pooled reads into R1/R2 lists by a substring match on
  `_R1_` or `_1.`, so `S_1.2_R2_clean.fastq.gz` went to R1 and the R1/R2 lists no longer
  matched. The role now comes from the last R1/R2/singleton token of
  the file name (`bin/read_role.sh`, and the same rule in Groovy), tested by
  `tests/bin/test_read_role.sh`
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
  `--max_memory`)
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
