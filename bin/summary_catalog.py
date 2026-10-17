"""What a BugBuster run publishes, stage by stage: the catalog behind the README and LLM_GUIDE.md of
bin/results_summary.py. Paths are relative to the run's --output directory.

Each stage lists the params that turn it on, the processes whose software_versions.yml entries show it ran,
what 00_summary/ takes from it, and its key published files: (path pattern, content, how the sample is
identified, units). Keep it in step with config/modules.config (publishDir) when a module changes.
"""

STAGES = [
    {
        "key": "qc", "title": "Read QC and host/contaminant removal",
        "params": "--quality_control, --read_aligner (bwamem3 | bowtie2)",
        "processes": ["FASTP", "QFILTER", "BWAMEM3_DECONTAMINATE", "BOWTIE2_DECONTAMINATE", "READS_REPORT"],
        "summary": "samples_qc (read pairs per step, depth class, low-depth flag), figure qc_read_fate",
        "files": [
            ("01_quality_control/summary/Reads_report.csv", "reads per step (raw, fastp, clean); with a taxonomic "
             "profiler also its DB, Unclassified and Classified columns", "column Id",
             "single reads (both mates); Classified: kraken2 0-100, sylph/sourmash 0-1"),
            ("01_quality_control/fastp/<sample>/<sample>.fastp.json", "fastp report (read lengths, quality)",
             "directory and file name", "-"),
            ("clean_reads/<sample>/", "clean reads (only with --store_clean_reads)", "directory", "FASTQ"),
        ],
    },
    {
        "key": "taxonomy", "title": "Taxonomic profiling",
        "params": "--taxonomic_profiler (sylph | kraken2 | sourmash | none)",
        "processes": ["SYLPH_PROFILE", "SYLPH_TAXPROF", "KRAKEN2_KRAKEN2", "BRACKEN_BRACKEN", "SOURMASH",
                      "TAXONOMY_PHYLOSEQ"],
        "summary": "taxonomy_long, taxonomy_<rank>_wide, taxonomy_species_prevalence, watchlist; taxonomy heatmaps",
        "files": [
            ("02_taxonomy/tables/<profiler>_<db>_otu_table.tsv", "abundance per taxon (rows) and sample (columns)",
             "column names", "sylph: relative abundance % of the classified community; kraken2: Bracken reads; "
             "sourmash: abundance-weighted hashes"),
            ("02_taxonomy/tables/<profiler>_<db>_tax_table.tsv", "lineage of each taxon id (Kingdom ... Species)",
             "-", "-"),
            ("02_taxonomy/tables/sylph_<db>_read_counts.tsv", "sylph: estimated read pairs per genome",
             "column names", "read pairs"),
            ("02_taxonomy/sylph/<sample>_sylph_<db>.tsv", "raw sylph profile: ANI, coverage per genome",
             "file name", "%"),
            ("02_taxonomy/kraken2/, 02_taxonomy/bracken/", "kraken2 reports and Bracken estimates", "file name",
             "reads"),
            ("02_taxonomy/sourmash/*.with-lineages.csv", "sourmash gather results", "file name", "hashes"),
            ("02_taxonomy/phyloseq/", "phyloseq objects (with --create_phyloseq_rds)", "-", "-"),
        ],
    },
    {
        "key": "kmarvel", "title": "Read-level ARGs: K-MARVEL (CARD)", "params": "--kmarvel_prediction",
        "processes": ["KMARVEL"],
        "summary": "arg_long (method kmarvel), arg_gene_rpkm and the other ARG matrices when K-MARVEL is the primary "
                   "method, arg_method_comparison",
        "files": [
            ("05_arg_prediction/read_level/kmarvel/<sample>/<sample>_kmarvel.tsv", "ARG calls with CARD annotation, "
             "supporting reads, coverage, identity, mutations", "file name", "supporting read pairs"),
            ("05_arg_prediction/read_level/kmarvel_summary/KMARVEL_*.tsv", "run-level long table and matrices",
             "column sample / column names", "read pairs; drug-class matrix = number of ARGs"),
        ],
    },
    {
        "key": "rgi", "title": "Read-level ARGs: RGI bwt + k-mer pathogen of origin (CARD + WildCARD)",
        "params": "--rgi_prediction (--rgi_include_wildcard, --rgi_summary_min_coverage)",
        "processes": ["RGI_BWT", "RGI_KMER", "RGI_REPORT"],
        "summary": "arg_long (method rgi, genes at >= --rgi_summary_min_coverage), arg_method_comparison, "
                   "rgi_kmer_hosts / rgi_host_check",
        "files": [
            ("05_arg_prediction/read_level/rgi/<sample>/<sample>_rgi_bwt.allele_mapping_data.txt", "every reference "
             "sequence with at least one mapped read (no detection call)", "file name",
             "mapped reads (mates counted separately), % coverage"),
            ("05_arg_prediction/read_level/rgi/<sample>/<sample>_rgi_bwt.gene_mapping_data.txt", "the same per ARO",
             "file name", "mapped reads"),
            ("05_arg_prediction/read_level/rgi_kmer/<sample>/<sample>_rgi_kmer_<k>mer_analysis.gene.txt",
             "k-mer pathogen of origin per ARO", "file name", "reads"),
            ("05_arg_prediction/read_level/rgi_summary/RGI_summary_report.csv", "per sample: genes and reads, all "
             "hits and at the coverage cutoff (*_cov<N>)", "column Sample", "genes, mapped reads"),
        ],
    },
    {
        "key": "karga", "title": "Read-level ARGs: KARGA / KARGVA (MEGARes) normalised by ARGs-OAP",
        "params": "--read_arg_prediction",
        "processes": ["KARGA", "KARGVA", "ARGS_OAP", "ARG_NORM_REPORT"],
        "summary": "arg_karga_long, arg_karga_group_cpc, arg_karga_class_cpc, arg_kargva_long; samples_qc "
                   "args_oap_n16S / args_oap_ncell / karga_*; figure arg_karga_class_heatmap",
        "files": [
            ("05_arg_prediction/read_level/summary/KARGA_norm.csv", "KARGA genes with >= 90% coverage, with the "
             "sample's ARGs-OAP counts", "column sample",
             "copies_per_cell = k-mer depth / ARGs-OAP nCell; CPM_ARGs = k-mer depth per million single reads"),
            ("05_arg_prediction/read_level/summary/KARGVA_norm.csv", "KARGVA variants with >= 2 SNP k-mers and >= "
             "80% coverage; a sample with no hit has one row with blank ARG fields", "column sample",
             "as KARGA_norm.csv"),
            ("05_arg_prediction/read_level/karga/<sample>_all_reads_KARGA_mappedGenes.csv", "every MEGARes gene with "
             "a k-mer hit (unfiltered)", "file name", "% covered (text with %), k-mer depth"),
            ("05_arg_prediction/read_level/kargva/<sample>_all_reads_KARGVA_mappedGenes.csv", "every KARGVA variant "
             "with a k-mer hit (unfiltered)", "file name", "SNP k-mer hits as hits/total"),
            ("05_arg_prediction/read_level/args_oap/<sample>_args_oap_s1_out/metadata.txt", "ARGs-OAP stage one "
             "counts: nRead, n16S, nCell (stage two, SARG abundance, is not run)", "column sample (<sample>_tmp_reads)",
             "single reads; 16S copies and cells (fractional estimates)"),
        ],
    },
    {
        "key": "read_functional", "title": "Read-level functional profile",
        "params": "--read_level_functional (woltka | superfocus | humann | none)",
        "processes": ["WOLTKA_CLASSIFY", "SUPERFOCUS", "HUMANN", "AGGREGATE_READ_FUNCTIONS"],
        "summary": "function_read_long, function_read_<ontology>; samples_qc func_*; figure function_cog_heatmap",
        "files": [
            ("07_functional_annotation/summary/read_function_abundance.tsv", "per sample and term: native "
             "abundance and CPGE", "column sample_id",
             "abundance_native: Woltka/SUPER-FOCUS reads (mates separately, fractional), HUMAnN RPK; abundance_cpge: "
             "RPK per genome equivalent"),
            ("07_functional_annotation/summary/read_function_wide_<ontology>_{native,cpge}.tsv", "wide matrices per "
             "ontology", "column names", "as above"),
            ("07_functional_annotation/summary/read_annotated_fraction.tsv", "share of ORF-assigned reads with a "
             "term, per ontology", "column sample_id", "reads, fraction"),
            ("07_functional_annotation/summary/read_sample_summary.tsv", "reads assigned, MicrobeCensus AGS/GE, CPGE "
             "status, tool and DB versions", "column sample_id", "reads, bp"),
            ("07_functional_annotation/reads/<backend>/<sample>/", "per-sample backend outputs (Woltka: ORF table "
             "with column <sample>.wol, functions, summary, Bowtie2 log in read pairs)", "directory", "reads"),
        ],
    },
    {
        "key": "microbecensus", "title": "Average genome size (MicrobeCensus)",
        "params": "--microbecensus (runs with a contig or read functional branch)",
        "processes": ["MICROBECENSUS"],
        "summary": "samples_qc mc_avg_genome_size_bp, mc_genome_equivalents",
        "files": [
            ("07_functional_annotation/microbecensus/<sample>/<sample>.ags.tsv", "average genome size and genome "
             "equivalents (from 2M sampled reads)", "column sample_id", "bp; genome equivalents"),
        ],
    },
    {
        "key": "assembly", "title": "Assembly (MEGAHIT)",
        "params": "--assembly_mode (assembly | coassembly | none)",
        "processes": ["MEGAHIT", "BBMAP", "CONTIG_FILTER_SUMMARY"],
        "summary": "not summarised yet",
        "files": [
            ("03_assembly/per_sample/<sample>/<sample>_contigs.fa, 03_assembly/coassembly/coassembly_contigs.fa",
             "raw contigs", "directory", "-"),
            ("03_assembly/.../*_filtered_contigs.fa, *_contig.stats", "length-filtered contigs and their stats "
             "(bbtools stats.sh)", "file name", "contigs, bp"),
            ("pipeline_info/contig_filtering_summary.txt", "contigs before and after the length filter",
             "column Sample_ID", "contigs"),
        ],
    },
    {
        "key": "binning", "title": "Binning and MAGs (SemiBin2 / MetaBAT2 / COMEBin, MetaWRAP, CheckM2, GTDB-Tk)",
        "params": "--include_binning, --binners (MetaWRAP refinement with >= 2 binners)",
        "processes": ["SEMIBIN", "METABAT2", "COMEBIN", "METAWRAP", "CHECKM2_BATCH", "GTDB_TK_BATCH"],
        "summary": "not summarised yet",
        "files": [
            ("04_binning/quality/summary/{Mag_quality_summary,All_Mag_quality_table,Refined_Mag_quality_table}.csv",
             "per-sample assembly: bin quality (CheckM2)", "column Sample", "completeness / contamination %"),
            ("04_binning/taxonomy/summary/MAGs_tax_summary.csv", "per-sample assembly: GTDB-Tk classification",
             "column sample", "-"),
            ("04_binning/coassembly/summary/Bin_summary.csv", "coassembly: quality, GTDB ranks and mean depth per "
             "sample (<sample>_avgcov columns)", "column names", "%, mean depth"),
            ("04_binning/.../quality/checkm2/, .../taxonomy/gtdbtk/", "raw CheckM2 and GTDB-Tk tables", "file name",
             "%"),
        ],
    },
    {
        "key": "contig_arg", "title": "Contig taxonomy and ARGs (BLAST nt, BlobTools, DeepARG)",
        "params": "--contig_tax_and_arg (needs an assembly)",
        "processes": ["NT_BLASTN", "BLOBTOOLS", "DEEPARG_CONTIGS", "ARG_CONTIG_LEVEL_REPORT"],
        "summary": "not summarised yet",
        "files": [
            ("05_arg_prediction/contig_level/summary/Contig_tax_and_arg_prediction.tsv", "one row per contig x ARG "
             "with lineage, coverage and DeepARG class/probability", "column Sample", "raw coverage"),
            ("05_arg_prediction/contig_level/deeparg/<sample>/", "DeepARG outputs", "directory", "-"),
        ],
    },
    {
        "key": "bin_arg", "title": "Bin-level ARG clustering (DeepARG on bins, MMseqs2)",
        "params": "--arg_bin_clustering (needs binning)",
        "processes": ["DEEPARG_BINS", "CLUSTERING"],
        "summary": "not summarised yet",
        "files": [
            ("05_arg_prediction/bin_level/clustering/Arg_cluster_{90,95,99,100}_cluster.tsv", "representative -> "
             "member ARG clusters (no header)", "member names", "-"),
        ],
    },
    {
        "key": "contig_functional", "title": "Contig-level functional (Pyrodigal, eggNOG-mapper, dbCAN, featureCounts)",
        "params": "--contig_level_functional (needs an assembly), --functional_cazy",
        "processes": ["PYRODIGAL", "EGGNOG_MAPPER_ANNOTATE", "RUN_DBCAN", "FEATURECOUNTS_GENES", "AGGREGATE_FUNCTIONS"],
        "summary": "not summarised yet",
        "files": [
            ("07_functional_annotation/summary/function_abundance.tsv", "per sample and term: TPM, CPGE, native "
             "counts", "column sample_id", "TPM of assembled genes; CPGE"),
            ("07_functional_annotation/summary/function_wide_<ontology>_{tpm,cpge}.tsv", "wide matrices", "column names",
             "as above"),
            ("07_functional_annotation/summary/{annotated_fraction,ags_and_ge,gene_abundance,gene_annotations}.tsv",
             "annotated fraction, AGS/GE, per-gene tables", "column sample_id", "-"),
        ],
    },
    {
        "key": "mag_functional", "title": "MAG annotation (Bakta)",
        "params": "--mag_level_functional (needs binning with >= 2 binners)",
        "processes": ["BAKTA_BAKTA"],
        "summary": "not summarised yet",
        "files": [
            ("07_functional_annotation/mags/<sample>/<sample>_<bin>.{tsv,gff3,faa,...}", "Bakta annotation per refined "
             "bin", "file name", "-"),
        ],
    },
]


def stage_status(processes_ran, tables_found):
    """[(stage, status)]: 'ran', 'not run', or 'ran, no table found' (a stage whose processes ran but whose
    summarised table is missing from the results directory)."""
    out = []
    for st in STAGES:
        ran = any(p in processes_ran for p in st["processes"])
        if st["key"] in tables_found:
            status = "ran"
        elif ran and st["summary"] != "not summarised yet":
            status = "ran, no table found"
        elif ran:
            status = "ran"
        else:
            status = "not run"
        out.append((st, status))
    return out
