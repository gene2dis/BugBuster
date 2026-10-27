"""The pipeline's stages, for the README and LLM_GUIDE.md of bin/results_summary.py.

Each stage lists the params that turn it on, the processes whose software_versions.yml entries show it ran,
and what the summary takes from it. Keep it in step with main.nf when a stage or process is added.
"""

STAGES = [
    {
        "key": "qc", "title": "Read QC and host/contaminant removal",
        "params": "--quality_control, --read_aligner (bwamem3 | bowtie2)",
        "processes": ["FASTP", "QFILTER", "BWAMEM3_DECONTAMINATE", "BOWTIE2_DECONTAMINATE", "READS_REPORT"],
        "summary": "samples_qc (read pairs per step, depth class, low-depth flag), figure qc_read_fate",
    },
    {
        "key": "taxonomy", "title": "Taxonomic profiling",
        "params": "--taxonomic_profiler (sylph | kraken2 | sourmash | none)",
        "processes": ["SYLPH_PROFILE", "SYLPH_TAXPROF", "KRAKEN2", "BRACKEN", "SOURMASH",
                      "TAXONOMY_PHYLOSEQ"],
        "summary": "taxonomy_long, taxonomy_<rank>_wide, taxonomy_species_prevalence, watchlist; taxonomy heatmaps",
    },
    {
        "key": "kmarvel", "title": "Read-level ARGs: K-MARVEL (CARD)", "params": "--kmarvel_prediction",
        "processes": ["KMARVEL"],
        "summary": "arg_long (method kmarvel), arg_gene_rpkm and the other ARG matrices when K-MARVEL is the primary "
                   "method, arg_method_comparison",
    },
    {
        "key": "rgi", "title": "Read-level ARGs: RGI bwt + k-mer pathogen of origin (CARD + WildCARD)",
        "params": "--rgi_prediction (--rgi_include_wildcard, --rgi_summary_min_coverage)",
        "processes": ["RGI_BWT", "RGI_KMER", "RGI_REPORT"],
        "summary": "arg_long (method rgi, genes at >= --rgi_summary_min_coverage), arg_method_comparison, "
                   "rgi_kmer_hosts / rgi_host_check",
    },
    {
        "key": "args_oap", "title": "Read-level ARGs: ARGs-OAP (SARG)",
        "params": "--read_arg_prediction",
        "processes": ["ARGS_OAP"],
        "summary": "not summarised yet",
    },
    {
        "key": "read_functional", "title": "Read-level functional profile",
        "params": "--read_level_functional (woltka | superfocus | humann | none)",
        "processes": ["WOLTKA_CLASSIFY", "SUPERFOCUS", "HUMANN", "AGGREGATE_READ_FUNCTIONS"],
        "summary": "function_read_long, function_read_<ontology>; samples_qc func_*; figure function_cog_heatmap",
    },
    {
        "key": "microbecensus", "title": "Average genome size (MicrobeCensus)",
        "params": "--microbecensus (runs with a contig or read functional branch)",
        "processes": ["MICROBECENSUS"],
        "summary": "samples_qc mc_avg_genome_size_bp, mc_genome_equivalents",
    },
    {
        "key": "assembly", "title": "Assembly (MEGAHIT)",
        "params": "--assembly_mode (assembly | coassembly | none), --assembly_read_reduction (none | normalize), --megahit_k_list, --megahit_min_count",
        "processes": ["BBNORM", "MEGAHIT", "BBMAP", "CONTIG_FILTER_SUMMARY"],
        "summary": "not summarised yet",
    },
    {
        "key": "binning", "title": "Binning and MAGs (SemiBin2 / MetaBAT2 / COMEBin, Binette, CheckM2, GTDB-Tk)",
        "params": "--include_binning, --binners (Binette refinement with >= 2 binners)",
        "processes": ["SEMIBIN", "METABAT2", "COMEBIN", "BINETTE", "CHECKM2_BATCH", "GTDB_TK_BATCH"],
        "summary": "not summarised yet",
    },
    {
        "key": "contig_arg", "title": "Contig taxonomy and ARGs (BLAST nt, BlobTools, DeepARG)",
        "params": "--contig_tax_and_arg (needs an assembly)",
        "processes": ["NT_BLASTN", "BLOBTOOLS", "DEEPARG_CONTIGS", "ARG_CONTIG_LEVEL_REPORT"],
        "summary": "not summarised yet",
    },
    {
        "key": "bin_arg", "title": "Bin-level ARG clustering (DeepARG on bins, MMseqs2)",
        "params": "--arg_bin_clustering (needs binning)",
        "processes": ["DEEPARG_BINS", "CLUSTERING"],
        "summary": "not summarised yet",
    },
    {
        "key": "contig_functional", "title": "Contig-level functional (Pyrodigal, eggNOG-mapper, dbCAN, featureCounts)",
        "params": "--contig_level_functional (needs an assembly), --functional_cazy",
        "processes": ["PYRODIGAL", "EGGNOG_MAPPER_ANNOTATE", "RUN_DBCAN", "FEATURECOUNTS_GENES", "AGGREGATE_FUNCTIONS"],
        "summary": "not summarised yet",
    },
    {
        "key": "mag_functional", "title": "MAG annotation (Bakta)",
        "params": "--mag_level_functional (needs binning with >= 2 binners)",
        "processes": ["BAKTA_BAKTA"],
        "summary": "not summarised yet",
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
