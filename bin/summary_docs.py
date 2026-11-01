"""Documentation of bin/results_summary.py's outputs: the column dictionary (one source for the Excel README
sheet, README.md / README.docx and LLM_GUIDE.md), the tier explanations, and the README / guide writers
(Markdown, and DOCX from the Python standard library: no python-docx or pandoc needed).
"""

import re
import zipfile
from datetime import date
from xml.sax.saxutils import escape

import pandas as pd

from summary_catalog import STAGES


# ---------------------------------------------------------------------------------------------------------------
# Column dictionary. Per table: a description and its columns as (names, meaning). A name is a column name or a
# token: <metadata> (columns from --metadata), <sample> (one column per sample, in wide tables), <rank> (a
# taxonomic rank column), or re:<regex>. Tables named with <rank> / <ontology> cover a family of tables.
# tests/bin/test_results_summary.sh fails when an emitted column is not documented here.
# ---------------------------------------------------------------------------------------------------------------

SAMPLE_KEY = [("sample", "sample name (the samplesheet's sample column)"),
              ("batch", "results directory the sample came from (--batch-labels / --summary_batch_label)"),
              ("<metadata>", "columns of the --metadata table (renamed meta_<name> when they clash with a summary "
                             "column)")]

WIDE_CELL = ("<sample>", "one column per sample. 0 = not detected in a sample of adequate depth; NA_low_depth = "
                         "not detected in a low-depth sample (unknown, not absent); NA_not_run = the method did not "
                         "run on the sample")

WIDE_RPKG = ("<sample>", "one column per sample: RPKM x clean read pairs / 1e6 / MicrobeCensus genome equivalents. 0 = "
                             "not detected; blank = a call without an RPKM; NA_no_ags = no MicrobeCensus estimate for "
                             "the sample; NA_low_depth / NA_not_run as in the other wide tables")

ASSEMBLY_KEY = [("assembly", "the assembly: a sample name (--assembly_mode assembly), or coassembly"),
                ("batch", "results directory the assembly came from (--batch-labels / --summary_batch_label); with the "
                          "assembly, the key: every co-assembly is called coassembly")]

TABLE_DOCS = {
    "samples_qc": {
        "title": "Per-sample summary",
        "description": "One row per sample: read pairs per step, depth flags, taxonomy and ARG counts, and the "
                       "functional and cell-count estimates of the stages that ran.",
        "columns": SAMPLE_KEY + [
            ("raw_pairs, qc_pass_pairs, clean_pairs", "read pairs before QC, after fastp, after host/contaminant "
                                                      "removal (Reads_report.csv counts single reads; these are "
                                                      "reads / 2)"),
            ("classified_fraction", "fraction of reads the taxonomic profiler classified (sylph: skewed upward when "
                                    "--sylph_read_seq_id is not set and many reads come from genomes absent from "
                                    "the database)"),
            ("pct_removed_qc, pct_removed_decontam, pct_clean", "% of raw read pairs removed by fastp, removed by "
                                                                "host/contaminant filtering, and kept"),
            ("low_depth", "clean_pairs below --min-clean-pairs: absence calls in this sample are not informative"),
            ("depth_class", "clean read pairs binned at --min-clean-pairs and --depth-marks (<1M, 1M-10M, ...)"),
            ("mean_read_length", "mean read length after fastp (mean of R1 and R2), bp"),
            ("taxonomy_db", "taxonomic profiler and database label of the taxonomy tables"),
            ("n_species, shannon", "named species detected and Shannon index (natural log) of their relative "
                                    "abundance, renormalised over them; taxa without a species name are left out"),
            ("arg_status", "ok, NA_low_depth or NA_not_run (no output of the primary ARG method for the sample)"),
            ("n_arg_calls_tier1, n_arg_calls_tier2, n_arg_calls_tier3, n_arg_calls_tier4",
             "calls of the primary ARG method per tier (allele rows; tier 4 = mutation dependent, only "
             "resistance_snp_status card_snp_present rows carry a known resistance mutation)"),
            ("acquired_genes, acquired_rpkm", "distinct tier 1-2 genes of the primary ARG method and their summed "
                                              "gene RPKM"),
            ("n_genes_kmarvel, n_genes_rgi, n_genes_all_methods", "distinct genes per CARD-based ARG method (K-MARVEL, "
                                                                  "RGI at the coverage cutoff) and found by both"),
            ("args_oap_n16S, args_oap_ncell", "ARGs-OAP stage one estimates of 16S copies and of cells (from "
                                              "essential single-copy genes), pairs and singletons summed; nCell is "
                                              "the copies_per_cell denominator of the SARG tables"),
            ("sarg_n_subtypes, sarg_copies_per_cell", "ARGs-OAP: distinct SARG subtypes found and their summed copies "
                                                      "per cell; blank when ARGs-OAP did not run or had no 16S or cell "
                                                      "estimate"),
            ("mc_avg_genome_size_bp, mc_genome_equivalents", "MicrobeCensus average genome size (bp) and genome "
                                                             "equivalents (total bases / AGS), the CPGE denominator; "
                                                             "an estimate of cells, like args_oap_ncell but from "
                                                             "other marker genes, so the two differ"),
            ("func_reads_assigned, func_fraction_annotated", "read functional profile: reads the backend assigned "
                                                             "(Woltka: to WoL ORFs, mates counted separately) and "
                                                             "the share of them with any functional term"),
            ("pct_reads_on_genes", "% of the sample's reads (mates counted separately) featureCounts placed on "
                                   "predicted genes of its assembly (--contig_level_functional; assembly_read_use): "
                                   "the reads the contig-level TPM describes"),
            ("func_contig_genes_annotated_pct, func_contig_annotated_by_abundance_pct",
             "contig-level functional profile: % of the predicted genes with any eggNOG / dbCAN term, and their share "
             "of the sample's TPM"),
        ],
    },
    "taxonomy_long": {
        "title": "Taxonomy, long format",
        "description": "One row per sample and taxon (genome for sylph) with abundance > 0.",
        "columns": SAMPLE_KEY + [
            ("taxon_id", "taxon identifier of the profiler (sylph: GTDB genome accession)"),
            ("<rank>", "lineage of the taxon (GTDB names keep their suffixes, e.g. Enterococcus_B)"),
            ("relative_abundance", "% of the profiled (classified) community; sylph: taxonomic abundance"),
            ("est_read_pairs", "read pairs assigned to the taxon (sylph: sequence_abundance / 100 x clean read "
                               "pairs)"),
            ("sequence_abundance", "sylph: % of all reads, unclassified included"),
            ("adjusted_ani, ani_5_95, true_cov, eff_lambda", "sylph containment ANI, its 5-95% interval (given at "
                                                             "low coverage only; blank = coverage high enough that "
                                                             "sylph reports no interval), estimated true coverage, "
                                                             "and effective lambda (a number, or HIGH)"),
            ("species_call_uncertain", "the lower bound of sylph's ANI interval is below 95%: the genome present may "
                                       "be a different species of the same genus"),
        ],
    },
    "taxonomy_<rank>_wide": {
        "title": "Taxonomy per rank, wide",
        "description": "Relative abundance (%) of each taxon of one rank (sum of its species), one column per sample; "
                       "each column sums to 100 over the classified community.",
        "columns": [("<rank>", "taxon name at this rank ('unassigned <rank>' when the lineage has none)"),
                    ("<sample>", "relative abundance %; 0 = not detected; NA_low_depth = not detected in a low-depth "
                                 "sample. A sample with no taxon detected (or no profiler output) has no column")],
    },
    "taxonomy_species_prevalence": {
        "title": "Species prevalence",
        "description": "Per species: in how many samples it is detected and how abundant it is, over the samples "
                       "with taxonomy output. Low-depth samples count only where the species is detected; they are "
                       "left out of the denominator.",
        "columns": [
            ("<rank>", "species name"),
            ("n_samples_detected", "samples where the species is detected"),
            ("n_samples_evaluable", "samples with taxonomy output where absence is informative (adequate depth), "
                                    "plus low-depth samples where the species is detected: the denominator of "
                                    "pct_samples_detected"),
            ("re:^n_.+_detected$", "groups (--group-by) with at least one detection"),
            ("pct_samples_detected", "n_samples_detected / n_samples_evaluable x 100"),
            ("mean_abundance", "mean relative abundance % over the evaluable samples"),
            ("median_abundance_when_detected, max_abundance", "median over the samples where detected, and maximum, "
                                                              "relative abundance %"),
        ],
    },
    "watchlist": {
        "title": "Watchlist taxa",
        "description": "Relative abundance (%) of each watchlist taxon (assets/summary_watchlist.tsv or "
                       "--watchlist): species rows match the species, genus rows every species of the genus, GTDB "
                       "suffixes removed.",
        "columns": [("taxon, rank, note", "the watchlist entry and its source"),
                    ("matched_names", "profiler names that matched (GTDB names with suffixes)"),
                    ("n_samples_detected", "samples where the taxon is detected"),
                    ("<sample>", "relative abundance %; 0 = not detected; NA_low_depth = not detected in a low-depth "
                                 "sample. A sample with no taxon detected has no column")],
    },
    "arg_long": {
        "title": "ARG calls of the CARD-based methods, long format",
        "description": "Every K-MARVEL call and every RGI ARO at or above the coverage cutoff, one row per sample "
                       "and ARO (allele). The primary method's rows come first.",
        "columns": SAMPLE_KEY + [
            ("method", "kmarvel or rgi; the ARG matrices use the primary method only"),
            ("tier, category, tier_note", "reporting tier from the tier rule table (1 clinical priority, 2 other "
                                          "acquired, 3 intrinsic/efflux/regulator, 4 mutation dependent), the rule's "
                                          "category and note"),
            ("gene", "allele numbers collapsed (KPC-2 -> KPC; OXA alleles -> their CARD OXA subgroup)"),
            ("aro_full_name, aro_id, gene_family, drug_class, resistance_mechanism, model",
             "CARD annotation of the ARO: name, accession, AMR gene family, drug classes (';' separated), resistance "
             "mechanism and CARD model type"),
            ("identity, median_kmer_depth", "K-MARVEL: identity of the assembled gene to the CARD reference (%) and "
                                            "the median k-mer depth; blank for RGI"),
            ("assembly_coverage", "K-MARVEL: % of the reference covered by the assembled gene; RGI: % of the ARO's "
                                  "best-covered reference sequence covered by reads"),
            ("supporting_reads", "read pairs: K-MARVEL supporting read pairs; RGI: mapped reads / 2, summed over the "
                                 "ARO's reference sequences (fractional when the reads are odd)"),
            ("gene_length_nt", "reference length in nt: K-MARVEL: CARD reference protein x 3; RGI: Reference Length "
                               "of the best-covered sequence"),
            ("rpm, rpkm", "supporting read pairs per million clean read pairs; RPM per kb of gene_length_nt"),
            ("gene_rpkm", "RPKM of the gene in the sample, on its gene_representative row (K-MARVEL: the "
                          "best-supported allele; RGI: summed over all the gene's AROs with reads, those below the "
                          "coverage cutoff included, see rgi_gene_read_pairs)"),
            ("resistance_snp_status", "tier 4 only: card_snp_present (a CARD snps.txt resistance mutation is among "
                                      "K-MARVEL's mutation_data), no_card_snp, no_snp_definition, not_checked (no "
                                      "--kmarvel-db), not_checked_rgi (RGI bwt does not report mutations). CARD lists "
                                      "some common polymorphisms (E. coli MarR Y137H), so card_snp_present is "
                                      "necessary, not sufficient, evidence of resistance"),
            ("resistance_snps_found", "the CARD resistance mutation(s) found"),
            ("kmer_blast_discordant", "K-MARVEL: its k-mer assignment and its BLAST best hit name different AROs "
                                      "(blank for RGI)"),
            ("blast_predicted_aro_name", "K-MARVEL: the ARO of the BLAST best hit"),
            ("gene_representative", "True for the row that represents the gene in the sample (K-MARVEL: the "
                                    "best-supported allele, as alleles of one gene share their reads; RGI: the "
                                    "detected, best-covered ARO); matrices and figures use these rows only"),
            ("card_plasmid_associated", "CARD prevalence lists this ARO on NCBI plasmids (with --card_prevalence)"),
            ("candidate_hosts", "species detected in the same sample that CARD prevalence lists as carriers of this "
                                "ARO, with their abundance. Co-occurrence only: reads are not linked to a host"),
            ("rgi_mapped_reads", "RGI: mapped reads (mates counted separately) summed over the ARO's sequences"),
            ("rgi_best_ref", "RGI: source of the best-covered reference sequence (CARD, or WildCARD and the species "
                             "it came from)"),
            ("rgi_gene_coverage", "RGI: best coverage of any ARO of the gene"),
            ("rgi_gene_read_pairs", "RGI: read pairs of the gene, summed over all its AROs with reads (the read "
                                    "count behind gene_rpkm and arg_method_comparison's rgi_read_pairs)"),
            ("rgi_kmer_hosts, rgi_host_check", "RGI k-mer pathogen of origin (top hosts, reads) and whether the top "
                                               "host is a watchlist taxon (watchlist_host, non_watchlist_host, "
                                               "genus_level_host, no_kmer_host, no_kmer_data)"),
            ("mutation_data", "K-MARVEL: every amino-acid difference from the reference (not only resistance "
                              "mutations)"),
        ],
    },
    "arg_gene_rpkm": {
        "title": "ARG genes x samples (RPKM)",
        "description": "Primary ARG method: RPKM of each gene (gene_rpkm), one column per sample.",
        "columns": [("gene", "gene (allele numbers collapsed)"), ("tier, category", "tier and category (arg_long)"),
                    WIDE_CELL],
    },
    "arg_family_rpkm": {
        "title": "ARG gene families x samples (RPKM)",
        "description": "Primary ARG method: gene RPKM summed over each CARD AMR gene family.",
        "columns": [("gene_family", "CARD AMR gene family"), WIDE_CELL],
    },
    "arg_presence": {
        "title": "ARG presence / absence",
        "description": "Primary ARG method: 1 detected, 0 not detected, NA_low_depth / NA_not_run unknown.",
        "columns": [("gene", "gene"), ("tier", "tier"), WIDE_CELL],
    },
    "arg_drug_class_rpkm": {
        "title": "ARG drug classes x samples (RPKM)",
        "description": "Primary ARG method: gene RPKM per drug class, split by tier group. A multi-class gene counts in "
                       "each class, so the rows are not additive.",
        "columns": [("tier_group", "acquired (tiers 1-2), intrinsic/efflux (tier 3), or CARD mutation present "
                                   "(tier 4 rows with resistance_snp_status card_snp_present)"),
                    ("drug_class", "CARD drug class"), WIDE_CELL],
    },
    "arg_priority_panel": {
        "title": "Tier 1 genes x samples (RPKM)",
        "description": "Primary ARG method: the tier 1 (clinical priority) rows of arg_gene_rpkm. Header only when the "
                       "primary method found no tier 1 gene; tier 1 genes found only by another method are in "
                       "arg_method_comparison.",
        "columns": [("gene", "gene"), ("category", "tier 1 category"), WIDE_CELL],
    },
    "group_arg_detection": {
        "title": "Tier 1-2 genes per group",
        "description": "Per --group-by group and tier 1-2 gene of the primary method: detections and the pattern in "
                       "--order-by order.",
        "columns": [("re:.*", "")],  # replaced at run time: the first column is the --group-by column
    },
    "arg_method_comparison": {
        "title": "K-MARVEL vs RGI, per sample and gene",
        "description": "Genes (alleles collapsed) detected by either CARD-based method on the samples both ran on.",
        "columns": SAMPLE_KEY + [
            ("gene, tier, category, gene_family", "gene and its annotation"),
            ("detected_by", "K-MARVEL + RGI, K-MARVEL only or RGI only"),
            ("detected_kmarvel, detected_rgi", "True when the method detected the gene (RGI: >= the coverage cutoff)"),
            ("kmarvel_read_pairs", "K-MARVEL supporting read pairs of the gene's best-supported allele"),
            ("rgi_read_pairs", "RGI read pairs of the gene, summed over its AROs (rgi_gene_read_pairs)"),
            ("rgi_gene_coverage", "best coverage of any RGI ARO of the gene (0 = RGI mapped no read to it)"),
            ("rgi_below_cutoff", "a K-MARVEL gene RGI saw with reads below the coverage cutoff: partial coverage, "
                                 "not absence"),
            ("rgi_best_ref, rgi_kmer_hosts, rgi_host_check", "as in arg_long (blank when RGI mapped no read)"),
        ],
    },
    "arg_method_overlap": {
        "title": "K-MARVEL vs RGI overlap counts",
        "description": "Distinct genes, and sample x gene detections, by the methods that found them: overall, per "
                       "tier and per category.",
        "columns": [("unit", "genes, or detections (sample x gene)"), ("level", "all, tier or category"),
                    ("group", "the tier or category (all for level all)"),
                    ("re:.+ (only|\\+ .+)$", "count found by that combination of methods only"),
                    ("re:K-MARVEL only, RGI below cutoff", "of the K-MARVEL-only count: RGI saw reads below its cutoff"),
                    ("total", "all genes / detections of the row")],
    },
    "arg_sarg_long": {
        "title": "ARGs-OAP (SARG) subtypes, long format",
        "description": "ARGS_OAP_REPORT's SARG subtypes, one row per sample and subtype found. A sample's singletons "
                       "are added to its pairs. Gene-level rows are in the run's "
                       "args_oap_summary/args_oap_long.tsv.",
        "columns": SAMPLE_KEY + [
            ("sarg_type, sarg_subtype", "SARG type (antibiotic class, e.g. tetracycline) and subtype (type__gene "
                                        "group, e.g. tetracycline__tetQ)"),
            ("read_count", "reads assigned to the subtype (single reads; a read on a two- or three-component "
                           "reference counts 1/2 or 1/3)"),
            ("copies_per_cell", "SARG copies (each read weighted by aligned length / reference length) / ARGs-OAP "
                                "nCell (samples_qc args_oap_ncell); blank when nCell is 0"),
            ("copies_per_16S", "SARG copies / ARGs-OAP 16S copies (args_oap_n16S)"),
            ("ppm", "read_count per million single reads (ARGs-OAP nRead)"),
        ],
    },
    "arg_sarg_type_cpc": {
        "title": "ARGs-OAP SARG types x samples (copies per cell)",
        "description": "Copies per cell summed over each SARG type's subtypes. Blank: the sample had no 16S or cell "
                       "estimate, so nothing was normalised.",
        "columns": [("sarg_type", "SARG type"), WIDE_CELL],
    },
    "arg_sarg_subtype_cpc": {
        "title": "ARGs-OAP SARG subtypes x samples (copies per cell)",
        "description": "Copies per cell of each SARG subtype. Blank as in arg_sarg_type_cpc.",
        "columns": [("sarg_type, sarg_subtype", "SARG type and subtype"), WIDE_CELL],
    },
    "arg_contig_mobile": {
        "title": "Contig ARGs on plasmids, viruses and proviruses (geNomad)",
        "description": "DeepARG genes on the assembled contigs (--contig_tax_and_arg) counted by where geNomad "
                       "(--genomad_prediction) places them, per assembly and ARG class. geNomad's plasmid class "
                       "also takes integrative conjugative elements and transposon cargo (e.g. tetM, ermG-mefA, "
                       "cfxA2), which it has no class for: 'on plasmid' means mobile-element contig, not "
                       "necessarily a plasmid. The per-contig calls are in the run's "
                       "contig_level/summary/Contig_tax_and_arg_prediction.tsv (mobile_class, plasmid_score, "
                       "virus_score, provirus).",
        "columns": [
            ("assembly", "the assembly the contigs come from: a sample name, or coassembly"),
            ("batch", "results directory the assembly came from (--batch-labels / --summary_batch_label); with "
                      "the assembly, the key: every co-assembly is called coassembly"),
            ("arg_class", "DeepARG predicted ARG class"),
            ("arg_genes", "ARG genes of the class on the assembly's contigs"),
            ("arg_genes_on_plasmid", "of those, genes on contigs geNomad calls plasmid"),
            ("arg_genes_on_virus", "genes on contigs geNomad calls virus (whole contig)"),
            ("arg_genes_in_provirus", "genes overlapping a provirus region geNomad found inside a contig"),
            ("pct_mobile", "% of arg_genes on a plasmid or virus contig or in a provirus (each gene counted once)"),
        ],
    },
    "function_read_long": {
        "title": "Read-level functional profile, long format",
        "description": "AGGREGATE_READ_FUNCTIONS' read_function_abundance.tsv of every batch: one row per sample and "
                       "term. The terms of one ontology do not add up to the annotated reads (an ORF can carry "
                       "several terms).",
        "columns": SAMPLE_KEY + [
            ("backend", "woltka, superfocus or humann"),
            ("ontology", "ko, ec, cog (COG functional category letters), pfam, metacyc (Woltka); seed_level1-3 "
                         "(SUPER-FOCUS); ko, ec, metacyc (HUMAnN)"),
            ("accession, term_description", "term and its name (blank where the backend gives none; COG "
                                             "categories named from NCBI COG)"),
            ("abundance_cpge", "copies per genome equivalent: RPK / MicrobeCensus genome equivalents (blank without "
                               "MicrobeCensus, or for SUPER-FOCUS)"),
            ("abundance_native, native_unit", "abundance in the backend's unit (Woltka, SUPER-FOCUS: reads, mates "
                                              "counted separately and fractional for multi-hit reads; HUMAnN: RPK)"),
        ],
    },
    "function_read_<ontology>": {
        "title": "Read-level functional profile per ontology, wide",
        "description": "One row per term, one column per sample, in copies per genome equivalent (CPGE) when the "
                       "backend computes it, else in its native unit (the README and the sheet's unit say which). "
                       "Rows sorted by mean abundance.",
        "columns": [("accession", "term"), ("description", "term name (COG categories named from NCBI COG)"),
                    ("<sample>", "abundance; 0 = term not found; NA_no_ags = no MicrobeCensus estimate for the "
                                 "sample, so no CPGE; NA_low_depth / NA_not_run as in the other wide tables")],
    },
    "arg_gene_rpkg": {
        "title": "ARG genes x samples (RPKG)",
        "description": "arg_gene_rpkm per MicrobeCensus genome equivalent: supporting read pairs per kb of gene per "
                       "genome equivalent (an estimate of cells), so samples of different average genome size compare. "
                       "Written when MicrobeCensus ran.",
        "columns": [("gene", "gene (allele numbers collapsed)"), ("tier, category", "tier and category (arg_long)"),
                    WIDE_RPKG],
    },
    "arg_family_rpkg": {
        "title": "ARG gene families x samples (RPKG)",
        "description": "arg_family_rpkm per MicrobeCensus genome equivalent.",
        "columns": [("gene_family", "CARD AMR gene family"), WIDE_RPKG],
    },
    "arg_drug_class_rpkg": {
        "title": "ARG drug classes x samples (RPKG)",
        "description": "arg_drug_class_rpkm per MicrobeCensus genome equivalent. A multi-class gene counts in each "
                       "class, so the rows are not additive.",
        "columns": [("tier_group", "acquired (tiers 1-2), intrinsic/efflux (tier 3), or CARD mutation present "
                                   "(tier 4 rows with resistance_snp_status card_snp_present)"),
                    ("drug_class", "CARD drug class"), WIDE_RPKG],
    },
    "assembly_stats": {
        "title": "Assemblies",
        "description": "One row per assembly (MEGAHIT): BBMap statistics of all its contigs, and how many passed the "
                       "length filter that every contig-based step uses (binning, contig taxonomy, ARGs, functions). "
                       "MEGAHIT is not reproducible run to run: the same reads give slightly different contigs.",
        "columns": ASSEMBLY_KEY + [
            ("assembly_mode", "assembly (one per sample) or coassembly (all samples pooled)"),
            ("contigs, assembly_length_mb", "contigs MEGAHIT wrote and their total length (Mb)"),
            ("n50_kb, l50", "N50 (kb): half of the assembly is in contigs at least this long; L50: how many contigs "
                            "that takes (BBMap prints them the other way round, as N/L50 = count/length)"),
            ("max_contig_kb", "longest contig (kb)"),
            ("gc_pct", "GC content of the contigs (%)"),
            ("min_contig_length_bp, contigs_kept, pct_contigs_kept", "length filter (--bbmap_length) and the contigs "
                                                                     "it kept, as a number and % of contigs"),
        ],
    },
    "assembly_read_use": {
        "title": "Reads on the assembly",
        "description": "Per sample, the share of its reads featureCounts placed on the contigs of its assembly (a "
                       "co-assembly for every sample under --assembly_mode coassembly) and on their predicted genes. "
                       "Reads not on contigs are outside every contig-based table. NA_not_run without "
                       "--contig_level_functional (featureCounts is part of it).",
        "columns": SAMPLE_KEY + [
            ("assembly", "the assembly the sample's reads were aligned to: the sample itself, or coassembly"),
            ("pct_reads_on_contigs", "% of the sample's reads (mates counted separately) aligned to a filtered contig"),
            ("pct_reads_on_genes", "% of its reads assigned to a predicted gene: the denominator of contig-level TPM"),
        ],
    },
    "mags": {
        "title": "MAGs (refined bins)",
        "description": "One row per refined MAG (Binette, or the single binner's bins): CheckM2 quality, MIMAG tier, "
                       "GTDB-Tk lineage and, with --mag_level_functional, Bakta's rRNA, tRNA and gene counts. MAGs of "
                       "different per-sample assemblies are not dereplicated: one species can appear once per sample.",
        "columns": ASSEMBLY_KEY + [
            ("bin", "MAG name in its assembly (binette_binN)"),
            ("binner", "binette (refinement of >= 2 binners) or the single binner"),
            ("completeness_pct, contamination_pct", "CheckM2 estimates (%)"),
            ("checkm2_model", "CheckM2 model behind the completeness (specific neural network or general gradient "
                              "boost)"),
            ("genome_size_mb, contigs, contig_n50_kb, gc_pct, coding_density", "CheckM2: MAG size (Mb), contigs, "
                                                                               "contig N50 (kb), GC (%), coding "
                                                                               "density"),
            ("quality_tier", "MIMAG (Bowers et al. 2017): high = completeness > 90, contamination < 5, 23S + 16S + 5S "
                             "rRNA and >= 18 tRNA amino acids; medium = completeness >= 50, contamination < 10; low = "
                             "the rest. Stricter on rRNA/tRNA than the High/Mid/Low of the run's 04_binning plots"),
            ("tier_basis", "checkm2_rrna_trna: Bakta's rRNA / tRNA checked; checkm2_only: no Bakta, so at most medium"),
            ("domain, phylum, class, order, family, genus, species", "GTDB-Tk classification (GTDB names; blank where "
                                                                    "it stops)"),
            ("gtdb_closest_reference", "GTDB-Tk's closest reference genome (blank when there is none)"),
            ("rrna_16s, rrna_23s, rrna_5s", "Bakta found the rRNA gene (truncated genes at a contig end count); blank "
                                            "without Bakta"),
            ("trna_aa_types", "standard amino acids with a Bakta tRNA (pseudo tRNAs not counted)"),
            ("cds, pct_hypothetical, crispr_arrays", "Bakta: coding sequences, % of them hypothetical proteins, CRISPR "
                                                     "arrays"),
        ],
    },
    "mag_abundance": {
        "title": "MAG depth per sample",
        "description": "Read depth of each MAG in each sample and the share of the sample's sequenced bases it holds. "
                       "Under co-assembly every MAG has a row per sample (BEDTOOLS_BIN_DEPTH: the mean of its contigs' "
                       "mean depths, each contig weighted equally); under per-sample assembly a MAG has a row for its "
                       "own sample only (length-weighted contig coverage from the contig table, --contig_tax_and_arg; "
                       "NA_not_run without it).",
        "columns": [
            ("assembly", "the assembly the MAG comes from: a sample name, or coassembly"),
            ("batch", "results directory (--batch-labels / --summary_batch_label)"),
            ("bin", "MAG name (mags)"),
            ("sample", "sample whose reads give the depth"),
            ("<metadata>", "columns of the --metadata table"),
            ("mean_depth_x", "mean read depth on the MAG's contigs (raw, not normalised)"),
            ("pct_sample_bases", "depth x MAG size / bases sequenced in the sample (MicrobeCensus total bases, else "
                                 "clean pairs x 2 x mean read length) x 100: comparable between samples, unlike raw "
                                 "depth"),
        ],
    },
    "arg_contig_long": {
        "title": "Contig ARGs with their taxon and bin",
        "description": "One row per DeepARG hit on an assembled contig (--contig_tax_and_arg), with the contig's "
                       "Metabuli taxon, bin and geNomad class. This is linkage (the ARG and the taxon are on one "
                       "contig), unlike read-level candidate hosts. Counts are genes on contigs, not read abundances; "
                       "ARGs of low-coverage or unassembled genomes are missing. DeepARG calls are predictions: many "
                       "have < 80 % identity to their reference.",
        "columns": ASSEMBLY_KEY + [
            ("contig, contig_length_bp", "contig name and length"),
            ("contig_depth_x", "BlobTools mean read depth of the contig (raw; under co-assembly pooled over all "
                               "samples)"),
            ("arg_gene, arg_class", "DeepARG ARG and its predicted class"),
            ("deeparg_probability, identity_pct", "DeepARG probability and % identity to the best SARG-DeepARG hit"),
            ("host_rank, host_taxon", "lowest of the seven ranks Metabuli reached for the contig and its name "
                                      "(unclassified when none)"),
            ("phylum, genus, species", "Metabuli's phylum, genus and species of the contig (NCBI names)"),
            ("taxonomy_score", "Metabuli's score of the call"),
            ("bin, mag_tier, bin_gtdb_species", "MAG the contig is binned in, its MIMAG tier (mags) and GTDB species; "
                                                "blank for unbinned contigs"),
            ("mobile_class", "geNomad: chromosome, plasmid or virus (blank without --genomad_prediction)"),
        ],
    },
    "arg_contig_hosts": {
        "title": "Contig ARGs per host taxon",
        "description": "ARG genes on contigs per assembly, ARG class and host: the contig's genus, or the lowest rank "
                       "above it that Metabuli reached (unclassified when none). Gene counts, not abundances.",
        "columns": ASSEMBLY_KEY + [
            ("arg_class", "DeepARG predicted ARG class"),
            ("host_rank, host_taxon", "rank (genus or above) and name of the host"),
            ("arg_genes, arg_contigs", "ARG genes of the class on contigs of the host, and those contigs"),
            ("on_mag", "of the genes, those on a binned contig"),
            ("on_mobile", "of the genes, those on a geNomad plasmid or virus contig"),
        ],
    },
    "arg_contig_abundance": {
        "title": "Contig ARG genes per sample",
        "description": "Each DeepARG gene's read abundance in each sample, from featureCounts on its predicted gene "
                       "(--contig_level_functional): under co-assembly one row per gene and sample, under per-sample "
                       "assembly the gene's own sample only. TPM and CPGE as in the contig-level functional tables.",
        "columns": SAMPLE_KEY + [
            ("assembly", "the assembly the gene was predicted on: a sample name, or coassembly"),
            ("gene_id, contig", "Pyrodigal gene and its contig"),
            ("arg_gene, arg_class", "DeepARG ARG and its predicted class"),
            ("count", "reads featureCounts assigned to the gene (mates counted separately)"),
            ("tpm", "transcripts-per-million style share of the sample's reads on predicted genes (design doc 6.1)"),
            ("cpge", "copies per genome equivalent: RPK / MicrobeCensus genome equivalents (blank without it)"),
        ],
    },
    "arg_mag_long": {
        "title": "ARG genes of MAGs and their clusters",
        "description": "One row per ARG gene on a MAG (--arg_bin_clustering: DeepARG hits of the assembly split by "
                       "bin) and its MMseqs2 cluster at 90, 95, 99 and 100 % identity. A cluster is named after its "
                       "representative gene; genes of one cluster share that sequence identity.",
        "columns": ASSEMBLY_KEY + [
            ("bin", "MAG (mags)"),
            ("gene_id", "Pyrodigal gene"),
            ("arg_gene, arg_class", "DeepARG ARG and its predicted class"),
            ("cluster_90, cluster_95, cluster_99, cluster_100", "cluster at each identity, as "
                                                                "<assembly>:<bin>:<gene of the representative>"),
        ],
    },
    "arg_cluster_sharing": {
        "title": "ARG clusters shared by MAGs",
        "description": "Clusters at 90 % identity found in more than one MAG. The same ARG sequence in MAGs of "
                       "different samples or taxa is consistent with transfer but not evidence of it on its own. Only "
                       "ARGs on binned contigs are here.",
        "columns": [
            ("cluster", "cluster (arg_mag_long cluster_90)"),
            ("batch", "results directory (clusters are made within one run)"),
            ("arg_gene, arg_class", "DeepARG ARGs and classes of the members"),
            ("n_mags, n_assemblies", "MAGs and assemblies holding the cluster"),
            ("mags", "the MAGs, as <assembly>:<bin>"),
        ],
    },
    "function_contig_long": {
        "title": "Contig-level functional profile, long format",
        "description": "AGGREGATE_FUNCTIONS' function_abundance.tsv (eggNOG-mapper and dbCAN on the predicted genes of "
                       "the assemblies, --contig_level_functional): one row per sample and term. A gene with two "
                       "terms of an ontology counts fully for each, so an ontology's terms add up to more than the "
                       "genes. Reported apart from the read-level profile (function_read_*), never merged with it.",
        "columns": SAMPLE_KEY + [
            ("backend", "eggnog-mapper or run_dbcan"),
            ("ontology", "ko, cog (COG category letters), ec, pfam (Pfam names), cazy (eggNOG), cazy_dbcan (dbCAN "
                         "families and eCAMI subfamilies)"),
            ("accession, term_description", "term and its name: COG categories from NCBI COG, eggNOG CAZy families "
                                             "from eggNOG; blank for KO, EC, Pfam and dbCAN"),
            ("abundance_tpm", "TPM: RPK_i = count_i / (length_i / 1000); TPM_i = RPK_i / (sum over all j of RPK_j) * "
                              "1e6, summed over the term's genes. The share of the sample's functional pool, relative "
                              "to reads on predicted genes (pct_reads_on_genes); compositional"),
            ("abundance_cpge", "copies per genome equivalent: genome_equivalents = total_bases_sampled / "
                               "average_genome_size_bp (MicrobeCensus); CPGE_i = RPK_i / genome_equivalents. Not "
                               "compositional; blank without MicrobeCensus"),
        ],
    },
    "function_contig_<ontology>": {
        "title": "Contig-level functional profile per ontology, wide",
        "description": "One row per term, one column per sample, in copies per genome equivalent (CPGE) when "
                       "MicrobeCensus ran, else in TPM (the README and the sheet's unit say which). Rows sorted by "
                       "mean abundance. Not comparable row by row with function_read_<ontology>: other genes, other "
                       "reads, other units.",
        "columns": [("accession", "term"), ("description", "term name (COG categories, eggNOG CAZy families; blank "
                                                            "for KO, EC, Pfam and dbCAN)"),
                    ("<sample>", "abundance; 0 = term not found; NA_no_ags = no MicrobeCensus estimate for the "
                                 "sample, so no CPGE; NA_low_depth / NA_not_run as in the other wide tables")],
    },
}


def table_doc(name, group_by=None):
    """The TABLE_DOCS entry of a table name, families (<rank>, <ontology>) resolved."""
    if name in TABLE_DOCS:
        doc = TABLE_DOCS[name]
    elif re.fullmatch(r"taxonomy_\w+_wide", name):
        doc = TABLE_DOCS["taxonomy_<rank>_wide"]
    elif name.startswith("function_read_") and name != "function_read_long":
        doc = TABLE_DOCS["function_read_<ontology>"]
    elif name.startswith("function_contig_") and name != "function_contig_long":
        doc = TABLE_DOCS["function_contig_<ontology>"]
    else:
        return None
    if name == "group_arg_detection":
        doc = dict(doc, columns=[
            (group_by or "<group>", "the --group-by group"), ("gene, tier, category", "gene and its annotation"),
            ("n_detected, n_evaluable, n_samples", "samples of the group with the gene; samples where absence is "
                                                   "informative (or the gene was found); all samples of the group"),
            ("first_detected, last_detected", "first and last sample (by --order-by) with the gene"),
            ("pattern", "detection per sample in --order-by order (+ detected, - not detected, NA unknown)")])
    return doc


RANKS = {"Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain",
         "domain", "kingdom", "phylum", "class", "order", "family", "genus", "species", "strain"}


def _names(spec):
    return [n.strip() for n in spec.split(",")] if not spec.startswith("re:") else [spec]


def column_meaning(name, column, meta_cols, samples, group_by=None):
    """The documented meaning of a table column, or None when it is not documented."""
    doc = table_doc(name, group_by)
    if doc is None:
        return None
    for spec, meaning in doc["columns"]:
        for n in _names(spec):
            if n == column:
                return meaning
            if n == "<metadata>" and column in meta_cols:
                return meaning
            if n == "<sample>" and column in samples:
                return meaning
            if n == "<rank>" and column in RANKS:
                return meaning
            if n.startswith("re:") and re.fullmatch(n[3:], column):
                return meaning
    return None


def undocumented_columns(tables, meta_cols, samples, group_by=None):
    """[(table, column)] of emitted columns (index included) that TABLE_DOCS does not document."""
    out = []
    for name, (df, index) in tables.items():
        cols = list(df.index.names) + list(df.columns) if index else list(df.columns)
        for c in cols:
            if c is None:
                continue
            if column_meaning(name, str(c), set(meta_cols), set(samples), group_by) is None:
                out.append((name, str(c)))
    return out


def unused_doc_columns(tables, group_by=None):
    """[(table, documented name)] of TABLE_DOCS column names (tokens excepted) that no emitted column of the table
    matches: stale entries, or columns this run did not produce."""
    out = []
    for name, (df, index) in tables.items():
        doc = table_doc(name, group_by)
        if doc is None:
            continue
        cols = {str(c) for c in (list(df.index.names) if index else []) + list(df.columns) if c is not None}
        for spec, _ in doc["columns"]:
            for n in _names(spec):
                if n.startswith("<") or n.startswith("re:"):
                    continue
                if n not in cols:
                    out.append((name, n))
    return out


# ---------------------------------------------------------------------------------------------------------------
# ARG tiers (shared by the HTML reading guide and the README)
# ---------------------------------------------------------------------------------------------------------------

TIER_MEANING = {
    1: ("Clinical priority", "last-line or highly transmissible resistance, usually on mobile elements"),
    2: ("Other acquired", "acquired genes of common resistance (ESBL, AmpC, aminoglycoside-modifying enzymes, "
                          "tetracycline, macrolide, sulfonamide/trimethoprim genes, ...)"),
    3: ("Intrinsic, efflux or regulator", "chromosomal or housekeeping genes that most strains of a species carry; "
                                          "they say little about acquired resistance"),
    4: ("Mutation dependent", "genes that confer resistance only with a specific mutation (gyrA, parC, rpoB, ...)"),
}


def tier_categories(rules_path):
    """{tier: ['category (note)', ...]} from the tier rule file, in rule order, plus the built-in rules."""
    out = {1: [], 2: [], 3: [], 4: ["mutation_dependent (CARD variant models and AROs with CARD SNP definitions)"]}
    try:
        rules = pd.read_csv(rules_path, sep="\t", comment="#", dtype=str, keep_default_na=False)
    except (OSError, ValueError):
        rules = pd.DataFrame(columns=["tier", "category", "note"])
    if "note" not in rules.columns:  # optional, as in load_rules
        rules["note"] = ""
    seen = set()
    for tier, cat, note in zip(rules.get("tier", []), rules.get("category", []), rules["note"]):
        key = (int(tier), cat)
        if key in seen:
            continue
        seen.add(key)
        out.setdefault(int(tier), []).append(f"{cat} ({note})" if note else cat)
    out[3].append("efflux_permeability (efflux pumps, porins and overexpression models without a rule)")
    out[2].append("acquired_unreviewed (no rule matched)")
    return out


# ---------------------------------------------------------------------------------------------------------------
# Document model: [("h", level, text) | ("p", text) | ("ul", [text]) | ("table", [header], [[cell]]) |
# ("code", text)], rendered to Markdown and to DOCX
# ---------------------------------------------------------------------------------------------------------------

def _md_cell(text):
    return str(text).replace("|", "\\|").replace("\n", " ")


def render_markdown(blocks):
    out = []
    for b in blocks:
        kind = b[0]
        if kind == "h":
            out.append("#" * b[1] + " " + b[2])
        elif kind == "p":
            out.append(b[1])
        elif kind == "ul":
            out.append("\n".join(f"- {item}" for item in b[1]))
        elif kind == "table":
            head, rows = b[1], b[2]
            lines = ["| " + " | ".join(_md_cell(h) for h in head) + " |", "|" + "---|" * len(head)]
            lines += ["| " + " | ".join(_md_cell(c) for c in r) + " |" for r in rows]
            out.append("\n".join(lines))
        elif kind == "code":
            out.append("```\n" + b[1].rstrip("\n") + "\n```")
    return "\n\n".join(out) + "\n"


_W = 'xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main"'


def _run(text, bold=False, mono=False, size=None):
    props = ""
    if bold or mono or size:
        props = "<w:rPr>" + ("<w:b/>" if bold else "") + \
                ('<w:rFonts w:ascii="Consolas" w:hAnsi="Consolas" w:cs="Consolas"/>' if mono else "") + \
                (f'<w:sz w:val="{size}"/><w:szCs w:val="{size}"/>' if size else "") + "</w:rPr>"
    return f'<w:r>{props}<w:t xml:space="preserve">{escape(str(text))}</w:t></w:r>'


def _inline(text, size=None):
    """**bold** and `code` spans of a Markdown line as runs."""
    runs = []
    for part in re.split(r"(\*\*[^*]+\*\*|`[^`]+`)", str(text)):
        if not part:
            continue
        if part.startswith("**") and part.endswith("**"):
            runs.append(_run(part[2:-2], bold=True, size=size))
        elif part.startswith("`") and part.endswith("`"):
            runs.append(_run(part[1:-1], mono=True, size=size))
        else:
            runs.append(_run(part, size=size))
    return "".join(runs)


def _para(text, style=None, size=None):
    ppr = f'<w:pPr><w:pStyle w:val="{style}"/></w:pPr>' if style else ""
    return f"<w:p>{ppr}{_inline(text, size)}</w:p>"


def render_docx(blocks, path, title):
    """A .docx (WordprocessingML in a zip) with headings, paragraphs, bullet lists, tables and code blocks."""
    body = []
    for b in blocks:
        kind = b[0]
        if kind == "h":
            body.append(_para(b[2], style=f"Heading{min(b[1], 3)}"))
        elif kind == "p":
            body.append(_para(b[1]))
        elif kind == "ul":
            body += [_para(item, style="ListBullet") for item in b[1]]
        elif kind == "code":
            body += [f'<w:p><w:pPr><w:pStyle w:val="Code"/></w:pPr>{_run(line, mono=True, size=16)}</w:p>'
                     for line in b[1].rstrip("\n").split("\n")]
        elif kind == "table":
            head, rows = b[1], b[2]
            grid = "".join("<w:gridCol/>" for _ in head)

            def row(cells, header=False):
                tcs = "".join(f"<w:tc><w:tcPr><w:tcW w:w=\"0\" w:type=\"auto\"/></w:tcPr>"
                              f"<w:p>{_run(c, bold=header, size=16) if header else _inline(c, 16)}</w:p></w:tc>"
                              for c in cells)
                return f"<w:tr>{'<w:trPr><w:tblHeader/></w:trPr>' if header else ''}{tcs}</w:tr>"
            body.append('<w:tbl><w:tblPr><w:tblStyle w:val="TableGrid"/><w:tblW w:w="5000" w:type="pct"/>'
                        f'</w:tblPr><w:tblGrid>{grid}</w:tblGrid>{row(head, True)}'
                        + "".join(row([str(c) for c in r]) for r in rows) + "</w:tbl>")
            body.append("<w:p/>")
    document = (f'<?xml version="1.0" encoding="UTF-8" standalone="yes"?><w:document {_W}><w:body>'
                + "".join(body) + '<w:sectPr><w:pgSz w:w="11906" w:h="16838"/><w:pgMar w:top="1134" w:right="1134" '
                'w:bottom="1134" w:left="1134" w:header="709" w:footer="709" w:gutter="0"/></w:sectPr></w:body>'
                '</w:document>')
    border = "".join(f'<w:{s} w:val="single" w:sz="4" w:space="0" w:color="BFBFBF"/>'
                     for s in ("top", "left", "bottom", "right", "insideH", "insideV"))
    styles = (f'<?xml version="1.0" encoding="UTF-8" standalone="yes"?><w:styles {_W}>'
              '<w:docDefaults><w:rPrDefault><w:rPr><w:rFonts w:ascii="Calibri" w:hAnsi="Calibri" w:cs="Calibri"/>'
              '<w:sz w:val="21"/><w:szCs w:val="21"/></w:rPr></w:rPrDefault><w:pPrDefault><w:pPr>'
              '<w:spacing w:after="100" w:line="264" w:lineRule="auto"/></w:pPr></w:pPrDefault></w:docDefaults>'
              '<w:style w:type="paragraph" w:default="1" w:styleId="Normal"><w:name w:val="Normal"/></w:style>'
              + "".join(f'<w:style w:type="paragraph" w:styleId="Heading{i}"><w:name w:val="heading {i}"/>'
                        f'<w:basedOn w:val="Normal"/><w:next w:val="Normal"/><w:qFormat/><w:pPr><w:keepNext/>'
                        f'<w:spacing w:before="{300 - 60 * i}" w:after="80"/><w:outlineLvl w:val="{i - 1}"/></w:pPr>'
                        f'<w:rPr><w:b/><w:color w:val="1F4E79"/><w:sz w:val="{36 - 4 * i}"/></w:rPr></w:style>'
                        for i in (1, 2, 3))
              + '<w:style w:type="paragraph" w:styleId="ListBullet"><w:name w:val="List Bullet"/>'
              '<w:basedOn w:val="Normal"/><w:pPr><w:numPr><w:ilvl w:val="0"/><w:numId w:val="1"/></w:numPr>'
              '<w:spacing w:after="40"/><w:ind w:left="360" w:hanging="360"/></w:pPr></w:style>'
              '<w:style w:type="paragraph" w:styleId="Code"><w:name w:val="Code"/><w:basedOn w:val="Normal"/>'
              '<w:pPr><w:spacing w:after="0" w:line="240" w:lineRule="auto"/><w:shd w:val="clear" w:fill="F2F2F2"/>'
              '</w:pPr></w:style>'
              f'<w:style w:type="table" w:styleId="TableGrid"><w:name w:val="Table Grid"/><w:tblPr><w:tblBorders>'
              f'{border}</w:tblBorders><w:tblCellMar><w:left w:w="80" w:type="dxa"/><w:right w:w="80" w:type="dxa"/>'
              '</w:tblCellMar></w:tblPr></w:style></w:styles>')
    numbering = (f'<?xml version="1.0" encoding="UTF-8" standalone="yes"?><w:numbering {_W}>'
                 '<w:abstractNum w:abstractNumId="0"><w:lvl w:ilvl="0"><w:start w:val="1"/><w:numFmt w:val="bullet"/>'
                 '<w:lvlText w:val="\u2022"/><w:lvlJc w:val="left"/><w:pPr><w:ind w:left="360" w:hanging="360"/>'
                 '</w:pPr></w:lvl></w:abstractNum><w:num w:numId="1"><w:abstractNumId w:val="0"/></w:num>'
                 '</w:numbering>')
    content_types = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                     '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
                     '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
                     '<Default Extension="xml" ContentType="application/xml"/>'
                     '<Override PartName="/word/document.xml" ContentType="application/vnd.openxmlformats-'
                     'officedocument.wordprocessingml.document.main+xml"/>'
                     '<Override PartName="/word/styles.xml" ContentType="application/vnd.openxmlformats-'
                     'officedocument.wordprocessingml.styles+xml"/>'
                     '<Override PartName="/word/numbering.xml" ContentType="application/vnd.openxmlformats-'
                     'officedocument.wordprocessingml.numbering+xml"/>'
                     '<Override PartName="/docProps/core.xml" ContentType="application/vnd.openxmlformats-'
                     'package.core-properties+xml"/></Types>')
    rels = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
            '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/'
            'officeDocument" Target="word/document.xml"/>'
            '<Relationship Id="rId2" Type="http://schemas.openxmlformats.org/package/2006/relationships/metadata/'
            'core-properties" Target="docProps/core.xml"/></Relationships>')
    doc_rels = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
                '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/'
                'styles" Target="styles.xml"/>'
                '<Relationship Id="rId2" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/'
                'numbering" Target="numbering.xml"/></Relationships>')
    core = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties" '
            'xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:dcterms="http://purl.org/dc/terms/" '
            'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">'
            f'<dc:title>{escape(title)}</dc:title><dc:creator>BugBuster results_summary.py</dc:creator>'
            f'<dcterms:created xsi:type="dcterms:W3CDTF">{date.today().isoformat()}T00:00:00Z</dcterms:created>'
            '</cp:coreProperties>')
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("[Content_Types].xml", content_types)
        z.writestr("_rels/.rels", rels)
        z.writestr("word/_rels/document.xml.rels", doc_rels)
        z.writestr("word/document.xml", document)
        z.writestr("word/styles.xml", styles)
        z.writestr("word/numbering.xml", numbering)
        z.writestr("docProps/core.xml", core)
    return path


# ---------------------------------------------------------------------------------------------------------------
# README and LLM guide
# ---------------------------------------------------------------------------------------------------------------

UNITS = [
    ("read pairs", "every count in the summary tables except ARGs-OAP's read_count and ppm (single reads, as ARGs-OAP "
                   "counts them); Reads_report.csv and ARGs-OAP nRead count single reads (both mates), RGI and Woltka "
                   "count mates separately"),
    ("RPM, RPKM", "K-MARVEL / RGI: supporting read pairs per million clean read pairs; RPM per kb of reference gene. "
                  "Comparable between samples of different depth; not copies per cell"),
    ("CPGE", "functional profiles: reads per kb per MicrobeCensus genome equivalent (RPK / GE). Contig level: "
             "genome_equivalents = total_bases_sampled / average_genome_size_bp (MicrobeCensus); CPGE_i = RPK_i / "
             "genome_equivalents. Approximate average copies of that gene per community member; not compositional"),
    ("TPM", "contig-level functional profile: RPK_i = count_i / (length_i / 1000); TPM_i = RPK_i / (sum over all j of "
            "RPK_j) * 1e6. The share of the sample's functional pool attributable to that gene; compositional, and "
            "relative to the reads on predicted genes only (assembly_read_use pct_reads_on_genes)"),
    ("RPKG", "arg_*_rpkg: K-MARVEL / RGI supporting read pairs per kb of gene per MicrobeCensus genome equivalent "
             "(RPKM x clean read pairs / 1e6 / genome equivalents)"),
    ("relative abundance", "taxonomy: % of the profiled (classified) community; each sample sums to 100. What is "
                           "profiled differs by profiler: Kraken2 + Bracken species read estimates, sourmash's "
                           "abundance-weighted k-mer share of the matched genomes, sylph's taxonomic abundance; their "
                           "classified fractions (samples_qc) are not comparable between profilers"),
    ("MIMAG tier", "mags quality_tier: high / medium / low (Bowers et al. 2017); tier_basis says whether Bakta's "
                   "rRNA / tRNA genes were checked (without them a MAG is at most medium)"),
    ("copies per cell", "ARGs-OAP: SARG copies (reads weighted by the share of the reference they cover) / "
                        "ARGs-OAP cell estimate (nCell)"),
    ("NA_low_depth", "not detected in a sample below --min-clean-pairs: unknown, not absent"),
    ("NA_not_run", "the method did not run on the sample (its output is missing)"),
    ("NA_no_ags", "CPGE table: the sample has no MicrobeCensus estimate"),
]

PITFALLS = [
    "geNomad's plasmid class includes integrative conjugative elements and transposon cargo (arg_contig_mobile): "
    "an ARG 'on a plasmid' sits on a mobile-element contig, not necessarily a plasmid.",
    "Single reads vs pairs: Reads_report.csv, ARGs-OAP nRead, RGI mapped reads and Woltka counts are single reads "
    "(both mates); the summary's counts are read pairs (reads / 2), except ARGs-OAP's read_count and ppm.",
    "ARG units are not interchangeable: K-MARVEL / RGI RPKM (read pairs) and ARGs-OAP copies per cell (SARG copies "
    "over the ARGs-OAP cell estimate) measure different things. Compare samples within one unit only.",
    "Two cell estimates: ARGs-OAP nCell and MicrobeCensus genome equivalents (CPGE) come from different marker genes "
    "and differ by tens of percent; they are not interchangeable denominators.",
    "SARG (ARGs-OAP) and CARD (K-MARVEL/RGI) name genes differently; their calls are reported separately, not "
    "matched gene by gene.",
    "Alleles share reads (K-MARVEL): a gene is represented by one allele (gene_representative); never sum allele "
    "rows.",
    "RGI bwt lists every reference with at least one mapped read: the summary counts an RGI gene only at >= "
    "--rgi_summary_min_coverage of its best reference; the raw RGI tables are not filtered.",
    "Absence in a low-depth sample is NA_low_depth (unknown), not 0.",
    "Functional terms of one ontology overlap (an ORF carries several terms): the tables are not compositional and a "
    "column does not sum to the annotated reads.",
    "Taxonomic relative abundance is relative to the classified community; the classified fraction is in "
    "samples_qc.",
    "Tier 1 calls of RGI from WildCARD sequences may come from commensal relatives of the clinical gene: check "
    "rgi_kmer_hosts / rgi_host_check.",
    "TPM and CPGE are not interchangeable (the most common error with these tables): TPM is a share of the reads on "
    "predicted genes (compositional), CPGE is copies per genome equivalent (not compositional). Contig-level TPM "
    "describes only the reads on predicted genes, which can be well under half of a sample's reads (assembly_read_use).",
    "A gene with two terms of one ontology counts fully for each (intended): an ontology's terms add up to more than "
    "the gene total.",
    "Read-level (function_read_*) and contig-level (function_contig_*) functional tables come from different genes, "
    "reads and units: never merge or compare them row by row.",
    "Contig ARG tables count genes on assembled contigs, not read abundance: ARGs of low-coverage or unassembled "
    "genomes are missing, and DeepARG calls below 80 % identity are common. Under co-assembly the contig depth is "
    "pooled over all samples; arg_contig_abundance gives each sample's reads on the gene.",
    "MAG quality: the mags quality_tier is MIMAG; the High/Mid/Low classes of the run's 04_binning plots and tables "
    "are another, older scheme. MAGs of different per-sample assemblies are not dereplicated.",
    "Co-assembly MAG depth (mag_abundance) is the mean of its contigs' mean depths, each contig weighted equally; "
    "compare samples with pct_sample_bases rather than raw depth.",
    "An ARG cluster shared by MAGs of different samples or taxa (arg_cluster_sharing) is the same sequence at the "
    "stated identity; it is consistent with transfer, not evidence of it by itself.",
]


def stage_rows(stage_status_by_batch):
    """[[stage, status per batch ..., what the summary takes]] for the README / guide stage table."""
    labels = list(stage_status_by_batch)
    rows = []
    for i, st in enumerate(STAGES):
        statuses = [stage_status_by_batch[b][i][1] for b in labels]
        rows.append([st["title"]] + statuses + [st["summary"]])
    return labels, rows


def readme_blocks(ctx):
    """README.md / README.docx content."""
    s = ctx["summary"]
    tables, figs = ctx["tables"], ctx["figures"]
    blocks = [("h", 1, ctx["title"]),
              ("p", f"BugBuster run-level summary generated on {date.today().isoformat()}: {s['n_samples']} samples "
                    f"from {s['n_batches']} results director{'y' if s['n_batches'] == 1 else 'ies'}. Open "
                    "`results_summary.html` for the report with figures; the tables are in `tables/` (TSV) and in "
                    "`results_summary.xlsx` (one sheet each, plus a README sheet with every column). "
                    "`LLM_GUIDE.md` describes this folder for a language model or a script.")]
    labels, rows = stage_rows(ctx["stage_status"])
    blocks += [("h", 2, "What ran"),
               ("p", "Pipeline stages per results directory, and what this summary takes from each. 'not run' stages "
                     "have no output; 'ran, no table found' means the stage ran but its table was not in the results "
                     "directory."),
               ("table", ["Stage"] + labels + ["In this summary"], rows)]
    blocks += [("h", 2, "Samples"),
               ("ul", [f"{s['n_samples']} samples; {s['n_low_depth']} below {ctx['min_pairs']:,} clean read pairs "
                       "(low depth: absence there is unknown, NA_low_depth)"]
                + [f"{n} samples with >= {m:,} clean read pairs" for m, n in s.get("depth_marks", [])])]
    if s.get("arg_methods"):
        m = [f"CARD-based ARG methods: {', '.join(s['arg_methods'])}; the ARG matrices and figures use "
             f"**{s.get('arg_primary', '')}**."]
        if s.get("n_tier1") is not None:
            m.append(f"{s['n_tier1']} samples with a tier 1 (clinical priority) gene of {s.get('arg_primary', '')}; "
                     f"{s.get('n_tier1_genes', 0)} distinct tier 1 genes.")
        if s.get("n_tier1_genes_all_methods") is not None:
            m.append(f"{s['n_tier1_genes_all_methods']} tier 1 genes found by both CARD-based methods.")
        blocks += [("h", 2, "Antimicrobial resistance genes"), ("ul", m)]
    if s.get("attention"):
        blocks += [("h", 2, "Attention"), ("ul", [text for _, text, _ in s["attention"]])]
    blocks += [("h", 2, "Tables")]
    trows = []
    for name, (df, _) in tables.items():
        doc = table_doc(name, ctx.get("group_by")) or {"title": name, "description": ""}
        unit = ctx.get("table_units", {}).get(name)
        trows.append([f"`tables/{name}.tsv`", doc["title"] + (f" ({unit})" if unit else ""), f"{len(df):,}",
                      doc["description"]])
    blocks.append(("table", ["File", "Content", "Rows", "Description"], trows))
    blocks.append(("p", "Every column is described in the README sheet of `results_summary.xlsx` and in "
                        "`LLM_GUIDE.md`."))
    if figs.records:
        blocks += [("h", 2, "Figures"),
                   ("p", "Each figure is in `figures/` as PNG, PDF and SVG."),
                   ("table", ["File", "Figure", "How to read it"],
                    [[f"`figures/{r['key']}.*`", r["title"], r["caption"]] for r in figs.records])]
    blocks += [("h", 2, "Units and missing values"), ("table", ["Term", "Meaning"], [list(u) for u in UNITS])]
    if s.get("arg_methods") or ctx.get("sarg"):
        tier_rows = []
        for tier, cats in tier_categories(ctx.get("arg_priority", "")).items():
            label, meaning = TIER_MEANING[tier]
            tier_rows.append([f"Tier {tier}: {label}", meaning, "; ".join(cats)])
        blocks += [("h", 2, "ARG reporting tiers"),
                   ("p", "CARD-based calls (K-MARVEL, RGI) get a tier from an editable rule table "
                         f"(`{ctx.get('arg_priority_name', 'arg_priority_tiers.tsv')}`): rules match CARD gene names or "
                         "families, top to bottom, first match wins. ARGs-OAP (SARG) has no tiers."),
                   ("table", ["Tier", "Meaning", "Categories (rule notes)"], tier_rows)]
    methods = ctx.get("method_notes", [])
    if methods:
        blocks += [("h", 2, "Methods"), ("ul", methods)]
    blocks += [("h", 2, "Caveats"), ("ul", PITFALLS)]
    if ctx.get("versions_table"):
        blocks += [("h", 2, "Software versions"),
                   ("table", ["Process", "Tool", "Version"], ctx["versions_table"])]
    blocks += [("h", 2, "Inputs"), ("code", "\n".join(ctx["inputs"]))]
    return blocks


def llm_guide(ctx):
    """LLM_GUIDE.md: this summary folder (tables, figures, methods), for a language model or a script. The
    folder is what gets shared, so every path is relative to it and no local path is written."""
    s = ctx["summary"]
    tables = ctx["tables"]
    meta_cols, samples = set(ctx["meta_cols"]), set(ctx["order"])
    out = ["# BugBuster results summary: guide for a language model", "",
           "This file describes the BugBuster metagenomics summary in this folder (one run, or several runs merged "
           "as batches) so that a model or a script can read its tables and figures without guessing. Every path "
           "below is relative to this folder.", "",
           "## Run facts", "", "```yaml",
           f"generated: {date.today().isoformat()}",
           f"n_samples: {s['n_samples']}",
           f"n_batches: {s['n_batches']}",
           "batches:"]
    for b in ctx["batches"]:
        out.append(f"  - label: {b['label']}")
        out.append(f"    samples: [{', '.join(sorted(b['samples'], key=str))}]")
    out += [f"min_clean_pairs: {ctx['min_pairs']}",
            f"low_depth_samples: [{', '.join(s['low_depth_samples'])}]",
            f"arg_primary_method: {s.get('arg_primary', 'none')}",
            f"arg_methods_card: [{', '.join(s.get('arg_methods', []))}]",
            f"rgi_min_coverage_percent: {ctx['rgi_min_coverage']:g}",
            f"metadata_columns: [{', '.join(ctx['meta_cols'])}]",
            f"group_by: {ctx.get('group_by') or 'none'}",
            "```", "",
            "## How to use this folder", "",
            "Contents: `tables/` (one TSV per table), `figures/` (each figure as PNG, PDF and SVG), "
            "`results_summary.xlsx` (the same tables, one sheet each, plus a README sheet), `results_summary.html` "
            "(the report with figures and methods), `README.md` / `README.docx` (the same content for people) and "
            "this guide.", "",
            "1. Start from `tables/samples_qc.tsv` (one row per sample): depth, low-depth flag "
            "(`low_depth`) and per-method counts. Treat low-depth samples as uninformative for absence.",
            "2. Long tables (`*_long.tsv`) hold one row per sample and feature with every annotation column; wide "
            "tables hold one column per sample. Join on `sample` (summary tables) and, for metadata, on the columns "
            "listed in run facts.",
            "3. Wide-table cells are numbers or the codes NA_low_depth, NA_not_run, NA_no_ags (see Units). Convert "
            "them to missing, not to 0, before statistics.",
            "4. The per-sample tool outputs (raw hits, tool reports) are not in this folder; they stay in the "
            "pipeline's results directory of each batch.", "",
            "## Units and codes", ""]
    out += [f"- **{t}**: {m}" for t, m in UNITS]
    out += ["", "## Pitfalls", ""] + [f"- {p}" for p in PITFALLS]
    labels, rows = stage_rows(ctx["stage_status"])
    out += ["", "## Pipeline stages", "",
            "Status per batch: `ran`, `not run`, or `ran, no table found` (the stage ran but its table was not "
            "found). A stage that did not run has no tables or figures here.", ""]
    for i, st in enumerate(STAGES):
        status = "; ".join(f"{b}: {ctx['stage_status'][b][i][1]}" for b in labels)
        out += [f"### {st['title']}", "", f"- Status: {status}", f"- Turned on by: {st['params']}",
                f"- In this folder: {st['summary']}", ""]
    methods = ctx.get("method_notes", [])
    if methods:
        out += ["## Methods", ""] + [f"- {m}" for m in methods] + [""]
    out += ["## Tables (tables/)", "",
            "Each table is also a sheet of `results_summary.xlsx`. Columns are listed in file order.", ""]
    for name, (df, index) in tables.items():
        doc = table_doc(name, ctx.get("group_by")) or {"title": name, "description": ""}
        unit = ctx.get("table_units", {}).get(name)
        out += [f"### {name}.tsv", "", f"{doc['title']}{f' ({unit})' if unit else ''}: {doc['description']} "
                f"{len(df):,} rows.", "", "| Column | Meaning |", "|---|---|"]
        cols = (list(df.index.names) if index else []) + list(df.columns)
        sample_cols = [c for c in cols if c in samples]
        seen_sample = False
        for c in cols:
            if c in samples:
                if seen_sample:
                    continue
                seen_sample = True
                label = f"{sample_cols[0]} ... {sample_cols[-1]}" if len(sample_cols) > 1 else sample_cols[0]
                out.append(f"| {_md_cell(label)} (one per sample) | "
                           f"{_md_cell(column_meaning(name, c, meta_cols, samples, ctx.get('group_by')) or '')} |")
                continue
            meaning = column_meaning(name, str(c), meta_cols, samples, ctx.get("group_by")) or ""
            out.append(f"| `{_md_cell(c)}` | {_md_cell(meaning)} |")
        out.append("")
    if ctx["figures"].records:
        out += ["## Figures (figures/)", ""]
        out += [f"- `figures/{r['key']}.png` (.pdf, .svg): {r['title']}. {r['caption']}"
                for r in ctx["figures"].records]
        out.append("")
    out += ["## Recipes", "",
            "- Which samples are reliable? `tables/samples_qc.tsv`: `low_depth == False`; `depth_class` for depth strata.",
            "- Most abundant species of a sample: `tables/taxonomy_long.tsv`, filter `sample`, sort `relative_abundance`.",
            "- Clinically relevant ARGs: `tables/arg_long.tsv` with `tier == 1` and `gene_representative == True`; check "
            "`tables/arg_method_comparison.tsv` (`detected_by`) for agreement between K-MARVEL and RGI, and "
            "`rgi_host_check` for the likely host.",
            "- ARG load per sample in one unit: `tables/samples_qc.tsv` `acquired_rpkm` (CARD methods) or "
            "`sarg_copies_per_cell` (ARGs-OAP, SARG); never add the two.",
            "- A functional term across samples: `tables/function_read_<ontology>.tsv` (CPGE when available), or "
            "`tables/function_read_long.tsv` filtered by `ontology` and `accession`; the contig-level profile is "
            "`tables/function_contig_<ontology>.tsv` / `function_contig_long.tsv` (never mix the two).",
            "- Good-quality genomes: `tables/mags.tsv` with `quality_tier == 'high'` (or 'medium'), lineage in "
            "`phylum` ... `species`; their depth per sample in `tables/mag_abundance.tsv` (`pct_sample_bases`).",
            "- Which organisms carry an ARG (linkage, assembled contigs): `tables/arg_contig_long.tsv` (`host_taxon`, "
            "`bin`, `mobile_class`), summarised in `tables/arg_contig_hosts.tsv`; ARGs shared between MAGs in "
            "`tables/arg_cluster_sharing.tsv`.",
            "- Software versions: `README.md`, section Software versions (per batch when several were merged).", ""]
    return "\n".join(out)


def write_docs(outdir, ctx):
    blocks = readme_blocks(ctx)
    (outdir / "README.md").write_text(render_markdown(blocks))
    render_docx(blocks, outdir / "README.docx", ctx["title"])
    (outdir / "LLM_GUIDE.md").write_text(llm_guide(ctx))
