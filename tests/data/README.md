# Test fixtures

Every fixture the tests read is committed here. No test downloads data. This table says how each folder was made:

- **generated**: a `tests/bin/make_*` script rewrites it. Re-run the script after changing it, and check `git diff`. Most scripts are byte-stable; the exceptions are noted.
- **hand-written**: small files written for the tests named. Edit them directly.
- **captured**: copied from a real tool run. They cannot be regenerated without that tool's database.

Nothing here is real sample data. The CARD-layout and DeepARG-layout fixtures are synthetic: made-up names and sequences in the tools' file layouts (roadmap decision RV8: the public repository carries no CARD data).

| Folder | How | Source / generator | Used by |
| --- | --- | --- | --- |
| `arg_clustering/` | generated | `make_arg_cluster_fixtures.py`: two samples' contig-to-bin tables, Pyrodigal-style proteins and DeepARG hits (identical, 93 % and unrelated proteins) | `modules/arg_fasta_formatter`, `modules/clustering` |
| `args_oap/` | captured, edited | rows of a real ARGs-OAP 4.0.0 stage two output, renamed (read names synthetic), metadata written for them (see the `make_summary_fixtures.sh` header) | `test_report_scripts.sh`, `make_summary_fixtures.sh` |
| `bam/` | generated | `make_bam_fixtures.sh` | `modules/featurecounts_genes` |
| `bins/real_bins/`, `marker_only_bins/` | generated | `make_bin_fixtures.py` (one 25 kb seeded contig; `marker_only_bins/` is no longer read by a test) | `modules/bin_depth`, `pyrodigal`, `bwamem3_*`, `eggnog_mapper_search` |
| `bins/` (other folders, `overlap_depth.txt`) | hand-written | bin sets and tables for the tests that read them (overlapping bins, a single-bin batch, a bin without classification, a no-bin batch, short MetaBAT2 contigs) | `modules/checkm2_batch`, `gtdb_tk_batch`, `metabat2`, `bin_depth`, `gunc_run` |
| `checkv/` | hand-written | an empty gzip, a truncated one, and a `db/README.txt` with only the version line | `modules/checkv_endtoend` |
| `contigs/` | generated | `make_contig_fixtures.py` | `modules/bbmap`, `metabat2` |
| `decontamination/` | captured | phage lambda (NC_001416.1, the mock host) and phiX174 genomes (`.gz` and no-final-newline variants derived from them) | the QC and decontamination tests, `-profile test`, and as genomes in the sylph, metabuli, depth and BOWTIE2_SAMTOOLS fixtures |
| `derep/` | generated | `make_derep_fixtures.py` | `modules/galah_cluster`, `coverm_genome` |
| `functional/` | generated | `make_functional_fixtures.sh` | `modules/aggregate_functions`, `run_dbcan`, `test_aggregate_functions.sh`, `test_db_reformat_scripts.sh` |
| `gff/` | generated | `make_gff_fixtures.sh` | `modules/featurecounts_genes`, `test_gff2saf.sh` |
| `gunc/` | generated | `make_gunc_fixtures.sh` | `modules/gunc_run` |
| `humann/` | generated | `make_humann_fixtures.sh` (database, reads) and `make_humann_fixtures.py` (tables) | `modules/humann`, `test_humann_scripts.sh`, `test_db_reformat_scripts.sh` |
| `kmarvel/db_synthetic/`, `kmarvel_sim_R*.fastq.gz` | generated | `make_kmarvel_fixtures.py`: made-up proteins in CARD 4.0.1's layout, reads from one of them | `modules/kmarvel` |
| `kmarvel/sample*_kmarvel.*` | hand-written | K-MARVEL tables and logs in its output layout (sampleA: two ARGs, sampleB: none) | `modules/kmarvel_report`, `results_summary`, `test_report_scripts.sh`, `make_summary_fixtures.sh` |
| `kraken2/stub_db/` | hand-written | empty placeholders of a Kraken2 database directory | `main.nf.test` (kraken2 stub run) |
| `metabuli/` | generated | `make_metabuli_fixtures.sh` | `modules/metabuli_classify` |
| `multisample/` | hand-written | three assemblies (wrapped MEGAHIT headers, a tab in a header, none) and a group depth table | `modules/multisample_concat`, `multisample_split_depth` |
| `nonpareil/` | generated | `make_nonpareil_fixtures.sh` | `modules/nonpareil`, `test_results_summary.sh` |
| `reads/test_sample{1,2}_*` | generated | `make_test_read_fixtures.sh`: the first 500 pairs of the nf-core reads `assets/test_samplesheet.csv` points to (the only generator that needs the network) | every pipeline-level nf-test, the CI profile runs (`--input`), `test_resource_limits.sh` |
| `reads/singleton_test_*` | generated | `make_singleton_fixtures.py`: 10 pairs (4 from lambda) and 10 singletons (6 from lambda), and plain copies of the mates | the QC tests, `modules/*_decontaminate`, `bowtie2_samtools`, `calculate_depth`, `sylph_profile`, `count_reads` |
| `reports/` | hand-written | report-stage inputs in each tool's layout | `modules/contig_filter_summary`, `reads_report`, `test_report_scripts.sh` |
| `rgi/wildcard_mini/` | generated | `make_rgi_wildcard_fixture.sh` (synthetic CARD + WildCARD) | `modules/rgi_load_wildcard` |
| `rgi/` (other files) | hand-written | RGI bwt tables, an empty k-mer JSON, and a `localDB/` without k-mers | `modules/rgi_kmer`, `rgi_bwt`, `rgi_report` |
| `sourmash/two_genome_*`, `sm_test_R*.fastq.gz` | generated | `make_sourmash_fixtures.sh`: two seeded genomes (120 / 100 kb), their lineages, and reads from one. The zip holds sourmash's timestamps, so a rebuild differs in bytes but not in content | `modules/sourmash` |
| `sourmash/stub_*` | hand-written | placeholders for the `-stub` sourmash runs | `main.nf.test` |
| `stub_paths/` | generated by the tests | every path a stub run with the tests' settings publishes. Regenerate with `BUGBUSTER_UPDATE_STUB_PATHS=1 nf-test test tests/nonpareil.nf.test tests/coassembly_groups.nf.test` (`TestData.publishedPathsDiff`), review the diff, then run the other list tests (`mag_virus_quality`, `dereplicate`, `multisample`) | those five test files |
| `summary/` | generated | `make_summary_fixtures.sh`, `make_grouped_summary_tree.py` | `test_results_summary.sh`, `modules/results_summary` |
| `superfocus/` | generated | `make_superfocus_fixtures.sh` (database, reads), `make_superfocus_fixtures.py` (tables) | `modules/superfocus`, `test_superfocus_scripts.sh` |
| `sylph/` | generated | `make_sylph_fixtures.sh` (lambda + phiX database, simulated reads) | `modules/sylph_profile`, `sylph_taxprof`, `calculate_depth` |
| `taxonomy/*sylph_gtdb.sylphmpa` | captured | sylph-tax 1.9.2 on a real sylph 1.0.0 profile against the GTDB r232 database (simulated E. coli + B. subtilis reads). Cannot be regenerated without that 20 GB database | `test_taxonomy_scripts.sh`, `make_summary_fixtures.sh` |
| `taxonomy/` (other files) | hand-written | Kraken2 / Bracken reports, sourmash gather tables and per-sample reports in the tools' layouts | `modules/taxonomy_report`, `taxonomy_phyloseq`, `phyloseq_converter`, `test_taxonomy_scripts.sh` |
| `versions/` | hand-written | a `versions.yml` as local modules write it | `subworkflows/software_versions` |
| `woltka/` | generated | `make_woltka_fixtures.sh` (database, alignments), `make_woltka_fixtures.py` (tables) | `modules/woltka_*`, `test_aggregate_read_functions.sh` |
