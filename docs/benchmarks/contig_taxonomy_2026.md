# Contig taxonomy and bin-level ARGs

Benchmark run 2026-10-27 to 2026-10-28. It asks two questions about the contig-level ARG branch (`--contig_tax_and_arg`) and the bin-level ARG branch (`--arg_bin_clustering`):

- **(L2) Contig taxonomy:** BLASTN of every contig against nt, then BlobTools (today's NT_BLASTN → BLOBTOOLS), against Kraken2 on the contigs. Runtime, peak memory, database size, contigs classified per rank, and agreement on the ARG-bearing contigs that `Contig_tax_and_arg_prediction.tsv` reports. MMseqs2 taxonomy and CAT need databases that are not on the server; their sizes are listed, and they were not run.
- **(M8) Bin-level ARGs:** DEEPARG_BINS (Prodigal per bin, then DeepARG per bin) against the contig-level DeepARG hits (DEEPARG_CONTIGS on the shared Pyrodigal proteins) mapped to bins through bin membership.

No pipeline code changed. Every run used dev's own modules (Setup).

**Decision (owner, 2026-10-28):**

- **L2:** Kraken2 on the contigs with `k2_standard`, replacing BLASTN against nt. A GTDB database for Kraken2 was considered and not taken: the local GTDB r220 build needs about 500 GB in memory and has no viruses, plasmid-only sequences or eukaryotes. sylph does not classify single contigs (it estimates how much of each reference genome is in a sample).
- **M8:** derive the bin-level ARGs from the contig-level DeepARG run; drop PRODIGAL_BINS and DEEPARG_BINS.
- **Taxdump bug:** fixed together with these changes, not as a separate fix on the development branch.

**Outcome (recommendations):**

- **L2: Kraken2 on the contigs instead of BLASTN against nt, keeping BlobTools for the tables and plots.** BLASTN is the most expensive step of the branch by far: 29,093 CPU-s (8.1 CPU-h) for the smallest assembly (3.7 Mbp) and 201,739 CPU-s (56 CPU-h, 4.7 h at 16 threads) for a 36.6 Mbp one: every ~5 Mbp of contigs adds another pass over the 434 GB of nt, so the two larger assemblies would take about 140 and 180 CPU-h. Kraken2 took 670–1,400 CPU-s per assembly, database load included (less when the database is already in the page cache): 668 CPU-s on the 36.6 Mbp one, 300 times less than BLASTN. Kraken2 classified each assembly in 1–6 s once its database was loaded (1–5 min per task with the load). Where both tools assign a rank to an ARG-bearing contig, they agree on 99 % of superkingdoms, 97 % of phyla and orders, 95 % of families and genera and 89 % of species. Kraken2 assigns fewer contigs (78–84 % of the ARG-bearing contigs of samples A, B and C against 98 %; of all contigs, 80 % against 91 % in sample C and 45 % against 67 % in sample D), because it stops at the lowest common ancestor where BlobTools takes the best hit's species. Section (L2) lists the alternatives, two of which keep BLAST for the ARG-bearing contigs.
- **M8: derive the bin-level ARGs from the contig-level DeepARG run, and drop PRODIGAL_BINS + DEEPARG_BINS.** The two routes agree on 228 of 247 hits (same ARG and class every time). 17 of the 19 that differ are on genes the two routes called identically. They differ because **DeepARG scales its features over the whole input file** (min-max), so one protein's probability depends on what else is in the file. A bin with few candidate genes gets unreliable probabilities: a class A β-lactamase scores 0.227 (dropped) when its bin is run alone and 0.994 when all bins of the sample are run together. A single DeepARG run on all bins' proteins reproduced the contig route's in-bin calls exactly (91 of 91, probabilities within 0.014). The bin route also costs about 1 CPU-hour and 4–5 GB per task for the three samples with bins, on top of the contig route.
- **Fix whatever is decided (pipeline bug, live on dev today):** the default `--taxdump_files ncbi` downloads the current NCBI taxdump. NCBI replaced the `superkingdom` rank with `domain` in 2025, and with that taxdump BlobTools 1.1.1 writes `undef` as every contig's superkingdom. With the default settings, the `superkingdom` column of `Contig_tax_and_arg_prediction.tsv` is `undef` on every row that has a hit. The other ranks are unaffected.
- **Fix if Kraken2 is used with the current standard database:** the 2026-06-26 k2_standard needs 103 GB of memory. KRAKEN2's `process_high` first attempt (72 GB) fails with it, and only the retry at 144 GB fits. The read-level KRAKEN2 has the same problem when `--custom_kraken_db` points at this database.

## Summary

| Question | Measured | Recommendation |
|---|---|---|
| Contig taxonomy cost (L2) | BLASTN: 29,093 CPU-s for 3.7 Mbp, 201,739 CPU-s (56 CPU-h) for 36.6 Mbp, ~5 Mbp of contigs per pass over nt; 434 GB database. Kraken2: 221–1,393 CPU-s per task (database load included), classification itself 1–6 s; 104 GB database, 103 GB in memory | Kraken2 |
| Contig taxonomy results (L2), ARG-bearing contigs of all four samples (317) | classified at genus: BlobTools 278, Kraken2 190. Where both assign the rank: superkingdom 254/256, phylum 233/239, order 204/210, family 191/201, genus 178/187, species 111/125 agree (same 2026 taxonomy on both sides) | Kraken2, with the trade-off stated |
| MMseqs2 taxonomy, CAT (L2) | not run: no local database. GTDB R232 proteins 132 GB, UniRef90 32 GB, nr 200 GB (compressed downloads); CAT's prebuilt archives' sizes could not be read | not pursued unless the owner wants the download |
| Taxdump (bug) | current NCBI taxdump → `superkingdom` = `undef` on every contig with a hit | fix in 9b |
| Kraken2 memory | k2_standard 2026-06-26: 103 GB anon; `process_high` first attempt 72 GB | fix in 9b if Kraken2 is used |
| Bin-level ARGs (M8) | 228 of 247 hits on both routes, all with the same ARG and class; 17 of the 19 others are DeepARG's input-dependent scaling, 2 differ in the gene call; bin route ≈ 1.0 CPU-h and 4.2–5.2 GB per task extra | derive from the contig route |

The rest of this document gives the measurements behind each row.

## Setup

**Contigs.** The four asm assemblies of the test fixture (real reads, subsampled; Phase 7a's fixed MEGAHIT contigs ≥ 1 kb, since MEGAHIT does not give the same contigs run to run), with the pipeline's bwa-mem3 BAMs, and Phase 7b's Binette bins of the same contigs:

| Sample | Contigs ≥ 1 kb | Mbp | Binette bins (contigs in bins) | ARG-bearing contigs (DeepARG hits) |
|---|---|---|---|---|
| A | 25,366 | 117.8 | 22 (4,476) | 142 (190) |
| B | 16,094 | 93.8 | 17 (2,482) | 132 (166) |
| C | 10,542 | 36.6 | 6 (2,390) | 37 (41) |
| D | 2,078 | 3.7 | none | 6 (6) |

**Harness.** A scratch Nextflow workflow ran dev `1a987d4`'s PYRODIGAL, DEEPARG_CONTIGS, PRODIGAL_BINS, DEEPARG_BINS, SAMTOOLS_INDEX, BLOBTOOLS and ARG_CONTIG_LEVEL_REPORT modules unchanged, with their labels' first-attempt resources (8 CPUs, 36 GB for `process_medium`). BLASTN and Kraken2 ran outside Nextflow in their containers, with the module's command and the CPUs and memory a task gets, so the cost could be measured from the container's cgroup.

**BLASTN** (BLAST 2.17.0, the NT_BLASTN image and command: `blastn -task megablast -outfmt '6 qseqid staxids bitscore std' -max_target_seqs 1 -max_hsps 1 -evalue 1e-25`, 16 threads, 72 GB, default threading only; `-mt_mode 1` was dropped after Phase 5). The database is nt at `/mds_data/dbs/blast_db_nt`: 434 GB, volumes dated May–June 2024. Phase 5 had already measured sample D's whole assembly. The owner chose to add sample C's whole assembly and only the ARG-bearing contigs of samples A and B (their whole assemblies would have taken about 12–15 h each on 16 threads):

| Run | Contigs | Mbp |
|---|---|---|
| D, whole assembly (Phase 5) | 2,078 | 3.7 |
| C, whole assembly | 10,542 | 36.6 |
| B, ARG-bearing contigs only | 132 | 6.56 |
| A, ARG-bearing contigs only | 142 | 8.67 |

**BlobTools** (1.1.1, the module's image) ran on each BLAST output through the BLOBTOOLS module, twice: with an NCBI taxdump of March 2024 (it matches the nt release) and with the taxdump of the Kraken2 database (June 2026). Section (L2) explains why both.

**Kraken2** 2.17.2 (the pipeline's KRAKEN2 image, `community.wave.seqera.io/library/kraken2_coreutils_pigz:658cf70caab8666d`) with `--confidence 0.1` (the pipeline default) on each whole assembly, 16 threads, against `/mds_data/dbs/kraken2/k2_standard`, the standard database of 2026-06-26 (the owner's choice). **One Kraken2 task at a time:** the four samples ran one after another, and no other Kraken2 process ran on the host. Kraken2 2.17.2 loads its database with several threads, and three tasks loading at once were about twice as slow in Phase 6. BLASTN and harness tasks did run alongside it.

**Per-contig lineages.** Kraken2's per-sequence output was expanded to a lineage with the database's own `nodes.dmp` and `names.dmp`, at the six ranks BlobTools reports (superkingdom, phylum, order, family, genus, species; BlobTools has no class). Kraken2's `domain` rank stands in for superkingdom, and virus taxa are put under "Viruses". A contig counts as classified at a rank when the tool gives it a name there; BlobTools' `no-hit`, `undef` and `*-undef` count as unclassified.

**DeepARG** 1.0.4 (the module's image), with the pipeline's default arguments (`--min-prob 0.8 --arg-alignment-identity 50 --arg-alignment-evalue 1e-10 --arg-alignment-overlap 0.8 --arg-num-alignments-per-entry 1000 --model-version v2`) and the threads bounded as in PR #41.

## Machine and runs

- 4 × Xeon E7-4820 v4 (40 cores, 80 threads), 2 TB RAM; Docker; Nextflow 26.04.4.
- Other users kept the 1-minute load at 62–139 on 80 threads. The BLAST runs overlapped each other and the harness. Wall times are inflated.
- **Time and memory come from the container's cgroup**, sampled every 2 s (GB here are 2^30 bytes): CPU-seconds, peak cgroup memory (it includes page cache) and peak anonymous memory ("anon", the process memory). The harness tasks' figures come from Nextflow's trace (realtime, %CPU, peak RSS).

## (L2) Contig taxonomy

### Cost

**BLASTN against nt:**

| Run | Mbp | Wall | CPU-s | Mean cores | Peak anon | Peak cgroup |
|---|---|---|---|---|---|---|
| D, whole (Phase 5) | 3.7 | 93.3 min | 29,093 | 5.2 | 0.19 GB | 72 GB |
| D, whole, repeat (Phase 5) | 3.7 | 72.2 min | 30,442 | 7.0 | | 72 GB |
| C, whole | 36.6 | 281.1 min | 201,739 | 12.0 | 0.44 GB | 72 GB |
| B, ARG-bearing contigs | 6.56 | 66.2 min | 46,234 | 11.6 | 0.25 GB | 72 GB |
| A, ARG-bearing contigs | 8.67 | 84.0 min | 53,044 | 10.5 | 0.28 GB | 72 GB |

- **Cost grows with the amount of contig sequence, in steps.** Megablast takes the query in batches of about 5 Mbp, and each batch is a full pass over nt. One pass costs roughly 23,000–30,000 CPU-s here: one batch for D, two for B's and A's ARG-bearing contigs, and eight for C's whole assembly (201,739 CPU-s, 25,200 per batch), as the batch size predicts. At C's 5,500 CPU-s per Mbp, B's whole assembly (93.8 Mbp) would take about 520,000 CPU-s (140 CPU-h) and A's (117.8 Mbp) about 650,000 (180 CPU-h), roughly 12–15 h each at 16 threads. These were not run.
- **Memory is page cache.** BLASTN's own memory is 0.2–0.3 GB. The 72 GB cgroup fills with nt pages, and with 434 GB of nt nothing useful stays cached between runs (the two D runs cost the same).
- **BlobTools** took 40–86 s and 1.4–1.7 GB peak RSS per sample, mostly reading the taxdump and the BAM.

**Kraken2 on contigs** (one task at a time):

| Sample | Mbp | Wall | CPU-s | Peak anon | Peak cgroup | Classification (Kraken2's own timer) | Contigs classified |
|---|---|---|---|---|---|---|---|
| D | 3.7 | 300.8 s | 1,393 | 103.0 GB | 122.8 GB | 1.13 s | 944 / 2,078 (45.4 %) |
| C | 36.6 | 150.8 s | 668 | 103.3 GB | 131.1 GB | 1.07 s | 8,383 / 10,542 (79.5 %) |
| B | 93.8 | 66.8 s | 330 | 103.5 GB | 110.3 GB | 3.64 s | 11,230 / 16,094 (69.8 %) |
| A | 117.8 | 57.8 s | 221 | 103.6 GB | 104.3 GB | 6.17 s | 17,531 / 25,366 (69.1 %) |

- **The run is the database load.** Classification took 1–6 s. D ran first, with the database cold on disk; the later runs read it partly from the page cache, which is why the wall time falls as the samples get larger. A task in a real run loads the database once per sample, so its cost is the load, about 1–5 min.
- **The database needs 103 GB of memory.** It loads `hash.k2d` (103 GB) into memory. These runs were given 144 GB. KRAKEN2 and a contig-level Kraken2 at `process_high` (72 GB first attempt, then 144 and 216 GB) would fail the first attempt with this database. Kraken2's `--memory-mapping` would read the database from the page cache instead of loading it; it was not tested here.

**Database sizes:**

| Database | Size | Notes |
|---|---|---|
| nt (BLAST, May–Jun 2024) | 434 GB | NT_BLASTN's registry downloads the current nt |
| NCBI taxdump (`nodes.dmp`, `names.dmp`) | 0.4 GB (2024), 0.5 GB (2026) | |
| Kraken2 k2_standard 2026-06-26 | 104 GB (`hash.k2d` 103 GB; the download is an 80 GB `.tar.gz`) | includes its own taxdump |
| Kraken2 standard-8 (the registry's default for reads) | 8 GB | not tested on contigs |

### Contigs classified

**All contigs, Kraken2** (contigs with a name at each rank):

| Sample | Contigs | Superkingdom | Phylum | Order | Family | Genus | Species |
|---|---|---|---|---|---|---|---|
| A | 25,366 | 17,525 | 17,369 | 16,935 | 16,554 | 15,608 | 12,521 |
| B | 16,094 | 11,225 | 11,102 | 10,741 | 10,323 | 9,381 | 7,370 |
| C | 10,542 | 8,380 | 8,331 | 8,148 | 7,951 | 7,840 | 6,725 |
| D | 2,078 | 943 | 918 | 864 | 776 | 741 | 677 |

**All contigs of samples C and D, both tools** (BlobTools / Kraken2 / both / agree; superkingdom from BlobTools with the March 2024 taxdump, the other ranks with the 2026 taxdump, as in the next section):

| Rank | C (10,542 contigs) | D (2,078 contigs) |
|---|---|---|
| Superkingdom (2024 taxdump) | 9,533 / 8,380 / 8,323 / 8,180 (98 %) | 1,394 / 943 / 937 / 909 (97 %) |
| Phylum | 9,363 / 8,331 / 8,201 / 8,107 (99 %) | 1,364 / 918 / 898 / 878 (98 %) |
| Order | 9,068 / 8,148 / 7,936 / 7,877 (99 %) | 1,287 / 864 / 831 / 818 (98 %) |
| Family | 9,008 / 7,951 / 7,727 / 7,635 (99 %) | 1,283 / 776 / 747 / 725 (97 %) |
| Genus | 8,938 / 7,840 / 7,601 / 7,486 (98 %) | 1,239 / 741 / 702 / 668 (95 %) |
| Species | 9,553 / 6,725 / 6,680 / 6,313 (95 %) | 1,394 / 677 / 673 / 589 (88 %) |

- **BLAST finds a hit for more contigs:** 9,557 of C's 10,542 (91 %) against Kraken2's 8,383 (80 %), and 1,396 of D's 2,078 (67 %) against 944 (45 %). Few contigs are classified by Kraken2 only (57 in C, 6 in D, at superkingdom).
- **"Classified" means different things.** BlobTools gives each contig the full lineage of its taxid, chosen from the contig's own hits (here one hit per contig, since NT_BLASTN keeps `-max_target_seqs 1`). So every contig with a hit gets a species, whatever the identity. Kraken2 at `--confidence 0.1` stops at the lowest common ancestor of the k-mers it matched, and gives a species to two thirds to four fifths of the contigs it classifies (66–80 %), but to only about half of the ARG-bearing ones. On the ARG-bearing contigs of A and B, the best BLAST hit had a median identity of 99.6 %, but for a tenth of them it was below 82–85 %, and for another tenth (in part the same contigs) it covered less than 15–19 % of the contig. For those, BlobTools' species is a nearest relative, not a species call.

### ARG-bearing contigs (`Contig_tax_and_arg_prediction.tsv`)

The report lists each contig with a DeepARG hit, with BlobTools' lineage. Agreement is counted on those contigs, at each rank where both tools assign a name.

**The same taxonomy on both sides.** Against the March 2024 taxdump, many "disagreements" were renames between the 2024 and 2026 taxonomies: BlobTools' *Thomasclavelia* / Coprobacillaceae / Erysipelotrichales against Kraken2's *Clostridium* / Clostridiaceae / Eubacteriales (8 contigs in A, 6 in B), *Lachnoclostridium* against *Otoolea* or *Enterocloster*. Rerun with the Kraken2 database's 2026 taxdump, BlobTools agreed with Kraken2 on 91 instead of 77 orders in A, 89 instead of 81 in B and 22 instead of 16 in C. The table below uses the 2026 taxdump for every rank except superkingdom, which BlobTools cannot assign with it (next section), so that row uses the 2024 run.

| Rank | A (142) BlobTools / Kraken2 / both / agree | B (132) | C (37) | D (6) | All four: agree / both |
|---|---|---|---|---|---|
| Superkingdom (2024 taxdump) | 139 / 111 / 111 / 109 | 129 / 110 / 110 / 110 | 37 / 31 / 31 / 31 | 6 / 4 / 4 / 4 | 254 / 256 (99 %) |
| Phylum | 138 / 103 / 103 / 99 | 126 / 105 / 104 / 103 | 37 / 28 / 28 / 27 | 6 / 4 / 4 / 4 | 233 / 239 (97 %) |
| Order | 130 / 97 / 96 / 91 | 126 / 90 / 90 / 89 | 37 / 22 / 22 / 22 | 6 / 2 / 2 / 2 | 204 / 210 (97 %) |
| Family | 125 / 95 / 94 / 86 | 122 / 87 / 87 / 85 | 35 / 20 / 19 / 19 | 6 / 1 / 1 / 1 | 191 / 201 (95 %) |
| Genus | 117 / 87 / 85 / 80 | 123 / 83 / 83 / 79 | 34 / 19 / 18 / 18 | 4 / 1 / 1 / 1 | 178 / 187 (95 %) |
| Species | 139 / 57 / 57 / 49 | 129 / 51 / 51 / 46 | 37 / 16 / 16 / 15 | 6 / 1 / 1 / 1 | 111 / 125 (89 %) |

- **BLAST hit 311 of the 317 ARG-bearing contigs; Kraken2 classified 256.** At genus, BlobTools named 278 and Kraken2 190; at species, 311 and 125.
- **Where both name a rank, they agree** on 97 % of phyla, 95 % of genera and 89 % of species. Most of the remaining genus and species disagreements are between close relatives (*Phocaeicola vulgatus* / *P. dorei*, *Bacteroides* / *Phocaeicola*, *Lachnoclostridium* / *Enterocloster*). The phylum disagreements are few but far apart: one contig in each of A, B and C that BLAST placed in *Campylobacter* and Kraken2 in *Enterococcus*, and two phage hits (Uroviricota) in A that Kraken2 placed in Bacillota. ARGs often sit on mobile elements shared across phyla, so neither tool's answer is wrong by construction for these.

### The taxdump bug (BlobTools 1.1.1)

NCBI replaced the rank `superkingdom` with `domain` (and `realm` for viruses) in its 2025 taxonomy. BlobTools 1.1.1 looks for `superkingdom` by name. With the 2026 taxdump, every contig with a hit got `undef` as its superkingdom (the report built here, over the BLAST runs above, 12,980 rows: 11,307 `undef` and 1,673 `no-hit`; with the 2024 taxdump: 10,257 Bacteria, 594 Viruses, 331 Eukaryota, 99 Archaea, 26 `undef`, 1,673 `no-hit`). Phylum and the lower ranks were unaffected (D's largest phylum: 645 contigs with the 2024 taxdump, 644 with the 2026 one).

The pipeline's default `--taxdump_files ncbi` downloads `taxdump.tar.gz` from NCBI at run time (`config/databases.config`: "unversioned: latest state at download time"), so **every new run with the default settings has this bug**. The blob plots color by phylum and are not affected. Two fixes: pin a taxdump from before the change, or rename the rank in `nodes.dmp` inside BLOBTOOLS before `blobtools create` (`domain` → `superkingdom`, and the same for virus `realm`), which keeps the current taxonomy.

### Alternatives (not measured)

- **Keep BLASTN, only for the ARG-bearing contigs.** The report needs the taxonomy of these contigs only; the blob plots would lose their taxonomy (or take it from Kraken2). ARG-bearing contigs are long, so this is still 6.6–8.7 Mbp for A and B: 46,000–53,000 CPU-s each, against about 520,000 and 650,000 CPU-s for their whole assemblies (a 10–12× saving).
- **Kraken2 for every contig, BLASTN for the ARG-bearing ones.** BlobTools' plots and table from Kraken2's taxids; the report's lineage from BLAST. Costs the sum.
- **Kraken2 with a smaller database.** The registry's read-level default (`standard-8`, 8 GB) would load an 8 GB database instead of 103 GB; it classifies fewer sequences than the full standard database. Not tested on contigs.
- **BlobTools from Kraken2.** `blobtools create --hitsfile` takes `contig<TAB>taxid<TAB>score`; Kraken2's per-sequence output gives the contig and taxid (the score would be a constant), so BLOBTOOLS, BLOBPLOT and the report could stay as they are. Not tested.

### MMseqs2 taxonomy and CAT

Neither has a database on the server, and nothing was downloaded. Sizes of the downloads (HTTP headers, 2026-10-27):

| Database | Used by | Download |
|---|---|---|
| GTDB R232 protein reps (latest) | MMseqs2 `databases GTDB`, CAT GTDB build | 131.9 GB |
| GTDB R226 protein reps | | 93.9 GB |
| UniRef90 (2026-09) | MMseqs2 `databases UniRef90` | 32.1 GB |
| UniRef50 | MMseqs2 `databases UniRef50` | 8.8 GB |
| NCBI nr (`nr.gz`) | CAT nr build, MMseqs2 `databases NR` | 200.0 GB |
| `prot.accession2taxid.FULL.gz` | CAT nr build | 26.7 GB |
| CAT prebuilt `20240422_CAT_nr.tar.gz`, `20231120_CAT_gtdb.tar.gz` | CAT | not readable: the server's TLS certificate chain fails verification |

Built databases are larger than the downloads: MMseqs2 builds its own index of each, and a CAT nr build has been reported to reach 350 GB. These are protein-level classifiers (they need gene calls, which PYRODIGAL already makes) that align every protein against the reference set, against Kraken2's k-mer lookup. They were left out.

## (M8) Bin-level ARGs: DEEPARG_BINS against contig-level hits mapped to bins

### The two routes

- **Bin route** (today, `--arg_bin_clustering`): PRODIGAL_BINS calls genes in each bin with Prodigal 2.6.3 `-p single`; DEEPARG_BINS runs DeepARG once per bin.
- **Contig route** (today, `--contig_tax_and_arg`): PYRODIGAL calls genes on the whole assembly with `-p meta`; DEEPARG_CONTIGS runs DeepARG once per assembly. Here its hits were mapped to bins by contig name, and a bin-route hit matched a contig-route hit on the same contig and strand when the two genes overlap by at least half the shorter one.

### Results

| Sample | Bins | Bin route hits | Contig route hits in bins | Both | Same ARG and class | Same gene coordinates | Bin route only | Contig route only | Bins with the same ARG set |
|---|---|---|---|---|---|---|---|---|---|
| A | 22 | 124 | 125 | 123 | 123 | 115 | 1 | 2 | 19 / 22 |
| B | 17 | 102 | 91 | 89 | 89 | 80 | 13 | 2 | 10 / 17 |
| C | 6 | 17 | 16 | 16 | 16 | 16 | 1 | 0 | 6 / 6 |
| Total | 45 | 243 | 232 | 228 | 228 | 211 | 15 | 4 | 35 / 45 |

Each sample had the same distinct ARGs on both routes (A 70, B 38, C 11).

- **Where both routes call a hit, it is the same ARG and class** (228 of 228). 17 of them are on genes with different coordinates on the two routes (overlapping calls of the same gene).
- **Two of the 19 single-route hits are gene-calling differences.** In C, Prodigal called a 372-bp *vanS* gene that Pyrodigal did not call. In A, a class A β-lactamase is called by both with different starts (GTG at 5,006 against ATG at 5,069), and only the Pyrodigal protein passed.
- **The other 17 are on genes with identical coordinates on both routes**, with the same single DIAMOND hit on both (e.g. an *efrB* gene: 51.2 % identity, e-value 6.7e-89, bitscore 278 on both routes; probability 0.554 on the contig route, 0.859 on the bin route). Most bin-route-only hits are *efrA*/*efrB* at 50–53 % identity with probability 0.84–0.92. Two contig-route-only hits had high probabilities: *rosA* (0.999) and a class A β-lactamase (0.994, 79.7 % identity).

### Why: DeepARG scales its features over the input file

DeepARG 1.0.4 (`deepARG.py`, lines 70–72) applies scikit-learn's `MinMaxScaler().fit_transform` to the feature matrix of all proteins that pass its alignment filters in one run. Each protein's features are therefore scaled by the minimum and maximum over the other proteins in the same file, and its probability depends on them. To check, DeepARG was run on the proteins of the bin holding the class A β-lactamase above (sample B, 4,310 proteins):

| Input | Proteins | Passing DeepARG's filters | Probability of the class A β-lactamase |
|---|---|---|---|
| This bin's Prodigal proteins (as DEEPARG_BINS runs it) | 4,310 | 1 | 0.227 (below 0.8: dropped) |
| This bin's contigs, Pyrodigal proteins | 4,350 | 1 | 0.227 |
| All of sample B's bins' Prodigal proteins in one file | 56,167 | 132 | 0.994 |
| Sample B's whole assembly (DEEPARG_CONTIGS) | | | 0.994 |

The gene caller made no difference, and the input set made all of it. Run on all bins' proteins at once, DeepARG reproduced the contig route's hits in bins exactly: 91 of 91 on the same contig with the same best hit, probabilities within 0.014. So the bin route's results are less reliable than the contig route's, not just different: a bin with one or a few candidate genes is scaled against itself. The `--deeparg_min_prob` threshold is also not comparable between runs with different inputs.

### Cost

From the harness trace (8 CPUs per task):

| Task | A | B | C | D |
|---|---|---|---|---|
| PYRODIGAL (shared with the functional branch) | 4 min 34 s, 2.2 GB | 3 min 46 s, 1.8 GB | 1 min 17 s, 0.8 GB | 11 s, 0.1 GB |
| DEEPARG_CONTIGS | 2 min 2 s, ~377 CPU-s, 0.9 GB | 2 min 2 s, ~326 CPU-s, 0.9 GB | 2 min 31 s, ~231 CPU-s, 0.9 GB | 1 min 50 s, ~130 CPU-s, 0.9 GB |
| PRODIGAL_BINS | 2 min 26 s, ~533 CPU-s, 0.8 GB | 2 min 2 s, ~483 CPU-s, 0.8 GB | 37 s, ~141 CPU-s, 0.5 GB | (no bins) |
| DEEPARG_BINS | 4 min 2 s, ~1,234 CPU-s, 5.2 GB | 2 min 56 s, ~832 CPU-s, 4.8 GB | 1 min 19 s, ~385 CPU-s, 4.2 GB | (no bins) |

CPU-s is realtime × %CPU from the trace; memory is peak RSS. The bin route (PRODIGAL_BINS + DEEPARG_BINS) cost about 3,600 CPU-s (1.0 CPU-h) over the three samples with bins, and DEEPARG_BINS needs 4–5 GB against DEEPARG_CONTIGS' 0.9 GB.

### What the derived route needs

- PYRODIGAL and DEEPARG_CONTIGS would run whenever `--arg_bin_clustering` is on, not only under `--contig_tax_and_arg` or `--contig_level_functional`.
- A small step maps DEEPARG_CONTIGS' hits to bins (contig name from the gene name, bin membership from the bin FASTA headers) and writes what ARG_FASTA_FORMATTER → CLUSTERING read today: per-bin `.mapping.ARG` files and per-bin protein FASTA (from Pyrodigal's proteins).
- Gene names change from Prodigal's to Pyrodigal's (both `<contig>_<n>`, but numbered by each caller), and 17 of the 228 shared hits are on genes with different coordinates, so the bin-level protein sequences for those differ slightly.
