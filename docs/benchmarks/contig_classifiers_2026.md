# Contig classifiers beyond Kraken2

Benchmark run 2026-10-28 to 2026-10-29, after [the contig taxonomy benchmark](contig_taxonomy_2026.md) chose Kraken2 with `k2_standard` for the contigs (L2). It asks what other tools add to that, on the same contigs:

- **(L5) Contig taxonomy beyond Kraken2:** (a) each binned contig inherits its bin's GTDB-Tk classification; (b) Taxometer refines Kraken2's labels with the contigs' composition and coverage; (c) Metabuli classifies the contigs against its RefSeq database.
- **(L6) Plasmid and phage classification:** (d) geNomad calls each contig chromosome, plasmid or virus, and finds proviruses; (e) RGI's pathogen-of-origin k-mers on the ARG-bearing contigs.

For each tool: runtime, CPU-seconds, peak memory, database size, contigs classified per rank or class, agreement on the ARG-bearing contigs, and what wiring it into the pipeline would take. No pipeline code changed.

**Decision (owner):** pending.

**Outcome (recommendations):**

- **(c) Metabuli: the one new classifier worth wiring, and a candidate to replace Kraken2 on the contigs.** It needs 5–12 GB of memory with `--max-ram 32` against Kraken2's 103 GB, with the same calls at every `--max-ram` tried (32, 64 and 128 GiB). It takes about the same CPU time (443–911 CPU-s per assembly against Kraken2's 221–1,393), and its database is the same size on disk (110 GB, an 80 GB download, against 104 GB). It gives a genus to more contigs: 73.9 % of all contigs against Kraken2's 62.1 %, and 236 of the 317 ARG-bearing contigs against 190, with `--min-score 0.2 --min-sp-score 0.5`. The trade-off is precision: where BLAST also names one, its genera agree with BLAST on 91.8 % of the ARG-bearing contigs (Kraken2: 95.2 %) and its species on 77.7 % (Kraken2: 88.8 %). Without score cutoffs it classifies 97 % of contigs, but its low-score calls are mostly wrong (below 0.5, species agree with BLAST 15–44 % of the time), so the cutoffs are needed. They were chosen on these four samples with BLAST as the reference, so the agreement figures at those cutoffs are an upper bound. Its database also has protozoa and fungi (Kraken2's PlusPF genome set): 442 contigs are *Trichomonas*, a trichomonad that BLAST confirms and Kraken2's standard database cannot name.
- **(a) Bin inheritance: add the bin and its GTDB classification as separate columns, without filling the NCBI columns from them.** It costs nothing new (GTDB-Tk already runs when binning is on) and labels 18 % of the contigs of A–C, including 173 of their 311 ARG-bearing contigs. When the GTDB and NCBI genera agree, the higher ranks differ on many contigs only because the two taxonomies place the same genus differently (same genus, same family on 63.5 % of contigs). The labels are least reliable where they would help most: on binned contigs Kraken2 gives no genus, the bin's genus agrees with BLAST for 119 of 347.
- **(b) Taxometer: do not wire.** It labels almost every contig, but its additions do not hold up: of the genera it adds where Kraken2 had none, 30 % agree with BLAST (35 % at score ≥ 0.95). Of the 12 labels it changes that BLAST can check, 11 match Kraken2's original label. It cannot predict what its input barely contains, so the eukaryotic and phage contigs Kraken2 leaves unclassified become Bacteria.
- **(d) geNomad: worth adding as columns, and as an ARGs-on-plasmids/phages summary, if 1.1–5.5 CPU-h per assembly is acceptable.** 34 of the 317 ARG-bearing contigs (50 of 403 ARG genes) are plasmid calls, 1 is a virus, and 2 more ARG genes lie inside proviruses. Several of the "plasmid" ARGs (tetM, ermG–mefA–mel, cfxA2) usually sit on integrative conjugative elements or transposons, which geNomad has no class for. geNomad explains a large part of the superkingdom disagreements between Kraken2 and BLAST (phage contigs: 41 of 88 in C, 7 of 18 in D), but not the contigs Kraken2 leaves unclassified (2.7–6.8 % of those are plasmid or virus calls).
- **(e) RGI pathogen-of-origin k-mers on contigs: do not wire.** `rgi kmer_query` needs RGI's own hits (`rgi main` on the contigs), and those cover 78 of the 317 ARG-bearing contigs. 44 of its 84 predictions name no taxon; where they do, the genus is what the taxonomy tools already say (33 of 33 with Kraken2). It never calls a gene plasmid-only: 55 of 84 are "chromosome or plasmid".

## Summary

| Question | Measured | Recommendation |
|---|---|---|
| (c) Metabuli cost | 59–100 s, 443–942 CPU-s per assembly at 16 threads; 5.3–12.3 GB anon at `--max-ram 32`, 16–23 GB at the default 128; database 110 GB (80 GB download). Kraken2 (9a): 221–1,393 CPU-s, 103 GB | Metabuli fits `process_medium` with no retry |
| (c) Metabuli results, all contigs | genus: 39,971 of 54,080 contigs (73.9 %) with the cutoffs, 47,590 without; Kraken2 33,570 (62.1 %). Genus agrees with Kraken2 on 99.4 % of contigs both classify | Metabuli, with `--min-score 0.2 --min-sp-score 0.5` |
| (c) Metabuli results, ARG-bearing contigs (317) | genus 236 (Kraken2 190, BLAST 278); agrees with BLAST at genus 202/220, species 129/166 (Kraken2 178/187, 111/125) | as above, with the trade-off stated |
| (a) Bin inheritance | 9,348 of 51,990 contigs of A–C binned, 173 of 311 ARG-bearing; GTDB vs Kraken2 genus 94.7 %; on binned contigs without a Kraken2 genus, 119 of 347 agree with BLAST | separate `bin` and GTDB columns |
| (b) Taxometer | 98–278 s, 850–3,202 CPU-s, < 0.8 GB; added genera 559 of 1,875 agree with BLAST; changed labels 1 of 12 | do not wire |
| (d) geNomad | 346–2,014 s, 3,927–19,730 CPU-s, 17–18 GB; database 2.2 GB (0.84 GB download). ARG-bearing contigs: 34 plasmid, 1 virus, 2 ARG genes in proviruses | add, if the cost is accepted |
| (e) RGI k-mers | 259–501 s, 271–539 CPU-s, 16.5 GB (the k-mer database load); RGI hits on 78 of 317 contigs; 40 of 84 k-mer predictions name a taxon | do not wire |
| mob-suite, platon | not run: geNomad left no gap the decision needs (see (d)) | an option for 9b if plasmid replicon types are wanted |

The rest of this document gives the measurements behind each row.

## Setup

**Contigs and inputs.** The same as the contig taxonomy benchmark: the four asm assemblies of the test fixture (Phase 7a's fixed MEGAHIT contigs ≥ 1 kb), the pipeline's bwa-mem3 BAMs, Phase 7b's Binette bins with their GTDB-Tk 2.8.0 classifications (GTDB r232), and that benchmark's Kraken2, BLASTN + BlobTools and DeepARG outputs.

| Sample | Contigs ≥ 1 kb | Mbp | Binette bins (contigs in bins) | ARG-bearing contigs (DeepARG genes) |
|---|---|---|---|---|
| A | 25,366 | 117.8 | 22 (4,476) | 142 (190) |
| B | 16,094 | 93.8 | 17 (2,482) | 132 (166) |
| C | 10,542 | 36.6 | 6 (2,390) | 37 (41) |
| D | 2,078 | 3.7 | none | 6 (6) |

**Reference for agreement.** BLASTN against nt with BlobTools (the contig taxonomy benchmark's runs) covers every contig of C and D and the ARG-bearing contigs of A and B. BLAST is a reference, not the truth: on the contigs Kraken2 cannot classify, its best hit is often a distant relative (below 82–85 % identity for a tenth of the ARG-bearing contigs). All NCBI names are compared on 2026 taxdumps: Kraken2's and BlobTools' (the k2_standard taxdump, June 2026) and Metabuli's (April 2026). 32 contigs whose Metabuli and Kraken2 genera differ only because a phage genus was renamed between the two dumps are counted as disagreements in the tables and noted where it matters. BlobTools' superkingdom is derived from its phylum through the June 2026 taxdump, because BlobTools 1.1.1 writes `undef` there with the 2026 taxdump (the contig taxonomy benchmark used the March 2024 taxdump's column instead, so its superkingdom counts differ slightly). BlobTools has no class rank.

**Tools.** Every tool ran in its biocontainer under Docker, with 16 threads unless stated, outside Nextflow, so that time and memory could be measured from the container's cgroup:

- **Metabuli** 1.2.0 (`quay.io/biocontainers/metabuli:1.2.0--pl5321h0bb26bb_0`), `classify --seq-mode 3 --lineage 1`, against `refseq_standard` (RefSeq genomes of Kraken2's PlusPF set: archaea, bacteria, viruses, plasmids, human, protozoa, fungi; NCBI taxonomy of 2026-04-03; built 2026-04-03). Mode 3 (long reads) is the single-sequence mode for contigs; on D, mode 1 (single-end reads) gave the same calls from superkingdom to species, and mode 3 placed 87 of them one level further, at strain. Each assembly ran three times: with the default `--max-ram 128` (container limit 160 GB), with `--max-ram 32` (limit 48 GB, so an overshoot would show and not be killed), and with `--max-ram 32 --min-score 0.2 --min-sp-score 0.5`; A also ran with `--max-ram 64` under 72 GB. `--precise` presets were not tried. The `gtdb232` database was not present, so only `refseq_standard` was measured.
- **Taxometer** (vamb 5.0.4, `quay.io/biocontainers/vamb:5.0.4--pyhdfd78af_0`), `vamb taxometer` on CPU, trained per sample on the contigs (`-m 1000`: vamb's default of 2,000 bp would drop about half of them), the sample's BAM, and Kraken2's labels (every contig, the lineage cut at Kraken2's first empty rank, names prefixed by rank so no name appears at two ranks; 232 of the 54,080 contigs lost ranks below a gap). `--seed 1`, default softmax threshold 0.5; results are also given at a score of 0.95.
- **geNomad** 1.13.0 (`quay.io/biocontainers/genomad:1.13.0--pyhdfd78af_0`), `end-to-end --cleanup` with default settings, database v1.9 (geNomad's release table pairs database 1.9 with geNomad 1.10–1.13). Per contig: the aggregated chromosome, plasmid and virus scores, and the plasmid and virus calls of its summary (default score thresholds).
- **RGI** 6.0.8 (the pipeline's RGI image), CARD 4.0.1 with WildCARD and its 61-mer database. `rgi kmer_query` does not read contigs: it reads RGI's own hits. So `rgi main -t contig -a DIAMOND --local` ran on the ARG-bearing contigs, then `rgi kmer_query --rgi -k 61 -n 2 -m 10 --local` (the pipeline's RGI_KMER settings) on its JSON.
- **Bin inheritance** is a join: contig names from each bin's FASTA, the bin's classification from GTDB-Tk's summary.

**Comparing GTDB with NCBI names.** GTDB suffixes (`Enterococcus_B` → *Enterococcus*) and NCBI brackets (`[Clostridium] innocuum` → *Clostridium innocuum*) are removed, and GTDB placeholder species (`Corynebacterium sp900539985`) count as no species. What still differs is split into placement (same genus, different higher taxon) and different genera; section (a) lists both.

## Machine and runs

- 4 × Xeon E7-4820 v4 (40 cores, 80 threads), 2 TB RAM; Docker.
- **One task at a time.** All runs went one after another, never two of these tasks at once. Other users' work (mostly a metaWRAP binning run) kept the 1-minute load at 58–88 on 80 threads throughout, so wall times are inflated.
- **Time and memory come from the container's cgroup**, sampled every 2 s (GB here are 2^30 bytes): CPU-seconds, peak cgroup memory (it includes page cache) and peak anonymous memory ("anon", the process memory).
- The Metabuli database had been written to disk shortly before, so part of it was probably still in the page cache: no cold-disk time was measured. The first Metabuli task (D) shows 55.7 GB of cgroup memory, which is page cache charged to it; its process memory was 16.0 GB.

## (c) Metabuli

### Cost

| | A | B | C | D |
|---|---|---|---|---|
| `--max-ram 128` (default): wall, CPU-s, anon | 100 s, 929, 22.9 GB | 98 s, 845, 21.6 GB | 82 s, 604, 17.9 GB | 83 s, 544, 16.0 GB |
| `--max-ram 32`: wall, CPU-s, anon | 85 s, 911, 12.3 GB | 81 s, 794, 10.9 GB | 64 s, 553, 7.2 GB | 59 s, 443, 5.3 GB |
| `--max-ram 32` with the score cutoffs | 96 s, 942, 12.2 GB | 89 s, 841, 10.9 GB | 69 s, 603, 7.2 GB | 64 s, 503, 5.3 GB |
| `--max-ram 64` | 90 s, 920, 15.8 GB | | | |
| Kraken2, k2_standard (contig taxonomy benchmark) | 58 s, 221, 103.6 GB | 67 s, 330, 103.5 GB | 151 s, 668, 103.3 GB | 301 s, 1,393, 103.0 GB |

- **Metabuli honours `--max-ram`, and the cap costs nothing here.** The calls were identical, row for row, at 32, 64 and 128 GiB, and the capped runs were no slower. Metabuli does not hold its 110 GB database in memory; it reads it from disk in parts.
- **CPU time is close to Kraken2's, for a different reason.** Kraken2's CPU-s are almost all the load of its 103 GB hash (classification itself took 1–6 s), so they shrink when the hash is in the page cache. Metabuli's grow with the assembly (443 CPU-s for 3.7 Mbp, 911 for 117.8 Mbp).
- **Memory is the difference:** at `--max-ram 32` Metabuli fits `process_medium` (8 CPUs, 36 GB) on the first attempt, where Kraken2 needs the 144 GB retry of `process_high`, clamped to `--max_memory`.
- Database: 110 GB on disk (`diffIdx` 52 GB, `info` 64 GB, taxonomy 0.8 GB); the download is 80 GB.

### Contigs classified

Contigs given a name at each rank, of all contigs:

| Rank | A: Kraken2 / Metabuli with cutoffs / without | B | C (and BLAST) | D (and BLAST) |
|---|---|---|---|---|
| Superkingdom | 17,525 / 20,031 / 24,737 | 11,225 / 12,277 / 15,659 | 8,380 / 9,006 / 10,307 (9,363) | 943 / 1,307 / 1,987 (1,364) |
| Genus | 15,608 / 18,859 / 22,435 | 9,381 / 11,274 / 13,794 | 7,840 / 8,681 / 9,681 (8,938) | 741 / 1,157 / 1,680 (1,239) |
| Species | 12,521 / 14,875 / 20,283 | 7,370 / 9,009 / 12,242 | 6,725 / 7,566 / 8,793 (9,553) | 677 / 757 / 1,519 (1,394) |

Over all four samples, Metabuli with the cutoffs names a genus for 39,971 of 54,080 contigs (73.9 %; Kraken2 62.1 %), and a species for 32,207 (59.6 %; Kraken2 50.5 %).

**Why the cutoffs.** Metabuli's defaults (`--min-score 0`) classify almost every contig, and the score (Metabuli's DNA-level identity score) says how far to trust it. Its median is 0.945 over all classified contigs, but only 0.12 on the 14,695 contigs Kraken2 leaves unclassified that Metabuli classifies. Against BLAST, on every contig where both name the rank:

| Metabuli score | < 0.1 | 0.1–0.2 | 0.2–0.3 | 0.3–0.5 | 0.5–0.7 | ≥ 0.7 |
|---|---|---|---|---|---|---|
| Genus agrees with BLAST | 48 % | 53 % | 61 % | 69 % | 86 % | 98 % |
| Species agrees with BLAST | 15 % | 34 % | 40 % | 44 % | 58 % | 94 % |

`--min-score 0.2` drops calls below 0.2, and `--min-sp-score 0.5` stops calls below 0.5 above species level. The pair was chosen from these tables, simulated, and then run: the real runs matched the simulation to within 26 contigs per sample. Since it was chosen on the same data it is measured on, the agreement at these cutoffs is an upper bound for other samples.

### Agreement

ARG-bearing contigs of all four samples (317), where both tools name the rank:

| Rank | Kraken2 vs BLAST | Metabuli with cutoffs vs BLAST | Metabuli without cutoffs vs BLAST | Metabuli with cutoffs vs Kraken2 |
|---|---|---|---|---|
| Superkingdom | 252/254 | 268/270 | 300/307 | 252/252 |
| Phylum | 233/239 | 259/264 | 289/300 | 235/237 |
| Order | 204/210 | 228/248 | 248/276 | 204/208 |
| Family | 191/201 | 216/235 | 227/257 | 196/198 |
| Genus | 178/187 (95.2 %) | 202/220 (91.8 %) | 211/239 (88.3 %) | 184/186 |
| Species | 111/125 (88.8 %) | 129/166 (77.7 %) | 150/220 (68.2 %) | 112/116 |
| Classified at genus | 190 | 236 | 266 | |

All contigs of C and D (BLAST covers them all): genus agrees with BLAST on 98.2 % of the contigs Kraken2 and BLAST both name (8,154 of 8,303) and 95.7 % for Metabuli with the cutoffs (8,852 of 9,253); species 93.9 % and 92.5 %.

**Where Kraken2 has nothing.** On the contigs of C and D that Kraken2 leaves without a genus but BLAST names one, Metabuli with the cutoffs names a genus for 306 of 537 in D (209 agree with BLAST, 68 %) and 698 of 1,337 in C (516 agree, 74 %). This is agreement with a nearest relative, not a measure of precision, for the reason given under Setup.

**Metabuli against Kraken2.** Where both name a genus, they agree on 33,108 of 33,317 contigs (99.4 %). 32 of the other 209 are one phage genus renamed between the two taxdumps (*Taduhovirus* → *Aurodevirus*, same taxid). The rest are different calls: *Prevotella* against *Segatella* (17), *Flavonifractor* against *Dysosmobacter* (8), *Bacteroides* against *Phocaeicola* or *Parabacteroides* (13), and smaller sets.

**Eukaryotes.** With the cutoffs Metabuli calls 543 contigs eukaryotic: 442 *Trichomonas vaginalis* (0.88 Mbp; 270 in C, 106 in A, 50 in B, 16 in D) and 99 *Homo sapiens*. Kraken2 names 109 contigs human and none trichomonad. BLAST ran on 286 of the trichomonad contigs: its best hit is *T. vaginalis* for 250, and the other 36 have no BLAST species. That makes them trichomonad DNA, not necessarily that species. This is the genome set, not the tool: k2_standard has no protozoa or fungi, and Kraken2's PlusPF database (not on the server) would include them.

**Phage contigs.** Metabuli does no better than Kraken2 on the 1,448 contigs geNomad calls viral: with the cutoffs it names 916 Bacteria, 41 Viruses and 31 Eukaryota, and leaves 455 unclassified (Kraken2: 830 Bacteria, 41 Viruses, 572 unclassified).

## (a) Bin inheritance (GTDB-Tk, GTDB r232)

Every contig in a Binette bin takes the bin's GTDB-Tk classification. All 45 bins of A–C have one down to species (GTDB placeholder species included).

| Sample | Contigs | Binned | Binned, no Kraken2 genus | ARG-bearing | ARG-bearing, binned | ARG-bearing, binned, no Kraken2 genus |
|---|---|---|---|---|---|---|
| A | 25,366 | 4,476 | 1,410 | 142 | 84 | 28 |
| B | 16,094 | 2,482 | 349 | 132 | 74 | 13 |
| C | 10,542 | 2,390 | 479 | 37 | 15 | 6 |

**Against Kraken2** (all binned contigs, where both name the rank): phylum 7,296 of 7,330 (99.5 %), class 90.6 %, order 82.3 %, family 63.4 %, genus 6,735 of 7,110 (94.7 %), species 4,907 of 5,309 (92.4 %).

**Placement, not disagreement.** On the 6,735 contigs where Kraken2 and GTDB give the same genus, the phylum always agrees, the class on 90.5 %, the order on 82.2 % and the family on 63.5 %. Those differences are the two taxonomies placing the same genus differently: *Odoribacter* and *Butyricimonas* in Odoribacteraceae (NCBI) or Marinifilaceae (GTDB), *Ruthenibacterium* in Oscillospiraceae or Ruminococcaceae, *Clostridium innocuum* in Clostridia or Bacilli (Erysipelotrichaceae), *Corynebacterium* in Corynebacteriaceae or Mycobacteriaceae, *Sutterella* in Betaproteobacteria or Gammaproteobacteria, *Prevotella* in Prevotellaceae or Bacteroidaceae.

**Different genera** (375 contigs): 37 are a naming split (NCBI's *Segatella*, separated from *Prevotella* in 2023, against GTDB r232's *Prevotella*). The others are contigs whose k-mers point to another genus than their bin's: *Bacteroides* contigs in a *Phocaeicola dorei* bin (146), *Sellimonas* in a *Merdimonas* bin (17), *Odoribacter* in a *Culturomica* bin (12), and smaller sets. Some are shared mobile DNA, some binning errors; the benchmark cannot tell which.

**Against BLAST:** on the ARG-bearing contigs, genus agrees on 131 of 155 (84.5 %) and species on 107 of 149; on all binned contigs of C, genus on 1,819 of 2,210 (82.3 %, *Segatella*/*Prevotella* included).

**Where Kraken2 gives no genus.** 2,238 binned contigs of A–C (47 of them ARG-bearing) get a genus only from their bin. These are the contigs where the label would help, and where it is weakest: the bin's genus agrees with BLAST on 119 of 347 such contigs, and with Metabuli (with the cutoffs) on 325 of 759. 277 binned contigs are geNomad virus calls and 254 plasmid calls, which a bin label describes poorly.

## (b) Taxometer

| | A | B | C | D |
|---|---|---|---|---|
| Wall, CPU-s, anon | 278 s, 3,202, 0.72 GB | 236 s, 2,195, 0.59 GB | 200 s, 1,649, 0.65 GB | 99 s, 850, 0.53 GB |

Taxometer's labels against its input (Kraken2), all four samples, at its default threshold (0.5) and at a score of 0.95:

| Rank | Kraken2 | Taxometer (0.5 / 0.95) | Added (0.5 / 0.95) | Added, BLAST names the rank: agree | Changed (0.5 / 0.95) | Changed, BLAST checks: agrees with Taxometer / with Kraken2 | Removed (0.5 / 0.95) |
|---|---|---|---|---|---|---|---|
| Phylum | 37,720 | 54,048 / 52,568 | 16,330 / 14,876 | 1,061 of 1,683 / 1,011 of 1,514 | 2 / 0 | 0 / 0 | 2 / 28 |
| Family | 35,604 | 53,746 / 50,174 | 18,145 / 14,660 | 950 of 1,859 / 861 of 1,519 | 7 / 0 | 1 / 0 | 3 / 90 |
| Genus | 33,570 | 53,096 / 46,741 | 19,532 / 13,380 | 559 of 1,875 (30 %) / 456 of 1,317 (35 %) | 246 / 167 | 1 / 11 (at 0.5); 0 / 5 (at 0.95) | 6 / 209 |
| Species | 27,293 | 52,736 / 44,290 | 25,452 / 17,678 | 1,097 of 3,578 / 1,039 of 2,582 | 279 / 167 | 0 / 11; 0 / 5 | 9 / 681 |

On the ARG-bearing contigs, Taxometer adds a genus to 120 (0.5) or 74 (0.95), and 29 of 85 and 25 of 55 of those agree with BLAST.

- **Its additions are a guess from the sample's own labelled contigs.** It can only predict labels that are frequent in its input. Of its superkingdom additions that BLAST can check, 510 of 1,601 disagree, and 509 of those are contigs BLAST calls eukaryotic (246, mostly the trichomonad) or viral (263), which Taxometer calls Bacteria (503), Archaea (4) or Eukaryota (2).
- **On the same hard contigs Metabuli does better.** Contigs of C and D without a Kraken2 genus but with a BLAST genus: Taxometer's genus agrees on 195 of 517 (D) and 341 of 1,291 (C), Metabuli with the cutoffs on 209 of 306 and 516 of 698.
- **Its changes are almost always wrong where they can be checked**, and they are few.

## (d) geNomad

### Cost

| | A | B | C | D |
|---|---|---|---|---|
| Wall, CPU-s, anon | 2,014 s, 19,730, 18.4 GB | 1,529 s, 15,244, 18.2 GB | 823 s, 8,647, 17.7 GB | 346 s, 3,927, 17.4 GB |

160–240 CPU-s per Mbp of contigs on A–C (D, at 1,060 per Mbp, shows a fixed part), 1.1–5.5 CPU-h per assembly; memory does not grow with the assembly. On A the neural-network step alone took about 20 min. Database v1.9: 2.2 GB unpacked, a 0.84 GB download. nf-core has `genomad/endtoend` (label `process_high`, geNomad 1.12.0, which also takes database 1.9) and `genomad/download`.

### Calls

| Sample | Contigs | Plasmid | Virus | Provirus regions (contigs) | ARG-bearing contigs: plasmid / virus | ARG genes on plasmid contigs | ARG genes inside a provirus |
|---|---|---|---|---|---|---|---|
| A | 25,366 | 402 | 605 | 61 (60) | 9 / 1 | 16 of 190 | 2 (vanU, twice) |
| B | 16,094 | 382 | 437 | 50 (49) | 17 / 0 | 24 of 166 | 0 |
| C | 10,542 | 209 | 360 | 7 (6) | 7 / 0 | 9 of 41 | 0 |
| D | 2,078 | 18 | 46 | 0 | 1 / 0 | 1 of 6 | 0 |

Classes on the 34 ARG-bearing plasmid contigs (a contig counts once per class): MLS 10, glycopeptide 9 (*van* operon genes), beta-lactam 4 (all cfxA2), tetracycline 4, multidrug 4 (qacA, qacG, efrA/efrB, marA/ompR), phenicol 2, and one each of bacitracin, fosfomycin and aminoglycoside. The one virus-called ARG-bearing contig carries bcrA. 18 ARG-bearing contigs carry a provirus region, but only two ARG genes (vanU, in A) lie inside one. geNomad's own AMR markers flag 5 of the 34 ARG-bearing plasmid contigs.

**"Plasmid" includes other mobile elements.** tetM, ermG–mefA–mel and cfxA2 are typically carried on integrative conjugative elements and transposons (Tn916-like, CTnDOT-like, MTn4555-like), and geNomad has no class for those: a contig with a conjugation module scores as plasmid. So read "plasmid" as "mobile element, probably plasmid". Replicon typing (mob-suite) could separate true plasmids, but whether to add geNomad does not depend on it, so it was not run.

**Virus taxonomy** is ICTV's, not NCBI's: 1,329 of the 1,448 virus calls are *Caudoviricetes* (tailed phages), 33 *Herviviricetes*, 76 unclassified.

### Against the taxonomy results

- **Disagreements between Kraken2 and BLAST** (C and D, all contigs): at superkingdom, 41 of 88 in C and 7 of 18 in D are virus calls; at genus, 26 of 115 in C (17 plasmid, 9 virus) and 1 of 34 in D. Phage contigs carry genes shared with their hosts, so one tool calls the host and the other the phage.
- **Contigs Kraken2 leaves unclassified are mostly chromosomal:** plasmid or virus calls are 2.7 % of them in D (31 of 1,134), 6.8 % in C (147 of 2,159), 3.5 % in B and 4.2 % in A. They are chromosomal DNA of genomes not in the database, not unclassified mobile elements.
- **Classifiers call phages Bacteria.** Of the 1,448 virus calls, Kraken2 names 41 as Viruses and 830 as Bacteria; Metabuli the same 41 and 916. Both databases contain phages within bacterial genomes (prophages), and for most phage contigs that is the nearest reference.

## (e) RGI pathogen-of-origin k-mers

| | A | B | C | D |
|---|---|---|---|---|
| Wall, CPU-s, anon (`rgi main` + `kmer_query`) | 501 s, 539, 16.5 GB | 365 s, 399, 16.5 GB | 266 s, 283, 16.5 GB | 259 s, 271, 16.8 GB |
| RGI hits (Strict / Perfect) on ARG-bearing contigs | 56 / 4 on 34 of 142 | 38 / 11 on 32 of 132 | 10 / 1 on 11 of 37 | 1 / 0 on 1 of 6 |

CPU-s are about equal to wall time (one thread at work), and D, with a single RGI hit, still took 259 s: most of the time is fixed cost, the 61-mer database load among it.

- **Coverage:** RGI's Strict and Perfect hits cover 78 of the 317 ARG-bearing contigs, so the k-mers say nothing about the other 239. 84 of the 121 hits got a k-mer prediction.
- **Taxonomy:** 44 of 84 predictions are "Unknown taxonomy". Where they name a taxon (40), the genus agrees with Kraken2 on 33 of 33, with Metabuli (with the cutoffs) on 34 of 37, with BLAST on 36 of 39 and with the bin's GTDB genus on 27 of 28. On contigs, the k-mers add nothing a contig classifier does not already give.
- **Chromosome or plasmid:** 55 of 84 are "chromosome or plasmid", 28 "chromosome", 1 "no genomic info", none plasmid-only. geNomad calls 21 of the 55 and 5 of the 28 plasmid. It cannot replace geNomad here.

mob-suite and platon (both local) were not run: geNomad answered the question the decision needs (how many ARGs sit on plasmids or phages), and RGI's k-mers do not.

## Wiring (what each option would take)

All changes would go into Phase 9b, with L2 and M8.

**(c) Metabuli on the contigs**, replacing Kraken2 there, or alongside it:
- Process: a local METABULI_CLASSIFY module (no nf-core module), `classify --seq-mode 3 --lineage 1 --threads ${task.cpus} --max-ram` from the task's memory (32 for `process_medium`'s 36 GB), `--min-score` and `--min-sp-score` as parameters (defaults 0.2 and 0.5), image `metabuli:1.2.0--pl5321h0bb26bb_0`. `process_medium` at 8 CPUs would take about twice the wall time measured here at 16.
- Database: a `metabuli_ref_db` registry entry (`refseq_standard`, 80 GB download, 110 GB unpacked; `gtdb232`, 634 GB download, 744 GB unpacked, not measured) and `--custom_metabuli_db`.
- Columns in `Contig_tax_and_arg_prediction.tsv`: the seven ranks from Metabuli's lineage (with class, which BlobTools lacks), plus `taxonomy_rank` and `taxonomy_score`. BlobTools' tables and plots would take Metabuli's per-contig taxid as their hits file, as planned for Kraken2.
- Summary: classified contigs per rank per sample, and the share of eukaryotic and unclassified contigs.
- If Metabuli replaces Kraken2 on the contigs, the 9b memory fix for KRAKEN2 on contigs is no longer needed, and the read-level KRAKEN2 keeps its own database.

**(a) Bin inheritance:**
- No new tool: a step after GTDB-Tk (when `--include_binning` is on) joins bin membership (the bin FASTA headers, already used by the M8 change) with GTDB-Tk's summary.
- Columns: `bin`, `bin_gtdb_classification` (GTDB names, untouched), next to the NCBI columns rather than merged into them.
- Summary: ARG-bearing contigs per bin species.

**(d) geNomad:**
- Process: nf-core `genomad/endtoend` (bump to 1.13.0, or keep its 1.12.0; both take database 1.9) under a new flag, or always with `--contig_tax_and_arg`. Its `process_high` (16 CPUs, 72 GB) is more memory than needed (18 GB measured); `process_medium` would fit at about twice the wall time.
- Database: `genomad_ref_db` registry entry (v1.9, 0.84 GB download) through nf-core `genomad/download`, and `--custom_genomad_db`.
- Columns: `mobile_class` (chromosome / plasmid / virus), `plasmid_score`, `virus_score`, `provirus` (coordinates when a provirus region overlaps the ARG gene), `virus_taxonomy`, `plasmid_conjugation_genes`.
- Summary: ARG genes on plasmid and virus contigs and inside proviruses, per sample and per class, with the "plasmid includes ICEs and transposons" caveat.
- geNomad's scores use marker frequencies across its input, so running it only on the ARG-bearing contigs could change them; that was not tested.

**(b) Taxometer** and **(e) RGI k-mers on contigs:** not recommended. Taxometer would need a vamb process per sample (the BAM and the Kraken2 or Metabuli labels as input; under 1 GB, 850–3,200 CPU-s). RGI on contigs would need `rgi main -t contig` plus `kmer_query --rgi` per sample (16.5 GB for the k-mer database).
