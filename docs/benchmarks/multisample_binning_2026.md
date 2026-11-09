# Multi-sample binning and run-level dereplication

Benchmark run 2026-11-08 (roadmap Phase 18a, items B4 and B5). No pipeline code changed. It asks two questions about per-sample mode, the default, in which each sample's assembly is binned with that sample's coverage only:

- **B4, multi-sample binning:** does giving each assembly the coverage of every sample in the run give better MAGs, which way of doing it is best, and what does it cost as the number of samples grows?
- **B5, run-level dereplication:** how dRep and galah cluster the run's MAGs at 95 % ANI (species) and 99 % (strain), which representatives they choose, and how to get each representative's abundance in every sample, with multi-sample binning on and with it off.

Both are run-level results. A study split across batches dereplicates its MAGs downstream, from the merged results, not in the pipeline (roadmap decision "Study-level analyses", 2026-11-05).

**Decision (owner, 2026-11-08):** add both as recommended below, optional and off by default. Multi-sample binning in the concatenation form (MetaBAT2 gets the extra depth columns; SemiBin2 as today), cross-mapped inside the groups of a metadata column, with a launch warning when a group has more samples than a set number. COMEBin with the same columns is measured in Phase 18b before it is enabled. Dereplication (galah, then CoverM) may run without multi-sample binning.

**Outcome (recommendations):**

- **B4: add multi-sample binning as an option, off by default, in the cheapest form measured.** In that form, each sample is aligned once more, to the concatenation of all assemblies in its group, and MetaBAT2 gets one depth column per sample from those alignments. SemiBin2 keeps its pre-trained model and the sample's own BAM, as today. On these four samples the refined bins went from 27 to 29 CheckM2 high-quality (> 90 % / < 5 %) and from 3 to 4 MIMAG high. Sequence in passing bins rose from 158.3 to 161.2 Mbp, and Σ(completeness − 5 × contamination) over all bins from 3,631 to 3,731. All 31 species were kept. The cost is one extra alignment per sample: 2.1 CPU-hours and 3.9 GiB of BAMs for the four samples, on top of today's 1.7 CPU-hours for alignment, ~40 % of today's alignment-plus-binning CPU when Binette is counted on the same basis (Cost). Binning CPU is unchanged. This is the same split the pipeline already uses under co-assembly (pre-trained SemiBin2, per-sample depth columns for MetaBAT2; [binning_2026.md](binning_2026.md), (d) and (g)).
- **The gain is modest and comes from MetaBAT2.** No variant found a genome that today's run misses. Every refined bin of every variant matches one of today's bins (Jaccard index ≥ 0.5 over contig base pairs), and the GTDB species are the same 31 (30 in (iii), which lost one). Multi-sample coverage completes existing bins. Raw MetaBAT2 high-quality bins go from 14 to 24–25, while SemiBin2 moves from 24 to 26–27, and only when self-trained.
- **Not variant (ii) as specified (N × N alignments + self-trained SemiBin2).** It gave the same refined quality (29 high, Σ 3,745) for 28.3 CPU-hours against 5.8 today. SemiBin2's self-training is 18.6 of those, and the N × N alignments grow with the square of the samples: 400 alignments at 20 samples, 2,500 at 50.
- **Not variant (iii), SemiBin2 `multi_easy_bin`.** Same quality (30 high), but the binning is one task: 4 h 41 min of wall time at N = 4, training one model per sample one after another. Its wall time grows with the number of samples, and Nextflow cannot spread it over tasks. It also lost one species (*Bacteroides ovatus*). (The roadmap's "one model" is not what `multi_easy_bin` does: SemiBin2 2.6.0 trains a separate self-supervised model for each sample, and refuses `--train-from-many` in multi-sample mode.)
- **Groups:** cross-mapping happens inside each group of a metadata column, as Phase 15's co-assembly groups do. A group bounds the concatenation's size, which drives index memory and splits reads between copies of the same genome. It also bounds the alignment count of the N × N alternative.
- **B5: dereplication can run on its own (with multi-sample binning on or off), with galah.** dRep and galah gave the same clusters in every case: 45 MAGs → 31 at 95 % ANI, 32 at 99 %. With `--quality-formula completeness-5contamination`, galah chose exactly dRep's representatives, in seconds and with no writes outside its work directory. galah's default formula (`Parks2020_reduced`) penalises contig count and in one cluster passed over a MIMAG-high bin for a medium one.
- **B5 abundance: CoverM on the representatives.** Its mean depth matches the binning depth (BIN_DEPTH's length-weighted definition) to within 0.5 % in three samples and 2 % in the fourth. It costs ~1.3 CPU-hours and 14.5 GB for four samples. With the recommended B4, CoverM is needed: the concatenation's depths are competitive, 0.32–0.97 of a sample's true coverage (section B4). Only the N × N variant has every sample's depth on every assembly; with it, BIN_DEPTH would give the abundance without CoverM.

## Summary

| Question | Measured (4 samples) | Recommendation |
|---|---|---|
| Multi-sample coverage, any form | 0 new genomes, same 31 species (30 in (iii)); CheckM2 high 27 → 27–30; Σ quality +50 to +114 | worth it only cheaply; optional |
| (ii) N × N + self-trained SemiBin2 | 29 high, 28.3 CPU-h (SemiBin2 18.6), 16 BAMs 14.9 GiB | no |
| (iii) `multi_easy_bin` | 30 high, one 4.7 h task, one species lost | no |
| N × N + pre-trained SemiBin2 (diagnostic) | 27 high, Σ 3,744, 10.5 CPU-h, 14.9 GiB | no (N² alignments for the same gain) |
| **Concatenation coverage for MetaBAT2 + pre-trained SemiBin2 (diagnostic)** | **29 high, 4 MIMAG high, Σ 3,731, 7.3 CPU-h, 7.2 GiB** | **yes, optional, inside groups** |
| Dereplication | dRep = galah partitions; galah ~2–3 s, dRep 11–21 s | galah, formula `completeness-5contamination`; runs alone |
| Abundance | CoverM vs BIN_DEPTH: ratio 1.00–1.02, Spearman ≥ 0.988 | CoverM on the representatives |

## Setup

**Contigs and reads.** The Phase 11a per-sample reference (four samples, dev `c38ff92`), so every variant bins the same contigs and MEGAHIT's run-to-run differences do not enter. The clean reads are the pipeline's QC output. Samples are named as in [binning_2026.md](binning_2026.md):

| Sample | Contigs ≥ 1 kb | Mbp | Clean read pairs |
|---|---|---|---|
| A | 25,366 | 117.8 | 14.3 M |
| B | 16,094 | 93.8 | 27.7 M |
| C | 10,542 | 36.6 | 2.9 M |
| D | 2,078 | 3.7 | 1.4 M |

**Harness.** A scratch Nextflow workflow imported the pipeline's own BWAMEM3_INDEX, BWAMEM3_SAMTOOLS, CALCULATE_DEPTH, METABAT2, SEMIBIN and BINETTE modules from dev `3743f39`, unchanged. The repository's `nextflow.config` and `config/modules.config` applied. Three scratch processes ran `multi_easy_bin`: one concatenated the assemblies (`SemiBin2 concatenate_fasta`, contig names prefixed `<sample>:`), one ran `SemiBin2 multi_easy_bin --self-supervised --separator : --random-seed 11` with 4 threads, and one split the concatenation's depth table back into one table per assembly. Every variant was refined by Binette per assembly from SemiBin2 + MetaBAT2, with its CheckM2 reports, as the pipeline does.

**Variants** (all bin the same per-sample contigs):

| Variant | Alignments | MetaBAT2 depth columns | SemiBin2 |
|---|---|---|---|
| (i) today | each sample to its own assembly (N) | 1, own sample | pre-trained `global`, own BAM (the 11a run) |
| self (control) | N, the same BAMs | 1 | self-trained, own BAM |
| (ii) cross | every sample to every assembly (N × N) | N | self-trained on N BAMs (`train_self`; a pre-trained model takes one coverage column) |
| (iii) multi | every sample to the concatenation (N), plus today's N | N, from the concatenation | `multi_easy_bin`, self-supervised |
| hybrid (diagnostic) | N × N | N, from the N × N BAMs | pre-trained `global`, own BAM |
| hybrid-concat (diagnostic) | concatenation (N) plus today's N | N, from the concatenation | pre-trained `global`, own BAM |

The self control separates the effect of self-training from that of coverage. The two diagnostics test whether MetaBAT2's extra columns alone, with SemiBin2 left as it is, carry the gain. They were added once (ii) and (iii) were in.

**Determinism.** The four diagonal BAMs (sample to own assembly) are byte-identical in content to Phase 11a's (`samtools view | md5sum`). MetaBAT2 on those BAMs gave bins identical to 11a's on every sample (the pipeline seeds it, `--seed 3`). SemiBin2 was seeded (`--random-seed 11`, the pipeline's `--semibin_random_seed`). The differences below are therefore the variants', not run-to-run noise. Binette was not rerun for a repeat.

**Quality.** CheckM2 from Binette's reports. MIMAG tiers use the pipeline's rules (`bin/summary_mags.py`). For high quality they need > 90 % / < 5 % plus 16S, 23S and 5S rRNA and ≥ 18 tRNA amino acids. The rRNA and tRNA came from Bakta 1.12.2 (`bakta:1.12.2--pyhdfd78af_0`, database `updated_db`) with every step but tRNA and rRNA skipped, run on each distinct bin above 90 % / 5 %. On today's bins this gives 3 high / 42 medium, as in 11a's summary. GTDB-Tk 2.8.0 (`gtdbtk:2.8.0--pyhdfd78af_0`, R232, 16 CPUs, `--pplacer_cpus 4`) classified the 236 distinct genomes of all variants in three batches (172 genomes in 9.5 min and 45.5 GB, then 32 and 32 for the diagnostics): all bacterial, all to species.

**Matching bins.** Each bin is paired with today's bin of the same assembly that shares the most sequence (Jaccard over contig base pairs). Genomes count as distinct by their sequence (md5 of the sorted contigs).

**Cost.** Nextflow's trace for the binning tasks (CPU-hours = realtime × %CPU). Each B5 tool ran in one container with its cgroup sampled every 2 s. The host (80 threads) was shared with other work at load 35–65, so wall times are indicative.

## Machine and runs

One 80-thread host, 2 TB RAM, docker. Bench run 14:38–19:59 UTC (72 tasks, none failed). The diagnostics, B5 and GTDB-Tk ran after it.

## B4. Multi-sample binning

### Refined bins

Per assembly (sample D gave no bin in any variant):

| Variant | Bins A / B / C | CheckM2 high A / B / C | MIMAG high A / B / C |
|---|---|---|---|
| (i) today | 22 / 17 / 6 | 13 / 12 / 2 | 2 / 1 / 0 |
| self | 22 / 17 / 5 | 12 / 10 / 1 | 1 / 1 / 0 |
| (ii) cross | 22 / 17 / 6 | 13 / 15 / 1 | 2 / 1 / 0 |
| (iii) multi | 22 / 16 / 6 | 14 / 14 / 2 | 1 / 1 / 0 |
| hybrid | 22 / 17 / 6 | 12 / 14 / 1 | 3 / 0 / 0 |
| hybrid-concat | 22 / 17 / 6 | 13 / 14 / 2 | 2 / 2 / 0 |

In total:

| Variant | Bins | CheckM2 high | MIMAG high / medium | Median completeness / contamination | Σ (compl. − 5 × cont.) | Mbp | GTDB species (in high bins) |
|---|---|---|---|---|---|---|---|
| (i) today | 45 | 27 | 3 / 42 | 94.4 / 1.01 | 3,631 | 158.3 | 31 (20) |
| self | 44 | 23 | 2 / 42 | 94.8 / 0.72 | 3,597 | 152.7 | 31 (18) |
| (ii) cross | 45 | 29 | 3 / 42 | 96.9 / 0.79 | 3,745 | 157.2 | 31 (21) |
| (iii) multi | 44 | 30 | 2 / 42 | 97.8 / 0.78 | 3,681 | 152.6 | 30 (22) |
| hybrid | 45 | 27 | 3 / 42 | 96.8 / 0.78 | 3,744 | 159.3 | 31 (20) |
| hybrid-concat | 45 | 29 | 4 / 41 | 96.7 / 1.23 | 3,731 | 161.2 | 31 (21) |

- **Self-training alone makes things worse.** With one coverage column, self-trained SemiBin2 gave 4 fewer high-quality bins than the pre-trained model (23 vs 27). The pre-trained per-sample baseline is already strong. Measured against the self control, multi-sample coverage is worth +6 or +7 high-quality bins; against today, +2 or +3.
- **MIMAG high hardly moves** (2–4 in every variant). It depends on the rRNA operons, which short-read assemblies rarely keep in a bin. The CheckM2 count shows the differences.
- **The raw binners** (Binette's CheckM2 reports; bins / ≥ 50 % and < 10 % / high):

  | Variant | MetaBAT2 | SemiBin2 |
  |---|---|---|
  | (i) today | 86 / 38 / 14 | 93 / 45 / 24 |
  | self | 86 / 38 / 14 | 95 / 43 / 23 |
  | (ii) cross | 65 / 44 / 25 | 89 / 43 / 27 |
  | (iii) multi | 68 / 43 / 24 | 92 / 41 / 26 |

  MetaBAT2 nearly doubles its high-quality bins with N columns, whether the depths come from N × N alignments or from the concatenation. This is the pattern 7a found under co-assembly (11 → 23).

### The same genomes, completed

Every bin of every variant matches one of today's bins at a Jaccard index ≥ 0.5 over contig base pairs (none new). The share at ≥ 0.95 ("same") is 25–26 of 45 for the two diagnostics, 27 of 44 for self, and 12–13 for cross and multi. Bins move most when self-trained SemiBin2 gets N coverage columns ((ii) and (iii)); self-training on the sample's own BAM alone moves them about as much as the diagnostics. Anyone comparing MAGs from runs with and without the option will see bin boundaries change.

Bins crossing the 90 % / 5 % line against today:

| Variant | Gained | Lost |
|---|---|---|
| (ii) cross | B *B. stercoris* 74.5 → 97.9 %; B *P. distasonis* 85.4 → 99.9 %; B 89.1 / 3.3 → 96.9 / 2.4; A 86.2 → 90.3 % | C *C. innocuum* contamination 4.4 → 6.0 %; A *M. torques* 99.1 → 87.2 % |
| (iii) multi | the same three in B, two in A | the same *M. torques*; a B bin 96.1 → 69.7 % |
| hybrid | B *B. stercoris*, B 89.1 → 96.8 % | the same two as cross |
| hybrid-concat | B *B. stercoris* 74.5 → 99.3 %, B 89.1 → 97.0 % | none |

The *C. innocuum* bin of sample C sits on the 5 % contamination line in every variant (4.4–6.0 %). *B. ovatus* (sample B, 79 % complete, 664 contigs) has no counterpart in (iii).

### Cost (four samples)

| Variant | Alignments | Alignment CPU-h | BAMs (GiB) | SemiBin2 CPU-h (longest task) | Depth + MetaBAT2 CPU-h | Binette CPU-h | Total CPU-h | Peak RSS |
|---|---|---|---|---|---|---|---|---|
| (i) today | 4 | 1.73 | 3.3 | 0.84 (15 min, 8 CPUs) | 0.06 | 3.02 | 5.8 | 15.8 GB (alignment) |
| self | 4 | 1.73 | 3.3 | 18.0 (2 h 9 min) | 0.05 | 2.45 | 22.3 | 19.7 GB (SemiBin2) |
| (ii) cross | 16 | 6.98 | 14.9 | 18.6 (2 h 16 min) | 0.12 | 2.44 | 28.3 | 19.6 GB (SemiBin2) |
| (iii) multi | 4 + 4 | 3.87 | 7.2 | 16.7 (one task, 4 h 41 min) | 0.06 | 2.42 | 23.2 | 20.0 GB (SemiBin2) |
| hybrid | 16 | 6.98 | 14.9 | 0.84 | 0.12 | 2.45 | 10.5 | 16.0 GB (alignment) |
| hybrid-concat | 4 + 4 | 3.87 | 7.2 | 0.84 | 0.06 | 2.35 | 7.3 | 16.7 GB (alignment) |

The totals include the indexes: 0.10 CPU-hours for the four assemblies, and 0.09 CPU-hours and 2.3 GB for the 252 Mbp concatenation. Today's Binette ran in 11a under a different load; the bench's own Binette runs took 2.35–2.45 CPU-hours. Counted with the bench's Binette (≈ 2.4), today is ~5.2 CPU-hours, and hybrid-concat costs 2.1 CPU-hours (~40 %) more, the extra alignments. Variants (iii) and hybrid-concat count today's N alignments as well. The pipeline needs each sample's BAM against its own assembly anyway: BIN_DEPTH, featureCounts and the other downstream steps read it.

**An alignment's cost is set by the sample's reads, not by the reference.** Sample B took 0.96–1.09 CPU-hours against any one assembly, and 1.25 against the four-assembly concatenation (CPU-hours per alignment):

| Reads of → aligned to | A | B | C | D | concatenation |
|---|---|---|---|---|---|
| A | 0.54 | 0.59 | 0.56 | 0.46 | 0.66 |
| B | 1.09 | 1.01 | 1.08 | 0.96 | 1.25 |
| C | 0.13 | 0.13 | 0.13 | 0.11 | 0.17 |
| D | 0.05 | 0.05 | 0.05 | 0.05 | 0.06 |

### How the cost grows with N

Taking these samples' mean per alignment (0.43 CPU-hours and a 0.83 GiB BAM against the own assembly, 0.44 and 0.97 against another; 0.54 and 0.98 against the concatenation). Real samples with more reads scale both up in proportion:

| Samples | Today: alignments, CPU-h, BAM GiB | N × N | N × N in groups of 5 | Concatenation, one group | Concatenation in groups of 5 |
|---|---|---|---|---|---|
| 4 | 4, 1.7, 3.3 | 16, 7.0, 15 | (one group) | 8, 3.9, 7.2 | (one group) |
| 20 | 20, 8.6, 17 | 400, 175, 385 | 100, 44, 94 | 40, ≥ 19, 36 | 40, 19, 36 |
| 50 | 50, 22, 42 | 2,500, 1,094, 2,420 | 250, 109, 236 | 100, ≥ 48, 91 | 100, 48, 91 |

- N × N grows with the square of N. Groups of g samples make it g × N; each sample's BAMs wait in the work directory until binning has read them.
- The concatenation keeps the count at 2N for any grouping. Without groups it has two problems that grow with N, neither measured beyond N = 4:
  - **Index memory.** The reference grows with every assembly: 2.3 GB to index 252 Mbp here. Fifty assemblies of this size would be ~3 Gbp, about a human genome.
  - **Split reads.** A read from a genome assembled in k samples is placed in only one of the k copies. Against the four-assembly concatenation, the median depth of binned contigs was 0.32–0.97 of the same sample's depth from the N × N alignments. For sample C's assembly, whose genomes all recur in A or B, it was 0.33–0.36 in every sample.

  MetaBAT2 bins on the depth profile across samples, which still worked at N = 4. These depths are not each sample's coverage, though, and must not be used as abundance. Groups bound both problems.
- Binning cost barely depends on N in the recommended form: SemiBin2 still bins each assembly from its own BAM. CALCULATE_DEPTH reads N BAMs against the group's concatenation and MetaBAT2 gets N columns, so both grow with N, but they took minutes (0.06 CPU-hours for N = 4, as today). With (iii), the single `multi_easy_bin` task would train N models one after another (4 h 41 min for N = 4 at 4 threads).

## B5. Run-level dereplication and abundance

### Clusters and representatives

dRep 4.0.3 (`dereplicate --genomeInfo <CheckM2> -comp 50 -con 10 -sa 0.95|0.99 --skip_plots`, skani for both steps) and galah 0.6.2 (`cluster --checkm2-quality-report --min-completeness 50 --max-contamination 10 --ani 95|99`, skani), on each variant's refined MAGs:

| Variant | MAGs | Clusters 95 % / 99 % | Same partition (dRep, galah) | Clusters with > 1 MAG | Same representative, galah default / `completeness-5contamination` |
|---|---|---|---|---|---|
| (i) today | 45 | 31 / 32 | yes / yes | 10 | 29 / 31 of 31 |
| self | 44 | 31 / 32 | yes / yes | 10 | 30 / 31 of 31 |
| (ii) cross | 45 | 31 / 32 | yes / yes | 10 | 29 / 31 of 31 |
| (iii) multi | 44 | 30 / 31 | yes / yes | 10 | 29 / 30 of 30 |
| hybrid | 45 | 31 / 32 | yes / yes | 10 | 29 / 31 of 31 |
| hybrid-concat | 45 | 31 / 32 | yes / yes | 10 | 30 / 31 of 31 |

- The 95 % clusters are the 31 GTDB species. Ten species were binned in two or three samples (six in two, four in all three), eight of them in both A and B; the 45 MAGs hold 31 species.
- 99 % splits one cluster: *Clostridium_AQ innocuum*, whose sample B MAG separates from those of A and C.
- **Representatives.** dRep scores completeness − 5 × contamination + contamination × (strain heterogeneity / 100) + 0.5 log10 N50 + (centrality − S_ani) (its defaults; the strain heterogeneity term is zero with CheckM2, which does not report it). galah's default `Parks2020_reduced` also subtracts 5 × contigs / 100. On today's MAGs galah chose a 43-contig *Parabacteroides distasonis* bin (85.4 % / 0.35 %, MIMAG medium) where dRep chose a 271-contig one (93.4 % / 1.70 %, MIMAG high). The other differences were between near-equal bins. With `--quality-formula completeness-5contamination`, galah chose exactly dRep's representative in every cluster of all six variants (at 95 %). galah's `dRep` formula needs CheckM1 strain heterogeneity and stops with CheckM2 input.

Cost (cgroup, 8 CPUs):

| | 45 MAGs | 178 MAGs (today, self, (ii) and (iii) together) |
|---|---|---|
| dRep | 11–16 s, 15–18 CPU-s, 0.2–0.5 GB | 21 s, 41 CPU-s, 0.37 GB |
| galah | ~3 s with container start, under 1 CPU-s (below the 2 s sampler) | 2.2 s, 13.7 CPU-s, 0.23 GB (mean of 10 runs) |

The 178 MAGs gave the same 31 clusters.

### Abundance

CoverM 0.8.0 `genome` (strobealign, `--min-read-percent-identity 97 --min-covered-fraction 0`, mean, trimmed mean, covered fraction, relative abundance) on the 31 representatives, with every sample's reads:

| Run | Wall | CPU-h | Peak | Reads mapped A / B / C / D |
|---|---|---|---|---|
| (ii) cross representatives | 6.8 min | 1.43 | 14.5 GB | 71 / 71 / 50 / 17 % |
| (i) today's representatives (B4 off) | 5.7 min | 1.24 | 14.4 GB | 73 / 72 / 51 / 17 % |

Against the binning depth (length-weighted mean of the jgi contig depths at ≥ 97 % identity, BIN_DEPTH's definition), from the N × N depth tables. Each representative was taken against its own assembly, depth > 1 only:

| Sample | Representatives | CoverM / binning depth, median (IQR) | Spearman |
|---|---|---|---|
| A | 28 | 1.004 (1.00–1.01) | 1.000 |
| B | 24 | 1.005 (1.00–1.01) | 0.994 |
| C | 12 | 1.004 (1.00–1.01) | 0.988 |
| D | 4 | 1.020 (1.01–1.03) | 0.992 |

- With B4 off, the binning has each MAG's depth in its own sample only. Abundance across the run's samples then needs CoverM: N alignments to the dereplicated set, ~70 % of today's alignment CPU here, plus seconds for galah.
- With the recommended B4 (concatenation), the depths are competitive (B4, How the cost grows with N) and CoverM is still needed.
- Only with N × N alignments does BIN_DEPTH already hold every representative's depth in every sample, matching CoverM.

### Images

Checked inside each image (`--version` / `--help` of every command used); the depot.galaxyproject.org sifs of the same builds exist and report the same versions, with gzip present:

| Tool | Image | Inside |
|---|---|---|
| galah | `quay.io/biocontainers/galah:0.6.2--hfa8f182_0` | skani 0.4.0 |
| CoverM | `quay.io/biocontainers/coverm:0.8.0--h750ce8b_0` | strobealign 0.17.0, minimap2 2.31-r1302, samtools 1.23.1 |
| dRep | `quay.io/biocontainers/drep:4.0.3--pyhdfd78af_0` (2.5 GB) | skani 0.3.2, Mash 2.3, MUMmer, Prodigal 2.6.3; `dRep --version` is not an option (version from pip) |

All three ran with a read-only root filesystem (`--read-only`, `/tmp` as tmpfs), the way Singularity/Apptainer mounts the image (design doc Q21). dRep tries to write matplotlib and fontconfig caches under `$HOME`. It warns when it cannot and gives the same clusters. galah 0.6.3 (2026-11-02) is less than two weeks old and was not used.

## Recommendations

**B4, multi-sample binning: optional, off by default, per-sample mode only.** Under co-assembly MetaBAT2 already has one column per sample of the group.

- `--multisample_binning` (boolean, default false). Each sample is aligned once more, to the concatenation of its group's assemblies (contig names prefixed with the sample). The concatenation's jgi depth table is split per assembly, and MetaBAT2 bins each assembly with one depth column per sample of the group. SemiBin2 keeps the pre-trained model on the sample's own BAM. BIN_DEPTH, the abundance tables and every other step that uses the per-sample BAM are unchanged. It needs `--include_binning` with MetaBAT2 among the binners. COMEBin, if selected, could take the same columns; that was not measured here, and Phase 18b measures it before it is enabled (owner).
- `--multisample_group_by <column>` of `--sample_metadata`: samples are cross-mapped only within their group, the way Phase 15's `--coassembly_group_by` forms co-assembly groups. It would reuse the same rules: label sanitising, startup errors for a missing table, column or value, metadata rows for other batches ignored. Without it, the run is one group. A one-sample group is binned as today. Phase 15's helper (`coassemblyGroups()`) is tied to `--coassembly_group_by` and co-assembly mode and would need to take the param as an argument.
- The manual should give the cost: one extra alignment per sample, against a reference that grows with the group. It should warn that large groups need memory to index the concatenation and lower each copy's depth. A size above which to warn was not measured here (N = 4). The owner chose a launch warning by the number of samples in a group; Phase 18b sets the number.
- Not chosen: N × N alignments, as an alternative coverage source, give no better bins for g × N alignments. Self-trained SemiBin2 on N BAMs costs ~20× the binning CPU for the same result. `multi_easy_bin` is one serial task.

**B5, dereplication: optional, off by default, and allowed without B4.** The roadmap's section 4 said "Run-level dereplication only with multi-sample binning (Phase 18)". The measurements say it does not depend on it, and the owner allowed it on its own (2026-11-08).

- `--dereplicate_mags` (boolean, default false). It needs `--include_binning` with Binette, which supplies the CheckM2 reports. galah 0.6.2 clusters the refined MAGs of all the run's assemblies with `--quality-formula completeness-5contamination` and Binette's thresholds (`--bin_min_completeness` / `--bin_max_contamination`, 50 / 10). CoverM 0.8.0 `genome` then maps every sample's reads to the representatives (`--min-read-percent-identity 97`, as jgi). Representatives, clusters and abundance (mean depth, relative abundance, covered fraction) go to the summary as new tables.
- `--dereplicate_ani` (default 95): 95 for species. 99 split one species here.
- Documented as covering this run's samples only. A study processed in batches dereplicates downstream, from the merged results.
- dRep is not needed: it gives the same clusters and representatives, but its image is 2.5 GB and it is ~4–10× slower (11–16 s against ~3 s on one variant's 45 MAGs, 21 s against 2.2 s on 178).

**For Phase 18b:** with the recommended B4, the abundance across samples comes from CoverM, not from BIN_DEPTH. The concatenation's depths serve only MetaBAT2.

## Caveats

- Four samples, one of them too shallow to bin (D). Ten species were shared between samples, most between A and B. A run with more samples of the same environment may gain more from extra coverage columns than this one did; this benchmark does not show that.
- The differences in CheckM2 high-quality bins between the variants with multi-sample coverage (27–30) come down to one or two bins near the 90 % / 5 % line. The consistent signals are the aggregate score (Σ +100–114 over today for every variant but (iii) and self) and the median completeness (+2.3–3.5 points, from the unrounded medians).
- The two diagnostics were added after (ii) and (iii). They reuse the bench's MetaBAT2 bins and 11a's SemiBin2 bins, with Binette run per assembly the same way.
- The extrapolation to 20 and 50 samples assumes these samples' mean read count and alignment cost. Index memory and read splitting in the concatenation are linear guesses beyond N = 4.
- SemiBin2's pre-trained model was `global`, as in Phase 11a; the pipeline's per-sample default is `human_gut`.
- "Today" is the Phase 11a reference configuration: SemiBin2 and MetaBAT2, refined by Binette. The pipeline's default `--binners` is `semibin` alone, which runs neither MetaBAT2 nor Binette. Multi-sample binning needs MetaBAT2 and dereplication needs Binette, so both need `--binners` with at least two binners (MetaBAT2 among them for B4); the gains above are against that configuration, not against the default.
- Wall times are from a shared host.

## Appendix: COMEBin with the multi-sample columns (Phase 18b)

**Question.** COMEBin reads the BAMs itself (one coverage column per BAM). With `--multisample_binning` the pipeline already aligns every sample of a group to the concatenation, so COMEBin could take the same N columns at no extra alignment cost. Does that improve the refined bins?

**Setup.** Samples C and D (the 11a contigs), the bench's concatenation BAMs from hybrid-concat. For each assembly the four concatenation BAMs were reheadered to that assembly's contig names (`samtools reheader`; COMEBin computes coverage with `bedtools genomecov` and ignores the other assemblies' contigs). COMEBin 1.1.0 (`comebin:1.1.0--hdfd78af_2`, the pipeline's COMEBIN module and its `config/modules.config` settings, 8 CPUs) ran in two variants:
- **own**: the sample's own BAM, as the pipeline runs it today;
- **multi**: the four reheadered BAMs.

Each was refined by Binette per assembly from SemiBin2 pre-trained `global` (11a), MetaBAT2 with the concatenation's columns (hybrid-concat) and the COMEBin variant. The baseline is hybrid-concat, which has no COMEBin. COMEBin was run once per variant; its run-to-run spread was not measured.

**Refined bins** (CheckM2 from Binette's reports; Σ = Σ(completeness − 5 × contamination)):

| Assembly | Refinement | Bins | > 90 % / < 5 % | ≥ 50 % / < 10 % | Σ | Mbp |
|---|---|---|---|---|---|---|
| C | hybrid-concat (no COMEBin) | 6 | 2 | 6 | 412.6 | 19.8 |
| C | + COMEBin own | 7 | 2 | 7 | 479.2 | 26.6 |
| C | + COMEBin multi | 7 | 2 | 7 | 497.0 | 27.1 |
| D | hybrid-concat (no COMEBin) | 0 | 0 | 0 | 0 | 0 |
| D | + COMEBin own | 1 | 0 | 1 | 61.7 | 0.7 |
| D | + COMEBin multi | 1 | 0 | 1 | 78.7 | 1.0 |

**COMEBin's own bins**, before refinement (Σ over every bin):

| Assembly | Variant | Bins | > 90 % / < 5 % | ≥ 50 % / < 10 % | Σ | Mbp |
|---|---|---|---|---|---|---|
| C | own | 16 | 1 | 7 | 525.1 | 35.9 |
| C | multi | 15 | 1 | 7 | 602.4 | 36.1 |
| D | own | 6 | 0 | 0 | 115.3 | 1.9 |
| D | multi | 8 | 0 | 1 | 141.4 | 3.2 |

**Cost** (Nextflow trace; CPU-hours = realtime × %CPU; host at load 28–43):

| Variant | C | D | CPU-h (both) | Peak RSS |
|---|---|---|---|---|
| own | 3 h 2 min, 22.3 CPU-h | 26 min, 3.1 CPU-h | 25.3 | 7.3 GB |
| multi | 3 h 5 min, 21.9 CPU-h | 31 min, 3.1 CPU-h | 25.0 | 7.2 GB |

Reheadering took seconds per BAM. Binette took 2–5 min per assembly (≤ 6.6 GB) in either variant.

**Reading.**
- Adding COMEBin to the refinement helps either way: on C, Σ +66.6 (own) and +84.4 (multi) over hybrid-concat, one more bin and +7 Mbp; on D, the first bin of that sample.
- The multi columns add less on top of that: Σ +17.8 on C (+3.7 %) and +17.0 on D, with the same bin counts and no change in high-quality bins. The raw COMEBin bins move the same way (Σ +77 on C, +26 on D).
- They cost nothing extra: COMEBin's time and memory are the same with four BAMs as with one, and the alignments are the ones `--multisample_binning` already makes.
- The gain is small, on two assemblies, from one run per variant. It could be inside COMEBin's run-to-run spread, which was not measured.

**Decision (owner, 2026-11-09):** COMEBin takes the group's concatenation BAMs under `--multisample_binning`.

**Implementation.** Rather than one reheadered copy of every BAM per assembly (N² BAM copies on disk), the pipeline gives COMEBin the assembly's contigs under their `<sample>:` names and the unmodified concatenation BAMs; the bins get the plain names back. COMEBin's coverage is the same either way. Only the contig names change.

### Confirmation of the pipeline's wiring

The `prefix` variant is the wiring above, run with the pipeline's COMEBIN module on the same inputs, then refined by Binette as before.

- **D twice.** Two runs gave byte-identical bins: COMEBin is deterministic for a given input (its seed is fixed).
- **But it is sensitive to contig names.** The prefix and the reheader variants give COMEBin the same coverage and differ only in the contig names, yet their bins differ, as two seeds would.

| Assembly | Refinement | Bins | > 90 % / < 5 % | ≥ 50 % / < 10 % | Σ | Mbp |
|---|---|---|---|---|---|---|
| C | + COMEBin own | 7 | 2 | 7 | 479.2 | 26.6 |
| C | + COMEBin multi (reheader) | 7 | 2 | 7 | 497.0 | 27.1 |
| C | + COMEBin multi (prefix) | 6 | 2 | 6 | 434.8 | 21.3 |
| D | + COMEBin own | 1 | 0 | 1 | 61.7 | 0.7 |
| D | + COMEBin multi (reheader) | 1 | 0 | 1 | 78.7 | 1.0 |
| D | + COMEBin multi (prefix) | 1 | 0 | 1 | 75.8 | 1.0 |

COMEBin's own bins with the prefix variant: C 12 bins, Σ 593.8, 36.0 Mbp; D 8 bins, Σ 146.8, 3.5 Mbp. Cost: C 2 h 25 min, 16.9 CPU-h; D 26 min, 2.8 CPU-h; peak RSS 7.2 GB.

**Reading.**
- The raw COMEBin gain holds in all four multi-sample runs: Σ +69 to +77 on C and +26 to +31 on D over the own BAM.
- The refined gain does not. On C the prefix variant loses an 80 %-complete bin that the other two refinements keep, and ends 44.4 below the own BAM. The spread between two runs with identical coverage (Σ 62.2 on C) is larger than the gain measured between own and multi (+17.8).
- One run per variant cannot separate the effect of the extra columns on the refined bins from COMEBin's sensitivity to the contig names.
