# COMEBin under co-assembly: pooled BAM or one BAM per sample

Benchmark run 2026-11-06 (roadmap Phase 15, B3). Under co-assembly, COMEBin gets the pooled BAM: the samtools merge of the per-sample alignments. Every sample is already aligned to its co-assembly on its own (Phase 7c), and MetaBAT2 and the MAG depth already read those per-sample BAMs. COMEBin accepts several BAMs and reads one coverage column from each. This benchmark asks whether giving it the per-sample BAMs changes its bins, or the bins Binette keeps.

**Decision (owner, 2026-11-06): switch.** Under co-assembly COMEBin now gets the per-sample BAMs, in the same PR (roadmap Phase 15); the benchmark itself changed no pipeline code.

**Outcome:**

- **Same cost.** COMEBin took 20.0 CPU-hours with the 3 per-sample BAMs and 21.9 with the pooled BAM, about 0.5 CPU-hours per Mbp of contigs. Peak RSS was 7.1 vs 7.2 GB.
- **Slightly better.** COMEBin's own passing bins:
  - per-sample BAMs: 9, of which 3 MIMAG high, 29.3 Mbp;
  - pooled BAM: 8, of which 2 high, 27.4 Mbp.

  After Binette with SemiBin2 + MetaBAT2 + COMEBin:
  - per-sample BAMs: 9 MAGs, 4 high, 30.1 Mbp;
  - pooled BAM: 8 MAGs, 4 high, 27.7 Mbp.
- **The gain is one borderline MAG.** It is medium quality (55.3 % complete, 7.2 % contamination, 1.97 Mbp) and shares 1 % of its sequence with any MAG of the pooled run. The other 8 MAGs are the same genomes in both runs: 79–100 % of each MAG's sequence is in one MAG of the other run.
- **COMEBin is still worth having on this group.** SemiBin2 + MetaBAT2 alone gave 6 MAGs (2 high, 20.1 Mbp), the same as the grouped pipeline run. Adding COMEBin gave 2 more high-quality MAGs and 2–3 more MAGs in all, whichever BAMs it read.
- **Recommendation:** give COMEBin the group's per-sample BAMs under co-assembly. This is the same input MetaBAT2 and the MAG depth already use, at no extra cost. The measured gain is small, and it comes from one group of 3 samples with one seed, so the main reason is consistency rather than the gain.

## Summary

| Binner or set | Bins | MIMAG high | Medium | Mbp in passing bins | CPU-h | Wall | Peak RSS |
|---|---|---|---|---|---|---|---|
| SemiBin2 (pre-trained `global`, pooled BAM) | 15 | 1 | 5 | 19.8 | | | |
| MetaBAT2 (per-sample depths) | 8 | 2 | 2 | 13.8 | | | |
| COMEBin, pooled BAM | 20 | 2 | 6 | 27.4 | 21.9 | 2 h 50 min | 7.2 GB |
| COMEBin, per-sample BAMs | 20 | 3 | 6 | 29.3 | 20.0 | 2 h 38 min | 7.1 GB |
| Binette, SemiBin2 + MetaBAT2 | 6 | 2 | 4 | 20.1 | 0.39 | 4 min 44 s | 5.9 GB |
| Binette, + COMEBin pooled | 8 | 4 | 4 | 27.7 | 0.32 | 3 min 54 s | 7.5 GB |
| Binette, + COMEBin per-sample | 9 | 4 | 5 | 30.1 | 0.41 | 4 min 47 s | 7.5 GB |

The quality classes come from CheckM2, through Binette's reports: Binette runs CheckM2 on every input bin and every final bin.
- **High:** completeness ≥ 90 % and contamination ≤ 5 %.
- **Medium:** completeness ≥ 50 % and contamination ≤ 10 %, and not high.

The "Bins" column counts every bin a binner wrote, passing or not. For the Binette rows it counts the bins Binette kept (completeness ≥ 50 %, contamination ≤ 10 %).

## Setup

- **Data:** one group of the grouped real-data run of Phase 15. The 5-sample test fixture (real reads, subsampled; clean reads) was split 3 + 2 by a metadata column. The benchmark uses the 3-sample group, the one with the lower depth:

  | Co-assembly | Contigs ≥ 1 kb | Mbp | N50 |
  |---|---|---|---|
  | 3-sample group (MEGAHIT on the 3 samples' reads) | 13,330 | 42.8 | 4.6 kb |

  The 5-sample co-assembly (182 Mbp) would have cost an estimated 360–580 CPU-hours for both variants, so the owner chose this group.
- **Inputs:** all fixed and taken from that run (branch `roadmap/p15-coassembly-groups` at `a1e1d6d`, Docker):
  - the filtered contigs;
  - the 3 per-sample bwa-mem3 BAMs and their samtools merge (the pooled BAM);
  - SemiBin2's bins, from the pre-trained `global` model on the pooled BAM, `--random-seed 11`, as the pipeline runs it under co-assembly;
  - MetaBAT2's bins, from the per-sample depths.
- **Harness:** a scratch Nextflow workflow ran the pipeline's own COMEBIN and BINETTE modules from that commit:
  - **COMEBin:** twice on the same contigs, once with the pooled BAM and once with the 3 per-sample BAMs. The module passes `-p ./`, and COMEBin reads every BAM in the task directory: its log said `BAMs=3` for the per-sample run.
  - **Binette:** run three times, on SemiBin2 + MetaBAT2, and on SemiBin2 + MetaBAT2 + each COMEBin variant.
- **Versions:**
  - COMEBin 1.1.0 (`comebin:1.1.0--hdfd78af_2`, the module's `realpath` shim, `-d cpu`, `-s 11`, first-attempt batch settings, temperature 0.15 for N50 < 10 kb);
  - Binette 1.2.1 with CheckM2 1.1.0 and the CheckM2 database `uniref100.KO.1`, completeness ≥ 50 %, contamination ≤ 10 %;
  - SemiBin2 2.6.0 and MetaBAT2 2.18, from the run.
- **Resources:** COMEBin and Binette at 8 CPUs each (the pipeline's COMEBin request), COMEBin with 32 GB.

## Machine and runs

- 4 × Xeon E7-4820 v4 (40 cores, 80 threads), 2 TB RAM; Docker; Nextflow 26.04.4.
- Both COMEBin runs ran at the same time, next to the two Phase 15 pipeline runs; the 1-minute load was 30–55 on 80 threads.
- **Time and memory:** CPU-hours are Nextflow's realtime × %cpu (7.6–7.7 of the 8 cores in both COMEBin runs). Peak RSS is Nextflow's.
- **Lower cost than 7a:** the 0.5 CPU-hours per Mbp here is half the 1.04 of Phase 7a ([binning_2026.md](binning_2026.md), section (f)). 7a ran at a load of 60–180, and its training epochs took about twice as long (92 vs 41–47 s).
- **No failures or retries:** no `FAILED.txt`, and both COMEBin runs finished on attempt 1.

## The refined MAGs, matched

Each MAG of the per-sample run is matched to the MAG of the pooled run that holds most of its sequence. The match column is the share of the per-sample MAG's bases found in that MAG.

| Per-sample run: MAG (completeness / contamination, size) | Pooled run: MAG (completeness / contamination) | Match |
|---|---|---|
| high (99.95 / 0.38, 3.67 Mbp) | high (99.95 / 0.38) | 1.00 |
| high (99.12 / 4.41, 4.17 Mbp) | high (99.12 / 4.41) | 1.00 |
| high (97.07 / 1.33, 2.72 Mbp) | high (96.99 / 3.87) | 0.96 |
| high (94.24 / 2.83, 2.83 Mbp) | high (93.12 / 1.51) | 0.86 |
| medium (99.06 / 6.30, 5.11 Mbp) | medium (98.65 / 6.19) | 0.99 |
| medium (81.60 / 0.59, 1.13 Mbp) | medium (81.62 / 0.62) | 0.98 |
| medium (82.27 / 1.11, 4.71 Mbp) | medium (76.51 / 3.89) | 0.79 |
| medium (72.41 / 4.37, 3.79 Mbp) | medium (84.74 / 5.43) | 0.93 |
| **medium (55.29 / 7.17, 1.97 Mbp)** | none (best 0.01) | 0.01 |

- **Every MAG of the pooled run is in the per-sample run.** Each one shares at least 50 % of its sequence with one MAG there.
- **The two variants trade completeness on two MAGs:**
  - one is 82.3 % complete with per-sample BAMs, 76.5 % with the pooled BAM;
  - the other is 72.4 % with per-sample BAMs, 84.7 % with the pooled BAM.
- **The high-quality MAGs are the same four genomes.** Their contamination is the same or lower with per-sample BAMs, except one (2.83 vs 1.51 %).
- **The MAG only the per-sample run has** passes the medium thresholds by a small margin, and has no counterpart in the SemiBin2 + MetaBAT2 set either.

## Recommendation

Give COMEBin the per-sample BAMs of its co-assembly under co-assembly. It costs no more than the pooled BAM, it gave as many or more passing bins before and after Binette, and it is the input the rest of the binning stage already uses: MetaBAT2's depths and the MAG depth come from the same BAMs. With the change, COMEBIN would take the group's per-sample BAMs, sorted by sample as CALCULATE_DEPTH gets them, instead of the merged BAM.

SemiBin2 keeps the pooled BAM, because its pre-trained model takes one coverage column (binning_2026.md, section (d)). BlobTools keeps it as well.

## Caveats

- **One group, one seed:** the gain (one borderline medium-quality MAG, +1.9 Mbp of 30) is small. Phase 7a saw Binette's classes move for one or two bins near a threshold between CheckM2 runs. Neither COMEBin's run-to-run spread under another seed nor a group with more samples was measured. A group with more samples, or with samples whose genomes differ more, would give COMEBin more columns to separate.
- **This group has little depth:** 27–66 % of each sample's reads map to the co-assembly, and its N50 is 4.6 kb. A deeper group may differ.
- **No GTDB-Tk in the harness:** matches are by contig, so "same genome" means the same contigs, not the same GTDB species.
- **"MIMAG high" here means CheckM2's thresholds only:** the rRNA and tRNA counts the pipeline takes from Bakta were not computed. Without Bakta, the pipeline's summary caps the tier at medium.
