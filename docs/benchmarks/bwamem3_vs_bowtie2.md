# bowtie2 vs bwa-mem3 in BugBuster

Benchmark of `--read_aligner bwamem3` (bwa-mem3 0.14.0) against `bowtie2`, then the default, run 2026-10-11 to 2026-10-13. The goal was to decide whether bwa-mem3 should replace bowtie2 in the pipeline's own alignment steps:

- host/phiX decontamination,
- reads→contigs mapping (contig depth, binning, featureCounts),
- co-assembly bin depth.

Woltka and HUMAnN keep bowtie2 with their own databases; RGI keeps its own aligner (KMA by default).

**Outcome:** bwa-mem3 is now the default. `--read_aligner bowtie2` remains as a fallback for CPUs without AVX2 and for nodes without ~75 GB RAM for the human index build.

## Summary

| Site | bowtie2 | bwa-mem3 | Wall-time speedup | CPU saving |
|---|---|---|---|---|
| Decontamination (5 samples, 52–167M pairs each) | 5h04 | 1h15 | 4.1× | 4.7× |
| Read→contig mapping (per-sample and co-assembly) | — | — | 4.4–5.1× (3.9–6.2× with the filter) | 6.5–7.0× |
| Co-assembly bin depth (35 bins, 5 samples) | 266 min | 6 min (7.7 min with the filter) | 43× (34×) | 31× (29×) |

**Decontamination is more accurate with bwa-mem3:**
- **Spike-in:** a mock community with 943k simulated human pairs added. bowtie2 leaked 0.157% of the human pairs and removed 3.02% of the mock. bwa-mem3 leaked none and removed 0.113%.
- **Real samples:** reads removed only by bwa-mem3 are 94–99% full-length human hits that bowtie2 kept. Reads removed only by bowtie2 are 20–30 bp chance matches.

**Bins, taxonomy and contig depth do not change:**
- Contig depth Spearman ≥ 0.996.
- MetaWRAP keeps 35 co-assembly bins, which pair one-to-one with bowtie2's. HQ/MQ counts shift by a few bins (with the filter: 21 / 13 / 1 below MQ, vs 18 / 17).
- GTDB-Tk species sets are identical.

**Unfiltered bwa-mem3 inflates low-coverage values.** bwa mem's local alignment places short or divergent hits from genomes the assembly lacks. In samples the assembly covers poorly, these inflated gene counts (Spearman 0.74 against bowtie2) and the depth of near-absent bins (~1.9×). The pipeline now filters these alignments: `--bwamem3_min_identity 0.95` and `--bwamem3_min_aligned_length 50`, with failing reads kept as unmapped. A co-assembly re-run with the filter in the pipeline agrees with bowtie2:
- gene-count Spearman 0.93–0.999,
- contig depth 0.998,
- bin depth within 1–2% at ≥ 1×,
- the same 35 GTDB species.

The filter costs 20–40% wall time on large mapping tasks; mapping stays 3.9–6.2× faster than bowtie2.

## Setup

**Samples:**
- 5 host-dominated human gut metagenomes, 2×151 bp NovaSeq X (bowtie2 overall alignment to human 69–99%).
- 1 MiSeq mock community (ERR2984773), with 943k human read pairs simulated with ART (HiSeq X profile) from GRCh38 primary chromosomes.

**Host index:** GRCh38.p14 (with alt loci) plus phiX174.

**Machine and runs:**
- 80 cores, 2 TB RAM; Nextflow 26.04.4, `-profile docker`.
- `process_high` (16 CPUs) for every alignment task.
- The two aligners ran one after the other on the same inputs. The second run used `-resume`, so only aligner-dependent tasks re-ran.
- Mapping times come from re-runs of only the mapping tasks, with no other pipeline tasks running.

**Inputs for assembly and binning:** both arms used the same clean reads. Two setups:
- Per-sample assemblies of 4 samples.
- A co-assembly of all 5, with ~39M pairs after subsampling the largest sample to 20M pairs.

Binning used SemiBin + MetaBAT2, refined by MetaWRAP, then CheckM2 and GTDB-Tk r232. Contig-level functional annotation was on.

**Aligner settings:**
- bowtie2 decontamination: the module defaults, `--local --ma 2 --mp 6,2 --score-min G,15,6 -k 1 -N 1 -L 20 -R 2 -i S,1,0.75 --un-conc-gz`, default `-X 500`.
- bowtie2 contig mapping: end-to-end defaults.
- bwa-mem3 decontamination: `-P -T 30`, with a pair kept only if both mates are unmapped.

## Speed and resources

| Task | bowtie2 | bwa-mem3 |
|---|---|---|
| Host + phiX index build (3.3 Gbp, 8 CPUs) | 66 min, 6 GB | 15.5 min, 57 GB |
| Decontamination: step wall time / CPU-h (5 samples) | 5h04 / 256 | 1h15 / 55 |
| Decontamination: peak RSS per task | 4.2 GB | 15.8 GB |
| Mapping, per-sample assemblies (4 samples): wall / CPU-h | 23.2 min / 10.7 | 4.6 min / 1.65 |
| Mapping, co-assembly (pooled + 5 per-sample BAMs): wall / CPU-h | 38.2 min / 24.5 | 8.6 min / 3.5 |
| Mapping: peak RSS per task | 13.8 GB | 17.2 GB |
| Co-assembly bin depth: wall / CPU-h | 266 min / 52 | 6 min / 1.7 |

**The bin-depth speedup is mostly a change of method.** The bowtie2 module built an index for each bin and mapped every read of the sample against it, 35 times per sample. The bwa-mem3 chain maps each sample once against all bins concatenated. That also changes what the number means: depth is now competitive (each read counts for its best bin), so a bin no longer recruits reads from related genomes.

The bin-depth times come from the main runs, not the dedicated mapping-only re-runs: a ~4.5 h bowtie2 re-run could not have changed a 43× gap. A job from another user (~22 cores) shared the machine during the decontamination runs, so absolute times are pessimistic; the ratios hold.

## Decontamination accuracy

### Spike-in

The mock had 7.75M pairs after fastp, plus 943,056 simulated human pairs.

| | Human pairs leaked | Mock pairs removed |
|---|---|---|
| bowtie2 | 1,481 (0.157%) | 234,056 (3.02%) |
| bwa-mem3 | 0 | 8,769 (0.113%) |

**bwa-mem3 needs `-P` (pairing off).** In paired mode, bwa mem skips the `-T` score floor for the pair it chooses. Random 19–22 bp host hits on overlapping mates then pass as proper pairs, which removed 17.5% of the mock. Simulated reads carry no human population variation, so the spike-in is an easy case for both aligners.

### Real samples

| Sample | Clean pairs, bowtie2 | Clean pairs, bwa-mem3 | Removed by bwa-mem3 only | Removed by bowtie2 only |
|---|---|---|---|---|
| AR-3 | 29,280,649 | 27,702,591 | 2,178,326 | 600,268 |
| BZ-10 | 1,637,616 | 226,925 | 1,451,734 | 41,043 |
| EG-9 | 4,445,971 | 2,935,647 | 1,625,902 | 115,578 |
| KP-10 | 1,356,426 | 1,415,758 | 153,815 | 213,147 |
| NH-10 | 16,392,314 | 14,301,179 | 2,525,866 | 434,731 |

20,000 pairs of each disagreement class were re-aligned to the host with both aligners.

**Removed by bwa-mem3 only: host that bowtie2 kept.**
- 94–99% are full-length host hits (≥ 90% of the read, ≥ 50 bp, ≥ 95% identity).
- bowtie2 never paired them concordantly, so `--un-conc` kept them.
- In 4 samples, 72–82% are FR pairs with 501–1000 bp inserts, just past bowtie2's default `-X 500`.
- In the fifth sample, most are overlapping or dovetailed short-insert pairs (`--no-dovetail`).

**Removed by bowtie2 only: not host.**
- None are full-length host hits.
- They are reads under 50 bp, or 20–30 bp local hits accepted by `--local --score-min G,15,6`.

## Downstream agreement

Gene counts are compared as Spearman correlations between the two arms.

| | Spearman | Median ratio (bwa-mem3 / bowtie2) |
|---|---|---|
| Contig depth, per-sample assemblies | 0.996 | 1.011–1.017 |
| Contig depth, co-assembly | 0.997 | 1.023 |
| Gene counts, per-sample assemblies | 0.985–0.997 | 1.00–1.04 |
| Gene counts, co-assembly, well-covered samples | 0.986–0.992 | 1.00–1.02 |
| Gene counts, co-assembly, poorly covered samples (EG-9, BZ-10, KP-10) | 0.930 / 0.814 / 0.738 | 1.00–1.13 |

Bin counts are given as HQ / MQ (CheckM2: HQ = completeness ≥ 90%, contamination < 5%; MQ = ≥ 50%, < 10%).

| MetaWRAP bins | bowtie2 | bwa-mem3 | GTDB-Tk species |
|---|---|---|---|
| AR-3 | 9 / 8 | 9 / 8 | 17 = 17 |
| EG-9 | 1 / 4 | 1 / 4 | 5 = 5 |
| NH-10 | 11 / 11 | 12 / 10 | 22 = 22 |
| Co-assembly (35 bins each) | 18 / 17 | 17 / 18 | 35 = 35 |

**Mapped reads with < 50 aligned bases:**
- bowtie2: 0.5–3.6%.
- bwa-mem3: 1.0–8.9% (per-sample assemblies) and up to 28% (a poorly covered sample against the co-assembly).

## The alignment filter

In the co-assembly, the bwa-mem3 BAMs were recomputed under different alignment filters. The filtered results were compared with today's unfiltered bowtie2 output.

**Bin depth**, as median ratios of bwa-mem3 to bowtie2 over 34 paired bins × 5 samples:

| bowtie2 bin depth | n | Unfiltered | MAPQ ≥ 1 | ≥ 95% identity, ≥ 50 bp |
|---|---|---|---|---|
| < 0.1× | 48 | 1.93 | 1.84 | 0.64 |
| 0.1–1× | 44 | 1.48 | 1.45 | 0.94 |
| 1–10× | 49 | 1.10 | 1.09 | 0.98 |
| ≥ 10× | 29 | 1.03 | 1.02 | 0.99 |

Below 0.1× the filtered value is lower than bowtie2's. The old per-bin method also over-recruits when the organism is absent, because nothing competes for the reads.

**Gene counts (featureCounts), Spearman:**

| Sample | Both unfiltered | bwa-mem3 filtered vs bowtie2 | Both filtered |
|---|---|---|---|
| AR-3 | 0.986 | 0.998 | 0.999 |
| BZ-10 | 0.814 | 0.931 | 0.993 |
| EG-9 | 0.930 | 0.988 | 0.996 |
| KP-10 | 0.738 | 0.959 | 0.990 |
| NH-10 | 0.992 | 0.999 | 0.999 |

**What the filter results show:**
- MAPQ does not remove these hits: such a hit is often the read's only placement.
- The identity/length filter does remove them.
- The gap left between filtered bwa-mem3 and bowtie2 comes from bowtie2's own divergent end-to-end hits: once both are filtered, they agree at ≥ 0.99.

**What is now in the pipeline:**
- `bin/filter_alignments.sh` runs in `BWAMEM3_SAMTOOLS` (contig depth, binning, featureCounts and co-assembly bin depth). Defaults are `--bwamem3_min_identity 0.95` and `--bwamem3_min_aligned_length 50`.
- A failing read stays in the BAM as unmapped, and `samtools fixmate` updates its mate.
- Reads clipped at a contig end still count when 50 bases align.
- Added after this benchmark: a read shorter than `--bwamem3_min_aligned_length` (a 2×50 library, or reads trimmed below 50 bp) needs 80% of its length aligned instead. With the 2×151 bp reads here it only concerns reads fastp trimmed below 50 bp (fastp keeps reads of 15 bp or more); the benchmark was not re-run with it.

### Validation with the filter in the pipeline

The co-assembly was re-run with the filter in `BWAMEM3_SAMTOOLS`. Everything downstream of mapping re-ran: binning, MetaWRAP, CheckM2, GTDB-Tk, featureCounts and bin depth. Each value below compares that run with today's unfiltered bowtie2.

| | Unfiltered bwa-mem3 | Filtered bwa-mem3 |
|---|---|---|
| Contig depth, Spearman / median ratio | 0.997 / 1.023 | 0.998 / 1.006 |
| Gene counts, Spearman (AR-3, BZ-10, EG-9, KP-10, NH-10) | 0.986, 0.814, 0.930, 0.738, 0.992 | 0.998, 0.931, 0.988, 0.959, 0.999 |
| Bin depth, Spearman / median ratio (33 paired bins × 5 samples) | 0.989 / 1.17 | 0.994 / 0.94 |
| Mapped %, difference from bowtie2 | +0.6 to +18 points | −0.5 to −1.8 points |
| Mapped reads with < 50 aligned bases | 1.0–28% | 0% (bowtie2: 0.5–3.6%) |
| MetaWRAP bins, HQ / MQ / below MQ (bowtie2: 18 / 17 / 0) | 17 / 18 / 0 | 21 / 13 / 1 |
| GTDB-Tk species (bowtie2: 35) | 35, identical | 35, identical |

**Bins:**
- 33 of 35 bins are near-identical to bowtie2's (≥ 0.93 bp Jaccard). Of these, one crossed 90% completeness (87% → 93%) and one at 98% completeness crossed 10% contamination (9.82% → 10.05%).
- 2 bins were refined differently by MetaWRAP (Jaccard 0.89 and 0.78), and both gained completeness: 84% → 93% and 81% → 100%.
- Net bin quality is slightly higher: median completeness is 93.3% against bowtie2's 92.6%.

**Cost:**
- Wall time rises 20–40% on the large mapping tasks: AR-3 3.7 → 5.1 min, pooled co-assembly 8.1 → 9.9 min. CPU time stays the same.
- bowtie2 takes 26.5 and 38.1 min for these tasks.
- Bin depth takes 7.7 min, against 266 min for bowtie2.

## Caveats

- One machine and one run per arm.
- Host-dominated gut samples only; low-host samples were not tested.
- The co-assembly was subsampled to ~39M pairs.
- Two pre-existing pipeline failures were worked around identically in both arms:
  - The eggNOG search ran at 256 GB.
  - SemiBin and MetaWRAP failures on zero-bin samples were ignored.
