# Co-assembly: read normalization, MEGAHIT settings and MetaHipMer2

Benchmark run 2026-10-23 to 2026-10-24 (roadmap Phase 5b, items A1–A3). It asks three questions about the co-assembly step:

- **A1:** does digital normalization of the assembler input (BBNorm) make the co-assembly cheaper without losing bins?
- **A2:** do MEGAHIT settings (`--min-count 3`, a shorter k-mer list) do the same?
- **A3:** is MetaHipMer2 2.2.2.0 a useful alternative assembler on this kind of machine?

No pipeline code changed. Every variant ran the pipeline's own MEGAHIT command, and its contigs went through the pipeline's own modules from the 1 kb filter to CheckM2.

**Outcome (recommendations, for the owner's decision):**

- **A1:** add normalization as an opt-in step, default off, with one mode (`normalize`, no error correction). At this depth it cost 11.7 CPU-hours to save 2.1 in the assembly. Error correction gave the same assembly for 2.3 more CPU-hours.
- **A2:** expose `--megahit_min_count` (default 2) and `--megahit_k_list`, and make the default k-mer list `21,41,61,81,101,121,141`: 12 % less assembly CPU time for 0.7 % less sequence, and no fewer bins. Keep `--min-count 3` as a knob for co-assemblies that run out of memory or time, not as the default.
- **A3:** do not add MetaHipMer2 as an assembler option for single-node runs.
- **MEGAHIT memory:** request memory from the size of the input for co-assembly: `max(72 GB, 6 GB per GB of gzipped pooled reads) × attempt`.

## Summary

| Variant | Assembly CPU-h | + normalization | Contigs ≥ 1 kb (Mbp) | N50 | Map-back | Refined bins | MIMAG high / medium |
|---|---|---|---|---|---|---|---|
| baseline 1 | 12.50 | | 27,864 (136.05) | 17,726 | 87.92 % | 23 | 11 / 12 |
| baseline 2 | 12.77 | | 27,864 (136.05) | 17,726 | 87.92 % | 24 | 10 / 14 |
| mc3 | 11.21 | | 25,360 (128.83) | 20,292 | 87.28 % | 23 | 12 / 11 |
| klist | 11.13 | | 27,633 (135.05) | 18,500 | 87.71 % | 24 | 15 / 9 |
| mc3klist | 10.03 | | 25,144 (127.93) | 21,344 | 87.09 % | 23 | 14 / 9 |
| norm | 10.51 | 11.65 | 26,298 (132.45) | 19,417 | 87.45 % | 25 | 12 / 12 |
| normecc | 10.52 | 13.97 | 26,290 (132.45) | 19,413 | 87.50 % | 24 | 11 / 11 |
| normklist | 9.70 | 11.65 | 26,059 (131.75) | 20,180 | 87.38 % | TBD | TBD |
| full depth (reference) | 16.9 | | 34,888 (182.41) | 19,512 | 91.24 % | 35 | 17 / 18 |

- **No setting changed the bins as much as depth did.** Every subsampled variant gave 23–25 refined bins and 22–24 MIMAG high plus medium bins. The full-depth assembly gave 35 and 35.
- **The shorter k-mer list is the only change that saves time at almost no cost.** Assembly CPU time fell 12 %. Total length fell 0.7 % and contigs ≥ 2.5 kb fell 2.2 %. Both k-list variants moved bins from MIMAG medium to high (14–15 high vs 10–12 elsewhere) without changing their total.
- **`--min-count 3` trades sequence for time.** CPU time fell 11 %, but it lost 5 % of the sequence in contigs ≥ 1 kb, 9 % of those contigs, and 0.6 points of map-back. Its bins are within the baseline spread.
- **Normalization costs more than it saves at this depth.** It kept 84 % of the pairs and saved 2.1 CPU-hours in MEGAHIT, but took 11.7 CPU-hours itself.
- **No MEGAHIT setting lowered MEGAHIT's peak memory.** Peak anonymous memory was 4.26 GB in every run without normalization, and 3.62 GB with it. `--min-count 3` lowered only the graph build's minimum (2.4 → 2.0 GB).
- **The binning chain costs more than the assembly.** Per assembly, SemiBin2 took 9–10 h at 2–3 of its 8 CPUs, 16–24 CPU-hours: more than MEGAHIT's 12.6. MetaWRAP took ~3 h, 10–11 CPU-hours, at 68.7 GB.

## Setup

**Reads:** a 5-sample co-assembly set (human gut, 2×151 bp NovaSeq X, clean reads after the pipeline's QC and host removal), 38.9 M pairs, ~10.9 Gbp, 4.79 GB gzipped. One of the five had been capped at 20 M pairs when the set was made, so this is a small co-assembly, not five full-depth samples.

**Subsample:** a MEGAHIT co-assembly of the full set took 1 h 18 min at 16 CPUs on an earlier run, and the machine was busy, so every assembly used the same 50 % subsample of each sample, pooled in the pipeline's order: 19.44 M pairs, 5.44 Gbp, 2.50 GB gzipped. It was made with BBTools 40.05 `ReformatReads` (`samplerate=0.5 sampleseed=11`, pairs kept together). BBTools 40.05's `reformat.sh` now runs a streaming class that ignores `samplerate` and writes every read (it also miscounts its input); the pipeline does not use `reformat.sh`.

**Full reads for everything after the assembly:** map-back and binning used the full 38.9 M clean pairs of the five samples, as the pipeline does, so the variants differ only in their contigs.

**Read reduction understates its gains here.** Normalization removes reads above the target depth. The deeper the input, the more it removes, so a small set, halved again, shows the least it can do. A deep co-assembly would gain more. The same goes for `--min-count 3`: at half depth, more true k-mers fall below a count of 3 (so its loss of sequence is overstated) and there are fewer error k-mers to remove (so its saving is understated).

**Machine and runs:**

- 4 × Xeon E7-4820 v4 (40 cores, 80 threads, AVX2), 2 TB RAM; Docker; Nextflow 26.04.4 for the binning chain.
- Every assembly ran in its own container at `process_high`'s first attempt, 16 CPUs and 72 GB, through `bin/megahit_run.sh` as the MEGAHIT module calls it (so `-m` = 90 % of 72 GB).
- Other users kept the 1-minute load at 104–232 on 80 threads the whole time. Each 16-CPU container got 4–5 threads on average, so wall times are about 3× what an idle machine gives. Several variants ran at the same time.
- Time and memory come from the container's cgroup (CPU-seconds, peak memory, peak anonymous memory), sampled every 2 s. "Peak anon" is the process memory; the cgroup peak also counts page cache.

**CPU-hours compare better than wall time, but not perfectly.** Under that load, threads share physical cores, which inflates CPU-seconds. The two baselines ran at the same time and differ by 2 %. The other variants ran at other times and loads (the baselines from 104 to 157, mc3 and klist from 130 to 201, mc3klist and the normalized ones from 185 to 229). So the 11–12 % savings of mc3 and klist are probably real but near the edge of what this machine can show. The normalization step's cost (+90 % of the assembly) is far outside that.

**Binning chain, per assembly:** BBMAP (1 kb filter) → bwa-mem3 index → bwa-mem3 mapping of the pooled full reads with the pipeline's filter (`--bwamem3_min_identity 0.95`, `--bwamem3_min_aligned_length 50`) → depth → MetaBAT2 and SemiBin2 (`train_self`, as in co-assembly mode) → MetaWRAP → CheckM2. The modules and `config/modules.config` came from dev `7d42827`. Map-back is the share of primary reads mapped in `samtools flagstat` of that BAM; reads failing the filter stay in it as unmapped.

**Bin classes:** the pipeline's high (completeness ≥ 90 %, contamination ≤ 5 %) and mid (≥ 50 %, ≤ 5 %) classes, and MIMAG high (> 90 %, < 5 %) and medium (≥ 50 %, < 10 %, not high). tRNA and rRNA were not checked, so MIMAG high here is completeness and contamination only.

**Variants:**

| Variant | Reads into MEGAHIT | MEGAHIT arguments |
|---|---|---|
| baseline 1, baseline 2 | subsample | none (defaults: `--min-count 2`, `--k-list 21,29,39,59,79,99,119,141`) |
| norm | subsample after `bbnorm.sh target=100 min=2 prefilter=t` | none |
| normecc | subsample after `bbnorm.sh target=100 min=2 prefilter=t ecc=t` | none |
| mc3 | subsample | `--min-count 3` |
| klist | subsample | `--k-list 21,41,61,81,101,121,141` |
| mc3klist | subsample | both |
| normklist | the norm reads | `--k-list 21,41,61,81,101,121,141` |
| MetaHipMer2 | subsample, uncompressed, `/1` `/2` read names | defaults |
| full depth (reference) | the full 38.9 M pairs, MEGAHIT from an earlier pipeline run | none |

normklist pairs the better normalization with the better MEGAHIT setting. It was chosen from the assembly numbers before the bins were ready: normecc gave the same assembly as norm for 2.3 more CPU-hours, and mc3 lost 5 % of the sequence in contigs ≥ 1 kb.

BBNorm ran with the BBTools image the pipeline already uses (`bbmap:40.05--h09cc210_0`), 16 threads, `-Xmx60g` in a 72 GB container.

## The noise floor: two baseline runs

MEGAHIT 1.2.9 is not reproducible run to run. The two baselines produced the same number of contigs with the same length distribution (137,354 contigs, 190,414,477 bp), but 116 contigs (0.24 % of bases) differ in sequence; in contigs ≥ 1 kb, 28 contigs and 0.31 % of bases. After the 1 kb filter, contig counts, total length, N50 and map-back are identical. CPU time differed by 2 % (12.50 vs 12.77 CPU-h).

The bins moved more than the contigs. From the same contig statistics, binning gave:

- 23 vs 24 refined bins;
- 11 vs 10 pipeline-high and 10 vs 12 pipeline-mid;
- 11 vs 10 MIMAG high and 12 vs 14 MIMAG medium;
- a median completeness of 90.2 vs 87.5 %.

SemiBin2 is not seeded and gave 53 vs 57 bins; MetaBAT2 gave 50 vs 49. So the noise floor is zero for the filtered contig statistics and 1–2 bins per class. A variant counts as different below only where it falls outside that.

## Assembly cost and contigs

| Variant | Wall (min) | CPU-h | %CPU | Peak anon (GB) | Cgroup peak (GB) | MEGAHIT minimum (GB) | Total (Mbp) | N50 | Longest | ≥ 1 kb | ≥ 2.5 kb |
|---|---|---|---|---|---|---|---|---|---|---|---|
| baseline 1 | 148.8 | 12.50 | 504 | 4.26 | 7.84 | 2.4 | 136.05 | 17,726 | 547,882 | 27,864 | 8,568 |
| baseline 2 | 142.0 | 12.77 | 539 | 4.26 | 7.82 | 2.4 | 136.05 | 17,726 | 547,882 | 27,864 | 8,568 |
| mc3 | 153.6 | 11.21 | 438 | 4.26 | 7.20 | 2.0 | 128.83 | 20,292 | 547,882 | 25,360 | 7,826 |
| klist | 164.8 | 11.13 | 405 | 4.26 | 7.83 | 2.4 | 135.05 | 18,500 | 639,717 | 27,633 | 8,379 |
| mc3klist | 139.7 | 10.03 | 431 | 4.26 | 7.20 | 2.0 | 127.93 | 21,344 | 639,717 | 25,144 | 7,631 |
| norm | 145.5 | 10.51 | 433 | 3.62 | 6.58 | 2.0 | 132.45 | 19,417 | 547,882 | 26,298 | 8,187 |
| normecc | 142.5 | 10.52 | 443 | 3.63 | 6.57 | 1.9 | 132.45 | 19,413 | 547,882 | 26,290 | 8,191 |
| normklist | 151.6 | 9.70 | 384 | 3.62 | 6.48 | 2.0 | 131.75 | 20,180 | 639,717 | 26,059 | 8,029 |
| full depth | 77.9 | 16.9 | 1305 | 8.5 (peak RSS) | | | 182.41 | 19,512 | 532,103 | 34,888 | 11,007 |

Total length, N50 and contig counts are after the pipeline's 1 kb filter. "MEGAHIT minimum" is the largest "Minimum memory required" in MEGAHIT's log (the graph build). The full-depth row comes from the earlier pipeline run's trace, on a less loaded machine.

Against the baselines (identical after the filter, so any change is outside the noise floor):

| Variant | CPU-h | Total length | ≥ 1 kb | ≥ 2.5 kb | Map-back |
|---|---|---|---|---|---|
| mc3 | −11 % | −5.3 % | −9.0 % | −8.7 % | −0.64 points |
| klist | −12 % | −0.7 % | −0.8 % | −2.2 % | −0.21 |
| mc3klist | −21 % | −6.0 % | −9.8 % | −10.9 % | −0.83 |
| norm | −17 % (+75 % with BBNorm) | −2.6 % | −5.6 % | −4.4 % | −0.47 |
| normecc | −17 % (+94 %) | −2.6 % | −5.7 % | −4.4 % | −0.42 |
| normklist | −23 % (+69 %) | −3.2 % | −6.5 % | −6.3 % | −0.54 |

Each setting raises N50 because it drops short contigs, not because the long ones grow. Halving the reads cost MEGAHIT only 26 % of its CPU time (16.9 → 12.6 CPU-h, on different loads), and 25 % of the sequence in contigs ≥ 1 kb.

**A1, the normalization step itself.**

| | Wall (min) | CPU-h | %CPU | Peak anon (GB) | Pairs kept | Bases kept |
|---|---|---|---|---|---|---|
| `bbnorm.sh target=100 min=2 prefilter=t` | 164.3 | 11.65 | 425 | 61.3 | 16,350,403 of 19,436,173 (84.1 %) | 84.9 % |
| same with `ecc=t` | 191.3 | 13.97 | 438 | 61.2 | 16,347,998 (84.1 %) | 84.9 % |

- BBNorm counts k-mers in two passes and takes as much of the Java heap as it is given, so 61 GB is the `-Xmx60g` allocation, not a measured need (its prefilter took 15.3 GB). A module should size `-Xmx` from the task's memory. Fewer passes (`passes=1`) were not tried.
- With error correction, 2,405 fewer pairs were kept and the filtered assembly was the same within MEGAHIT's own noise: the same total length, 8 fewer contigs ≥ 1 kb.
- At this depth, the step removed only 16 % of the pairs, so MEGAHIT saved 2.1 CPU-hours for the 11.7 the step took.

**Memory.** Neither MEGAHIT setting changed the peak (4.26 GB anon in all five runs without normalization). The peak is in the k-mer counting of the first k, which scales with the input; only read reduction lowered it (3.62 GB). `--min-count 3` lowered MEGAHIT's graph-build minimum from 2.4 to 2.0 GB.

## Map-back and bins

| Variant | Map-back | MetaBAT2 | SemiBin2 | Refined | Pipeline high / mid | MIMAG high / medium | Median completeness / contamination |
|---|---|---|---|---|---|---|---|
| baseline 1 | 87.92 % | 50 | 53 | 23 | 11 / 10 | 11 / 12 | 90.24 / 0.70 |
| baseline 2 | 87.92 % | 49 | 57 | 24 | 10 / 12 | 10 / 14 | 87.49 / 0.70 |
| mc3 | 87.28 % | 51 | 54 | 23 | 12 / 9 | 12 / 11 | 94.35 / 0.59 |
| klist | 87.71 % | 52 | 54 | 24 | 15 / 7 | 15 / 9 | 96.03 / 0.70 |
| mc3klist | 87.09 % | 49 | 49 | 23 | 14 / 7 | 14 / 9 | 95.17 / 0.70 |
| norm | 87.45 % | 50 | 52 | 25 | 12 / 9 | 12 / 12 | 94.75 / 0.56 |
| normecc | 87.50 % | 53 | 55 | 24 | 11 / 10 | 11 / 11 | 92.22 / 0.56 |
| normklist | 87.38 % | TBD | TBD | TBD | TBD | TBD | TBD |
| full depth | 91.24 % | 66 | 74 | 35 | 17 / 16 | 17 / 18 | 91.39 / 0.60 |

Map-back and bins are against the full reads. Completeness and contamination are CheckM2's, for MetaWRAP's refined bins.

- **Every variant is within the baseline spread for refined bins** (23–25 vs 23–24) and for MIMAG high plus medium (22–24 vs 23–24).
- **The two k-list variants are outside it for MIMAG high**, in their favour: 15 and 14 vs 10–11 for the baselines and 11–12 for every other variant. Their total stayed the same, so bins moved from medium to high rather than new bins appearing. Each variant ran once, so this is a lead to check on the Phase 5c run, not a measured gain.
- **Map-back fell with every setting**, by 0.2 (klist) to 0.8 points (mc3klist), against 0 between the baselines. That is roughly the share of sequence each setting dropped.
- **Depth is what matters for the bins.** The full-depth assembly gave 35 refined bins and 17 MIMAG high, and mapped 91.2 % of the reads.

**Binning cost per assembly** (from the Nextflow traces; same load caveat):

| Step | Wall | Cores used | CPU-h | Peak RSS |
|---|---|---|---|---|
| bwa-mem3 mapping (16 CPUs) | 32–44 min | 2.7–3.6 | ~2 | 17.4 GB |
| SemiBin2 `train_self` (8 CPUs) | 9 h 03 min – 10 h 15 min | 1.7–2.7 | 16–24 | 21 GB |
| MetaWRAP | 2 h 47 min – 3 h 31 min | 3.1–3.7 | 10–11 | 68.7 GB |
| CheckM2 (all three bin sets) | 21–41 min | 5–10 | ~3.5 | 14.2 GB |

SemiBin2's self-training is the largest single cost of a co-assembly, larger than MEGAHIT's. It used a quarter to a third of its 8 CPUs. Phase 7a (the binning benchmark) looks at this.

## A3: MetaHipMer2 2.2.2.0

**How it was installed.** There is no bioconda package. The developer's Docker Hub image (`robegan21/mhm2`) stops at v2.2.0.0-155 (June 2024), and the upstream Dockerfile clones the default branch rather than a tag, so the image was built locally from the `v2.2.2.0` git tag (commit `a7f17d1`, 2025-10-02) with its submodules: UPC++ 2023.9.0, GASNet 2024.5.0 (smp conduit), Debian bookworm-slim. CMake reported "Building for CPU only" (no CUDA or HIP). Its CPU-architecture tests failed on this Broadwell CPU, so the build has no `-march` flags (generic x86-64). `mhm2 -h` reports `MHM2 version 2.2.2.0- with upcxx-utils 0.4.0` (the git hash is empty because `.git` was not copied into the build context). It runs on this CPU without a GPU.

**What it needs that MEGAHIT does not:**

- **Uncompressed FASTQ.** It rejects gzip input. Decompressing the subsample took 1.3 min and 11.8 GB of disk; the full reads of a large co-assembly would need tens to hundreds of GB of scratch space.
- **Paired read names ending in `/1` and `/2`** (or an Illumina `1:N:…` comment). With two paired files, its reader finds the start of a record only from those markers. The pipeline's clean reads have neither, and the first run stopped at once with "Missing record info" and "Could not find the partition boundary". Rewriting the names took 4.6 min.
- **A memory setting from the task.** It sizes its k-mer tables from the free memory in `/proc/meminfo`, which inside a container is the host's (987 GB here), not the container's limit. It estimated 3.34 GB per rank of tables for k=21 and, seeing room, scaled them up 3× (the code's maximum). At 72 GB the container killed it after 4.8 min. In a pipeline it would need `--kcount-fraction-mem` set from the task's memory, and it would still not adapt to a scheduler's limit.
- **A shared-memory size.** Without `--procs` it starts one rank per core (80) and hangs at Docker's default `/dev/shm`; the runs here used `--procs 16`, `--shm-size` equal to the memory limit and `GASNET_PHYSMEM_MAX`.

**Cost.** With a 256 GB limit (and the same 1.84 GB shared heap per rank) it ran for 2 h 04 min and was stopped by the owner's decision while still at k=21, the first of its five k values (21, 31, 47, 63, 95, then scaffolding). By then it had used 12.0 CPU-hours, about what the whole MEGAHIT baseline used (12.5), and peaked at 114 GB (MEGAHIT: 4.3 GB anon, 7.8 GB cgroup). Its k=21 k-mer analysis alone took 40 min. No contigs, map-back or bins exist for it.

**Why stop it.** MetaHipMer2's speed comes from spreading one assembly over many nodes of a cluster, which MEGAHIT cannot do. On one node it was slower in CPU time, not only in wall time, so more cores would not close the gap, and it needed ~25× MEGAHIT's memory.

## Recommendations and defaults

**A1: opt-in, default `none`.** `--assembly_read_reduction none|normalize`, with `normalize` = `bbnorm.sh target=100 min=2 prefilter=t`, on the assembler input only. Leave out the error-correction mode: it gave the same assembly for 2.3 more CPU-hours. Size BBNorm's `-Xmx` from the task's memory (e.g. 85 %); 16 CPUs and 72 GB covered this set with room to spare. The case for it is a deep co-assembly, where it removes far more than 16 % of the reads. That case was not measured here: if the owner wants a default other than `none`, a deep real set should be measured first.

**A2: shorter k-mer list by default, min-count unchanged.**

- `--megahit_k_list` default `21,41,61,81,101,121,141`: 12 % less assembly CPU time, 0.7 % less sequence, 2.2 % fewer contigs ≥ 2.5 kb, 0.2 points less map-back, and no fewer bins (more MIMAG high in its one run). If the owner prefers outputs to stay MEGAHIT's own, the default can be empty (MEGAHIT's list) and the shorter list documented. The saving is ~1.5 CPU-hours on a 12.6 CPU-hour assembly, against a ~30 CPU-hour binning chain.
- `--megahit_min_count` default 2. Document 3 as the setting for co-assemblies that run out of memory in the graph build or run too long: it cut the graph-build minimum by 17 % and CPU time by 11 %, for 5 % less sequence at this (understated) depth.

**A3: do not adopt** for single-node runs (section above). If it is ever reconsidered for a cluster, the module would need uncompressed `/1` `/2` input, `--kcount-fraction-mem` from the task's memory, `--procs` and a shared-memory size.

**MEGAHIT's memory request.** Today MEGAHIT gets `process_high`'s 72 GB × attempt (`-m` = 90 % of it) for per-sample assembly and co-assembly alike. Known points:

| Input | Gzipped reads | Graph-build minimum | Measured peak |
|---|---|---|---|
| this subsample (19.4 M pairs) | 2.50 GB | 2.4 GB | 4.26 GB anon (7.8 GB cgroup) |
| this full set (38.9 M pairs) | 4.79 GB | not kept | 8.5 GB peak RSS |
| a 404 M-read co-assembly | ~25 GB (estimated, not measured) | ~86 GB at k=21 | not kept; killed at 72 GB on the first attempt after 2 h 19 min |

- Memory grows faster than the input: about 1.8 GB of peak per GB of gzipped reads on the small set, and at least 3.4 GB per GB (the graph-build minimum alone) on the large one. A co-assembly of many samples adds distinct k-mers with every sample.
- On this set, the peak was 1.8× the graph-build minimum. If that ratio holds, the 404 M-read co-assembly needs ~150 GB. The first attempt fails at 72 GB, and the second at 144 GB is at the edge.
- **Recommended:** for co-assembly, `memory = max(72 GB, 6 GB × gzipped pooled reads in GB) × attempt`, the same pattern as BWAMEM3_INDEX's size-based request (`resourceLimits` still caps it at `--max_memory`). It gives 72 GB here and ~150 GB for the 404 M-read co-assembly, and saves its 2 h 19 min failed first attempt. Per-sample assembly can keep 72 GB. A request does not make MEGAHIT faster; over-requesting only holds memory other tasks could use.

## Caveats

- One co-assembly set, subsampled to 50 %; one run per variant except the two baselines.
- The machine was heavily loaded by other users throughout; wall times are inflated, and CPU-hours carry some contention too (see Setup).
- MEGAHIT is not reproducible run to run, and SemiBin2 is not seeded: a difference smaller than the gap between the two baselines is noise.
- The memory rule rests on two measured small sets and one large set whose input size is estimated.
- MetaHipMer2 was not run to the end.
