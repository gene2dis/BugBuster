# Binning: tool updates, bin refiners and co-assembly binning

Benchmark run 2026-10-25 to 2026-10-26. It asks seven questions about the binning chain (CALCULATE_DEPTH → MetaBAT2 / SemiBin2 / COMEBin → MetaWRAP → CheckM2 → GTDB-Tk):

- **(a) Tool updates:** SemiBin2 2.1.0 → 2.6.0, MetaBAT2 2.15 → 2.18, COMEBin 1.0.4 → 1.1.0, GTDB-Tk 2.7.2 → 2.8.0, MMseqs2 15 → 18. Same bins or better, and do the modules still work?
- **(b) Refinement:** MetaWRAP 1.2 against Binette 1.2.1 and DAS Tool 1.1.7, on the updated binners' bins.
- **(c) GTDB-Tk at `--pplacer_cpus 4`:** how many genomes reach pplacer, and its peak memory.
- **(d) Co-assembly coverage:** binners given the pooled BAM against binners given one BAM per sample (each sample aligned to the co-assembly).
- **(e) CheckM2 on the raw binner sets:** what it costs, and whether the refiner already reports it.
- **(f) COMEBin:** keep or remove. BugBuster runs it on CPU only.
- **(g) SemiBin2 in co-assembly mode:** where `train_self` spends its 7–10 h, whether more CPUs help, and what a pre-trained model gives instead.

No pipeline code changed. Every run used the pipeline's own modules, with the updated images swapped in by config and a few workarounds in the harness only (Setup). Pins and fixes are for Phase 7b.

**Decision (owner, 2026-10-26):** all recommendations below approved, with the GTDB-Tk 2.8.0 biocontainer (skani 0.4.0) and SemiBin2's pre-trained `global` model for co-assembly. Phase 7b applies them, Phase 7c the co-assembly coverage.

**Outcome (recommendations):**

- **Refiner: replace MetaWRAP with Binette.** Same or more bins and MIMAG high bins than MetaWRAP on both sets, every GTDB species MetaWRAP found plus one, at about an eighth of the CPU time and a third of the memory (asm set: 3.0 vs 23.0 CPU-hours, 11.2 vs 34.5 GB). DAS Tool is cheaper still but lost bins and kept one below 50 % completeness.
- **CheckM2 on raw bins: drop it under Binette.** Binette runs CheckM2 on every input bin and every final bin and writes both reports. Its classes matched a separate CheckM2 run on all 176 raw and all 45 refined bins of the asm set; on the co-assembly sets one or two bins near the 90 % line fell in a different class (section (d)). The separate run costs 4.5 CPU-hours and 12 GB on the asm set.
- **Pins:** SemiBin2 2.6.0 (more high-quality bins: 23 vs 18 raw on the asm set, 24 vs 17 in the co-assembly), with `--random-seed`. MetaBAT2 2.18, with its bin headers cut to the contig name. MMseqs2 18 (identical clusters). GTDB-Tk 2.8.0 (identical calls); two images, which differ in skani version (section (a)); the owner chose the biocontainer.
- **COMEBin: keep as an option (not the default, as today), with a cost warning and 8 CPUs.** It took 38 CPU-hours on the sample measured in Phase 5, in both versions: about 370× SemiBin2. After Binette it added one genome, upgraded one bin to high quality and added 39 % more sequence in passing bins on that sample, and gave the only passing bin on the small one. COMEBin 1.1.0 does not run as built: it needs a `realpath` fix, `-d cpu` and its bin headers cut.
- **Co-assembly: a pre-trained SemiBin2 model on the pooled BAM, and per-sample depths for MetaBAT2.** After Binette this gave 36 bins with 25 MIMAG high, the same as the best self-trained run, with SemiBin2 at 0.65 CPU-hours instead of 24.8 (38× less; self-training cost more than the assembly). Per-sample coverage doubles MetaBAT2's high-quality bins (11 → 23); a pre-trained SemiBin2 model cannot use it (one coverage column only), and self-trained SemiBin2 gained nothing from it.
- **pplacer:** with `--pplacer_cpus 4`, one genome in a class-level tree took 60 GB, a few 124–131 GB, and fifty 141–142 GB. GTDB_TK_BATCH's 256 GB request covers it; 160 GB would cover every case measured. In the normal runs, 0 of 120 asm genomes and 2 of 219 co-assembly genomes reached pplacer.
- **Fix in 7b, whatever is decided:** COMEBin writes `FAILED.txt` and exits 0 on its first failure, so an out-of-memory kill is never retried and the user gets an empty binner silently.

## Summary

| Question | Measured | Recommendation |
|---|---|---|
| SemiBin2 2.1.0 → 2.6.0 | MIMAG high raw bins 18 → 23 (asm), 17 → 24 (coasm); same cost per sample; `train_self` 18.3 → 24.8 CPU-h; co-assembly mode breaks (`model.pt`) | pin, with `--random-seed` and the `model.pt` fix |
| MetaBAT2 2.15 → 2.18 | high 11 → 14 (asm), 11 → 11 (coasm); bin headers gain depth fields, which made MetaWRAP's bins wrong | pin, with headers cut |
| MMseqs2 15 → 18 | identical clusters at every identity level | pin |
| GTDB-Tk 2.7.2 → 2.8.0 | identical classifications (120 asm and 219 coasm genomes); skani 0.4.0 reads the R232 sketches | pin; biocontainer image (owner's choice) |
| Refiner (L1) | asm: Binette 45 bins, 26 high, 158.4 Mbp, 3.0 CPU-h, 11.2 GB vs MetaWRAP 44, 23, 152.9 Mbp, 23.0 CPU-h, 34.5 GB; coasm the same pattern | Binette |
| CheckM2 on raw bins (M7) | 4.5 CPU-h, 12 GB on the asm set; Binette's reports agree with it on every asm bin, and on all but 1–2 bins near the 90 % line in the co-assembly | take it from Binette |
| pplacer at 4 CPUs (M10) | 0 of 120 asm and 2 of 219 coasm genomes reach it; a class-level tree then takes 60 GB (1 genome) to 142 GB (53, forced) | keep the 256 GB request |
| Co-assembly coverage (M1) | per-sample depths: MetaBAT2 high 11 → 23; refined set unchanged (36 bins, 23–25 high) | per-sample depths for MetaBAT2 |
| SemiBin2 co-assembly mode | `train_self` 24.8 CPU-h, slower with more threads; pre-trained 0.65 CPU-h; with per-sample MetaBAT2, Binette gives 36 bins, 25 high either way | pre-trained model by default |
| COMEBin | 38 CPU-h on a 36.6 Mbp assembly (both versions); after Binette +1 genome, +1 high, +39 % Mbp on that sample; 1.1.0 broken as built; failures silent | keep as an option with a cost warning; fix failure handling |

The rest of this document gives the measurements behind each row.

## Setup

**Two sets, both from the test fixture (real reads, subsampled; clean reads after the pipeline's QC and host removal).** Assemblies are MEGAHIT's with its own settings, which are the Phase 5c defaults (no normalization, MEGAHIT's k-list, `--min-count 2`), filtered at 1 kb, with the pipeline's bwa-mem3 BAMs:

| Set | Contigs ≥ 1 kb | Mbp |
|---|---|---|
| asm, sample A | 25,366 | 117.8 |
| asm, sample B | 16,094 | 93.8 |
| asm, sample C (the sample measured for COMEBin in Phase 5) | 10,542 | 36.6 |
| asm, sample D | 2,078 | 3.7 |
| coasm, five samples pooled (the full-depth co-assembly of [assembly_2026.md](assembly_2026.md)) | 34,888 | 182.4 |

**Harness.** A scratch Nextflow workflow ran the pipeline's CALCULATE_DEPTH, METABAT2, SEMIBIN, COMEBIN and METAWRAP modules from dev `43d933b` on these fixed contigs and BAMs, then refined each binner combination with MetaWRAP, Binette and DAS Tool. Binette and DAS Tool were scratch modules written for this benchmark, with MetaWRAP's thresholds (completeness ≥ 50 %, contamination ≤ 10 %) for Binette and DAS Tool's own score threshold (0.5). Tool versions came from container overrides in config:

| | Current pins | Updated |
|---|---|---|
| SemiBin2 | `semibin:2.1.0` | `semibin:2.6.0--pyhdfd78af_0` |
| MetaBAT2 (and `jgi_summarize_bam_contig_depths`) | `metabat2:2.15--h4da6f23_2` | `metabat2:2.18--h38e344b_2` |
| COMEBin | `comebin:1.0.4--hdfd78af_0` | `comebin:1.1.0--hdfd78af_2` |
| GTDB-Tk | `gtdbtk:2.7.2--pyhdfd78af_1` | 2.8.0, two images (section (a)) |
| MMseqs2 | 15.6f452 | `mmseqs2:18.8cc5c--hd6d6fdc_0` |
| Binette | | `binette:1.2.1--pyh7e72e81_0` (CheckM2 1.1.0 inside) |
| DAS Tool | | `das_tool:1.1.7--r44hdfd78af_1` |

CheckM2 1.1.0 then ran once over every bin set of a batch (the pipeline's one-run pattern), and GTDB-Tk once over the distinct genomes of all refined sets (a bin with the same contigs as another is one genome).

**SemiBin2 was seeded.** It is not seeded in the pipeline: two unseeded 2.6.0 runs of the same sample differ by about one bin (15 vs 16, 38 vs 39). With `--random-seed 11`, reruns gave identical bins on all four samples. Every comparison below uses seed 11 for both versions, so old against new is not seed noise. The module would need `--random-seed` passed for this (7b).

**SemiBin2's pre-trained model.** The fixture's parameter files set `--semibin_env_model global`; the pipeline's default is `human_gut`. Every per-sample run and the co-assembly pre-trained-model run used `global`. `human_gut` was not tested; the co-assembly set is human gut, so it may do better there than `global` did.

**Workarounds in the harness only** (each is a 7b fix if the image is pinned; section (a)):

- METABAT2 cut each bin header to its first word (a no-op on 2.15).
- SEMIBIN in co-assembly mode used `model.pt` when there was no `model.h5`.
- COMEBIN 1.1.0 had a Python `realpath` on its PATH and `-d cpu`.

**Bin classes.** CheckM2 only. MIMAG high: completeness ≥ 90 % and contamination < 5 %. MIMAG medium: ≥ 50 % and < 10 %, not high. tRNA and rRNA were not checked, so "high" here is completeness and contamination only. "Passing" means high or medium.

## Machine and runs

- 4 × Xeon E7-4820 v4 (40 cores, 80 threads), 2 TB RAM; Docker; Nextflow 26.04.4. Each process got its pipeline label's first attempt: SemiBin2 and MetaBAT2 8 CPUs (`process_medium`); COMEBin and MetaWRAP 16 (`process_high`), and Binette and DAS Tool were given the same 16. CheckM2 and GTDB-Tk ran in their own containers at 16 CPUs.
- Other users kept the 1-minute load at 60–180 on 80 threads, and several runs overlapped. Wall times are inflated, often by 2–3×.
- **Time and memory come from the container's cgroup**, sampled every 2–3 s (GB here are 2^30 bytes): CPU-seconds and peak anonymous memory ("anon", the process memory). CPU-hours compare across runs better than wall time, but contention on shared cores inflates them too; differences under ~10 % are not meaningful here.
- **Nextflow's peak RSS double-counts.** For MetaWRAP it reports 68.7 GB, but the cgroup's anon peak was 34.5 GB; it sums the RSS of processes that share memory. All memory figures below are cgroup anon unless marked.

## (a) Tool updates

### What is inside the images

Each command the modules call was run with `--version` or `--help` inside the image.

| Image | Found inside |
|---|---|
| SemiBin2 2.6.0 | Python 3.14.8, PyTorch 2.14.1 (2.1.0: Python 3.12.2, PyTorch 2.1.2). Every subcommand and flag the module uses is still there: `single_easy_bin`, `generate_sequence_features_single`, `train_self`, `bin_short`. Both versions accept `--random-seed`. |
| MetaBAT2 2.18 | `metabat2` 2:2.18, `jgi_summarize_bam_contig_depths` 2.18; the same options as 2.15. |
| COMEBin 1.1.0 | Python 3.13.15, PyTorch 2.13.0 (a CPU-only build), numpy 2.5.2 (1.0.4: Python 3.7.12, PyTorch 1.10.2). New options: `-d` device (default `cuda`), `-s` seed. `run_comebin.sh --version` is not an option in either version, so the module's versions line records an error message. |
| MMseqs2 18 | `mmseqs` 18.8cc5c; `easy-cluster` keeps every option the module uses. |
| GTDB-Tk 2.8.0, biocontainer `gtdbtk:2.8.0--pyhdfd78af_0` | The only 2.8.0 biocontainer. **skani 0.4.0.** A depot.galaxyproject.org sif of the same build exists. |
| GTDB-Tk 2.8.0, Seqera Containers build `community.wave.seqera.io/library/gtdbtk_skani:d13cc293988742dd` | Built on request from bioconda `gtdbtk=2.8.0` and `skani=0.3.2`: skani 0.3.2, pplacer v1.1.alpha19, Prodigal 2.6.3, HMMER 3.4, FastTree 2.2.0, Python 3.13.17, and the shell tools GTDB_TK_BATCH calls. Its Singularity twin has no `find`, which GTDB_TK_BATCH calls, so the docker URI would be pinned for every engine (design doc Q20). |
| GTDB-Tk 2.7.2 (current) | skani 0.3.2, the same pplacer. |

**skani 0.4 reads the R232 reference sketches.** The plan required skani 0.3.x because 0.4 was expected not to read them. It did: `skani search` against the R232 sketch database found the same hit for a test bin with both versions (ANI 99.28), and the biocontainer gave the same classifications as the other two images on every genome below. So both 2.8.0 images work and gave identical results here; they differ in skani version (0.3.2 vs 0.4.0) as well as packaging. Only the Seqera Singularity build was checked for the commands GTDB_TK_BATCH calls; the depot sif of the biocontainer was not, and 7b checks it and runs it once under apptainer before pinning it.

### Raw binners

**asm set, four samples summed:**

| Binner | Bins | Passing | MIMAG high / medium | Mbp in passing bins | CPU-s (cgroup) | Peak anon |
|---|---|---|---|---|---|---|
| SemiBin2 2.1.0 | 103 | 45 | 18 / 27 | 153.4 | 3,328 | 1.53 GB |
| SemiBin2 2.6.0 | 93 | 45 | 23 / 22 | 159.0 | 3,087 | 0.88 GB |
| MetaBAT2 2.15 | 87 | 37 | 11 / 26 | 119.2 | 100 (+163 depth) | 0.50 GB |
| MetaBAT2 2.18 | 83 | 38 | 14 / 24 | 124.7 | 117 (+110 depth) | 0.50 GB |

**coasm set, pooled BAM:**

| Binner | Bins | Passing | MIMAG high / medium | Mbp in passing bins |
|---|---|---|---|---|
| SemiBin2 2.1.0 `train_self` | 76 | 35 | 17 / 18 | 111.7 |
| SemiBin2 2.6.0 `train_self` | 66 | 35 | 24 / 11 | 117.6 |
| MetaBAT2 2.15 | 66 | 32 | 11 / 21 | 92.3 |
| MetaBAT2 2.18 | 65 | 34 | 11 / 23 | 96.2 |

- **SemiBin2 2.6.0 makes better bins.** The same number pass, but more of them are high quality (asm 23 vs 18, coasm 24 vs 17), and it makes fewer bins that fail. Per sample it costs the same.
- **MetaBAT2 2.18 is slightly better or the same** (asm 14 vs 11 high). With the module's `--seed 3`, it also gives different bins from 2.15 on the same depth file (samples A and B: 41 vs 43 and 34 vs 36 bins).
- **SemiBin2 2.6.0's `train_self` costs more than 2.1.0's:** 24.8 vs 18.3 CPU-hours on the co-assembly (section (g)).

### What breaks in the modules with the updated images

These were worked around in the harness only. Each needs a fix in 7b if the image is pinned.

- **SemiBin2 2.6.0 in co-assembly mode.** `train_self` saves `model.pt` (2.1.0: `model.h5`). The module then calls `bin_short --model contig_output/model.h5`, which fails.
- **MetaBAT2 2.18 bin headers.** 2.18 appends tab-separated depth fields to each contig name in the bin FASTA headers (`>k141_16839<TAB>total_depth=11.75<TAB>…`; 2.15 writes the bare name). DAS Tool stopped on them. **MetaWRAP accepted them and gave wrong bins:** it matches contigs by the whole header, so a MetaBAT2 2.18 contig never matched the same contig in a SemiBin2 bin, and both copies were kept. With the raw 2.18 headers it returned 31 / 9 / 41 bins for samples B / C / A, against 17 / 5 / 22 with 2.15, and sample A's bins held 7,788 contig headers for 4,447 distinct contigs. Binette had no duplicates. Cutting the header to its first word in METABAT2 (`sed -i '/^>/s/[[:space:]].*//'`, a no-op on 2.15) fixed it; every 2.18 result in this document is with that fix. Anything else that reads bin headers whole needs the same.
- **COMEBin 1.1.0 does not run as built, in two ways.** (1) The image's `realpath` is busybox's, which has no `-e`, `-m` or `--`, and `run_comebin.sh` calls `realpath -e -- <contigs>`: it stops after 17 s. All three 1.1.0 builds share that base. (2) With a working `realpath`, the default device `cuda` fails at training on the CPU-only PyTorch ("CUDA training was requested, but the installed PyTorch is CPU-only"). The harness put a Python `realpath` first on the PATH and passed `-d cpu`.
- **COMEBin failures are silent (both versions).** The module runs `run_comebin.sh` under `set +e`, and on a non-zero exit writes `FAILED.txt` into an empty bin directory and exits 0. Nextflow sees a success, so attempts 2 and 3 (smaller batches, meant for memory failures) never run, and an out-of-memory kill (exit 137) is not retried either. Both 1.1.0 failures above ended this way on attempt 1, with no error in the run.

### MMseqs2 18

The pipeline's CLUSTERING module was run with both versions on the ARG proteins (Pyrodigal → DeepARG → the pipeline's formatter) of 68 Binette bins from both sets: 456 proteins. The clusters at 90 / 95 / 99 / 100 % identity were 230 / 236 / 246 / 254 with both, and the sorted cluster tables were identical at every level. Both took 41 s and 0.88 GB.

### GTDB-Tk 2.8.0

GTDB-Tk ran on the distinct genomes of every refined set, at 16 CPUs and `--pplacer_cpus 4`, with the R232 reference data:

| Set | Image | Wall | CPU-s | Peak anon | Classifications |
|---|---|---|---|---|---|
| asm, 120 genomes | 2.7.2 | 12.5 min | 774 | 39.8 GB | 120 to species |
| | 2.8.0 Seqera (skani 0.3.2) | 11.1 min | 691 | 39.9 GB | identical |
| | 2.8.0 biocontainer (skani 0.4.0) | 10.4 min | 646 | 39.3 GB | identical |
| coasm, 219 genomes | 2.7.2 | 48.8 min | 3,074 | 59.8 GB | 217 to species, 2 to genus |
| | 2.8.0 Seqera | 46.1 min | 2,908 | 59.7 GB | identical |
| | 2.8.0 biocontainer | 45.5 min | 2,839 | 59.6 GB | identical |

Classifications and ANI values were identical in all three images on both sets, so no classification changed. The coasm runs had more genomes, and took more memory because two of them went through pplacer (section (c)).

## (b) Refinement: MetaWRAP, Binette, DAS Tool

All three refined the same SemiBin2 2.6.0 + MetaBAT2 2.18 bins. "MetaWRAP, current tools" is the pipeline today: MetaWRAP on SemiBin2 2.1.0 + MetaBAT2 2.15.

**asm set, four samples summed:**

| Refiner | Bins | Passing | MIMAG high / medium | Mbp in passing bins | Median completeness / contamination | Wall (sum) | CPU-h | Peak anon |
|---|---|---|---|---|---|---|---|---|
| MetaWRAP, current tools | 44 | 44 | 19 / 25 | 149.4 | 90.2 / 0.94 | 546 min | 22.8 | 34.5 GB |
| MetaWRAP | 44 | 44 | 23 / 21 | 152.9 | 94.0 / 0.97 | 492 min | 23.0 | 34.5 GB |
| **Binette** | 45 | 45 | **26 / 19** | **158.4** | 96.2 / 1.23 | 40 min | 3.0 | 11.2 GB |
| DAS Tool | 43 | 41 | 23 / 18 | 144.0 | 96.1 / 0.91 | 17 min | 1.4 | 0.57 GB |

Per sample (bins from MetaWRAP / Binette / DAS Tool): A 22 / 22 / 22; B 17 / 17 / 16; C 5 / 6 / 4; D 0 / 0 / 1. Sample D's one DAS Tool bin is 36.8 % complete: DAS Tool scores bins with its own single-copy genes and applies no completeness filter.

**coasm set, pooled BAM, SemiBin2 `train_self`:**

| Refiner | Bins | Passing | MIMAG high / medium | Mbp in passing bins | Wall | CPU-h | Peak anon |
|---|---|---|---|---|---|---|---|
| MetaWRAP, current tools | 35 | 35 | 16 / 19 | 112.3 | 3 h 51 min | 14.7 | 34.6 GB |
| MetaWRAP | 35 | 35 | 21 / 14 | 116.5 | 3 h 50 min | 15.7 | 34.6 GB |
| **Binette** | 36 | 36 | **25 / 11** | **118.7** | 23 min | 1.8 | 11.9 GB |
| DAS Tool | 34 | 33 | 24 / 9 | 107.9 | 14 min | 1.05 | 0.58 GB |

The same pattern held for the other two co-assembly inputs (section (d)): Binette gave 36 bins each time, all passing, and MetaWRAP's cost was the same.

**Do the refiners find the same genomes?**

- Pairing bins by shared contigs (Jaccard ≥ 0.9) on the asm set: 37 of Binette's 45 bins have a MetaWRAP twin, and 29 have a DAS Tool twin. 31 of MetaWRAP's 44 have a DAS Tool twin.
- By GTDB species (per sample): MetaWRAP gave the same 44 species with the current and the updated binners. Binette gave those 44 plus one more (sample C). DAS Tool lost three and found two, one of them the 36.8 %-complete bin. All three together: 47.
- In the co-assembly, every MetaWRAP set (current tools, and the three inputs of section (d)) had the same 35 species. Binette had those 35 plus one, with each of the three inputs. DAS Tool lost one to three and found one.

**Binette chooses bins by CheckM2, and is judged here by CheckM2.** So its lead in CheckM2 classes is partly by construction. The GTDB-Tk result does not depend on CheckM2: it found every species MetaWRAP found, and one more. MetaWRAP chooses by CheckM 1 and DAS Tool by its own single-copy gene score.

**Cost.** Binette took an eighth to a ninth of MetaWRAP's CPU time and a third of its memory on both sets; DAS Tool a fifteenth and a sixtieth. MetaWRAP's cost grows with the number of binners (section (f): 2.2× with a third binner, where Binette's grew 1.4×).

## (c) GTDB-Tk: genomes reaching pplacer, and its memory

**In the normal runs, few genomes reach pplacer.** The ANI screen (skani against the species representatives) classified all 120 asm genomes, so none reached pplacer, and the run peaked at ~40 GB. It classified 217 of the 219 coasm genomes. The other two, one bacterial and one archaeal, went through pplacer, the bacterial one into a class-level tree. That raised the coasm runs' peak to ~60 GB (backbone 4 min, class-level tree 20 min, for the one genome).

**The worst case, forced.** To measure pplacer with more genomes, the 120 asm genomes were run again with 2.8.0 (Seqera image) and `--place_species`, which sends every genome through pplacer even when the ANI screen has classified it:

| Stage | Genomes | Wall | Peak anon |
|---|---|---|---|
| ANI screen, gene calling, alignment | 120 | 32 min | 39.8 GB |
| Backbone tree | 120 | 17 min | 37.3 GB |
| Class-level tree 2 | 53 | 69 min | 142.4 GB |
| Class-level tree 4 | 52 | 62 min | 141.2 GB |
| Class-level tree 6 | 8 | 45 min | 128.5 GB |
| Class-level tree 1 | 4 | 38 min | 130.8 GB |
| Class-level tree 5 | 2 | 49 min | 123.6 GB |
| Class-level tree 3 | 1 | 18 min | 58.7 GB |
| Whole run | 120 | 5 h 29 min (9.0 CPU-h) | 142.4 GB |
| *For comparison, coasm normal run:* class-level tree 4 | 1 | 20 min | 59.7 GB |

- **Forced placement changed no classification** (120 of 120 the same as the normal run).
- **Memory rises steeply with the first few genomes in a class-level tree, then levels off.** One genome took 59–60 GB (trees 3 and 4); two to eight genomes in other trees 124–131 GB; 52–53 genomes 141–142 GB. Phase 5 measured 142.5 GB at `--pplacer_cpus 4` for its class-level step.
- **A class-level tree took 38–69 min for two to fifty-three genomes, 18–20 min for one,** in the forced run, which placed genomes into six trees one after another (4 h 41 min in the six trees, 4 h 58 min with the backbone).
- **For 7b:** GTDB_TK_BATCH's request (`process_high_memory`, 256 GB) covers this, with room. A run whose genomes all pass the ANI screen needs ~40 GB and one with a few genomes in a class-level tree ~60–130 GB, but neither is known in advance; 160 GB would cover every case measured here.

## (d) Co-assembly coverage: pooled BAM or one BAM per sample

Today the binners of a co-assembly get one BAM: all samples' reads pooled and aligned once. The alternative aligns each sample to the co-assembly on its own (five bwa-mem3 runs, 22.5 min of task time in total and 3.09 GB of BAMs, against 28.6 min and 3.0 GB for the pooled alignment in Phase 5b). The binners then see one coverage column per sample, which tells apart genomes that co-occur in some samples and not others. Phase 7c makes these per-sample alignments anyway (issue #28), so the coverage itself would cost nothing extra.

A pre-trained SemiBin2 model takes only one coverage column ("can only be used for single-sample binning"), so with per-sample BAMs SemiBin2 has to self-train. The table has the three runs made, and a fourth that combines MetaBAT2's per-sample bins with the pre-trained SemiBin2 bins:

| Binner input | MetaBAT2 bins (high / medium) | SemiBin2 bins (high / medium) | SemiBin2 CPU-h | Binette (high / medium; Mbp) | DAS Tool (high / medium) | MetaWRAP (high / medium) |
|---|---|---|---|---|---|---|
| pooled BAM, SemiBin2 `train_self` | 65 (11 / 23) | 66 (24 / 11) | 24.8 | 36 (25 / 11; 118.7) | 34 (24 / 9) | 35 (21 / 14) |
| pooled BAM, SemiBin2 pre-trained | 65 (11 / 23) | 64 (22 / 11) | 0.65 | 36 (23 / 13; 120.5) | 33 (22 / 9) | 35 (21 / 13) |
| per-sample BAMs, SemiBin2 `train_self` | **55 (23 / 14)** | 67 (23 / 12) | 24.8 | 36 (24 / 12; 124.3) | 35 (22 / 11) | 35 (22 / 12) |
| **MetaBAT2 per-sample, SemiBin2 pre-trained on the pooled BAM** | 55 (23 / 14) | 64 (22 / 11) | 0.65 | **36 (25 / 11; 121.7)** | | |

- **Per-sample coverage doubles MetaBAT2's high-quality bins** (11 → 23; 96 → 120 Mbp in passing bins), with fewer bins in all (55 vs 65).
- **It barely moves SemiBin2** (23 / 12 vs 24 / 11), **or the refined set:** Binette gave 36 bins with 23–25 high in every case, and 1–6 Mbp more with per-sample coverage.
- **What it costs is the pre-trained model.** Per-sample coverage with SemiBin2 means `train_self`, 24.8 CPU-hours, against 0.65 for the pre-trained model on the pooled BAM.
- **The combination keeps both.** MetaBAT2 on the per-sample depths and the pre-trained SemiBin2 on the pooled BAM, refined by Binette, gave 36 bins, 25 high and 11 medium, 121.7 Mbp: as many high-quality bins as the best self-trained run, with SemiBin2 at 0.65 CPU-hours instead of 24.8. Binette took 15 min, 1.9 CPU-hours and 11.9 GB. It was run outside Nextflow, with the harness's Binette command on the two runs' published bins.
- Binette's own CheckM2 report and the separate CheckM2 run disagree on one or two bins near the 90 % line for these sets (for example 23 / 13 vs 25 / 11 for the pooled `train_self` run), so the four Binette results (23–25 high, 36 bins each) are within noise of each other. What separates them is cost.

## (e) CheckM2 on the raw binner sets

The pipeline runs CheckM2 once over the raw binner bins and the refined bins together (`ch_binner_bins.mix(ch_refined_bins)`).

**Cost on the asm set** (two CheckM2 runs at the same time, 16 CPUs each, load 110–130):

| Bins checked | Wall | CPU-h | Peak anon |
|---|---|---|---|
| raw SemiBin2 + MetaBAT2 bins, 176 | 33.0 min | 4.5 | 12.3 GB |
| Binette's refined bins, 45 | 19.1 min | 2.5 | 11.0 GB |

Most of the cost is fixed (loading the DIAMOND database), so CheckM2 on 4× more bins cost less than 2× more.

**Binette already reports it.** Binette runs CheckM2 1.1.0 (the pipeline's version, with the same database) on every input bin and every bin it builds, and writes `input_bins_quality_reports/*.tsv` and `final_bins_quality_reports.tsv`. Against the separate CheckM2 run:

- raw bins: 176 of 176 in the same class (high / medium / below), median completeness difference 0.07 points (largest 2.6), contamination difference at most 1.33;
- refined bins: 45 of 45 in the same class, largest difference 1.58 points.

So with Binette, the pipeline can take both reports from Binette and skip the separate CheckM2 run whenever refinement runs. It is still needed when only one binner is chosen (no refinement). With MetaWRAP (which uses CheckM 1) or DAS Tool (its own score), the separate run stays.

## (f) COMEBin

BugBuster runs COMEBin on CPU only: the module has no GPU wiring, and this server has no GPU. It was run on two asm samples, C (the one measured in Phase 5) and D, at `process_high`'s 16 CPUs. Version 1.1.0 ran with the two harness fixes (`realpath`, `-d cpu`).

**Cost per sample:**

| Version, sample | Wall | CPU-h | Cores used | Peak anon | Bins (high / medium) |
|---|---|---|---|---|---|
| 1.0.4, sample C | 7 h 49 min | 38.6 | 4.9 | 7.15 GB | 17 (2 / 5) |
| 1.1.0, sample C | 6 h 00 min | 38.2 | 6.4 | 2.73 GB | 12 (2 / 5) |
| 1.0.4, sample C, Phase 5 (pipeline trace) | 5 h 54 min | ~41 | | 7 GB | 16 |
| 1.0.4, sample D | 1 h 52 min | 5.9 | 3.2 | 4.0 GB | 7 (0 / 1) |
| 1.1.0, sample D | 1 h 40 min | 5.7 | 3.4 | 2.67 GB | 6 (0 / 1) |
| SemiBin2 2.6.0, sample C (8 CPUs) | 4.8 min | 0.10 | 1.3 | 0.55 GB | 16 (1 / 5) |
| MetaBAT2 2.18, sample C | 0.1 min | 0.01 | | 0.15 GB | 6 (1 / 2) |

- **1.1.0 costs the same CPU time as 1.0.4** (38.2 vs 38.6 CPU-hours on sample C). Its shorter wall time came from more cores under a different load, not from less work. It needed less memory: 2.7 vs 7.2 GB on sample C, 2.7 vs 4.0 GB on sample D.
- **On sample C, COMEBin costs ~370× SemiBin2's CPU time** for one more high-quality bin before refinement (2 vs 1) and the same number of medium ones.
- It used at most 6.4 of its 16 CPUs on average in these runs, and 6.9 in Phase 5's run of sample C. If it is kept, 8 CPUs is the measured request.
- Its CPU time grows with the assembly: 5.7–5.9 CPU-hours for sample D's 3.7 Mbp of contigs, 38 for sample C's 36.6 Mbp. At that rate sample A (118 Mbp) would take on the order of 100 CPU-hours; this was not measured.

**What it adds after refinement** (sample C; SemiBin2 + MetaBAT2 against SemiBin2 + MetaBAT2 + COMEBin, all 2.6.0 / 2.18 / 1.1.0):

| Refiner | Without COMEBin: bins (high / medium; Mbp) | With COMEBin | Refiner CPU-h without / with |
|---|---|---|---|
| **Binette** | 6 (1 / 5; 19.2) | **7 (2 / 5; 26.6)** | 0.41 / 0.56 |
| DAS Tool | 4 (1 / 3; 13.9) | 4 (1 / 3; 14.7) | 0.17 / 0.17 |
| MetaWRAP | 5 (1 / 4; 15.8) | not valid (below) | 3.0 / 6.5 |

- **With Binette, COMEBin added one genome and improved two bins on sample C.** The new bin (1,646 contigs, 84 % complete, 9 % contamination) shares no more than 2 % of its contigs with any bin made without COMEBin. One bin went from 61 % to 91 % complete (medium to high), and another from 56 % to 80 %. Sequence in passing bins rose by 39 % (19.2 → 26.6 Mbp).
- **On sample D, every refiner had a medium-quality bin with COMEBin (80 % complete, 1.0 Mbp) and none without it.**
- **MetaWRAP's result with COMEBin 1.1.0 is not valid.** COMEBin 1.1.0 writes each contig's whole assembly header into its bins (`>k141_3986 flag=1 multi=3.0000 len=1059`); 1.0.4 and the other binners write the bare name. MetaWRAP matches whole headers, so it kept COMEBin's copy of a contig next to SemiBin2's or MetaBAT2's: its 12 bins held 42.2 Mbp, more than the 36.6 Mbp assembly. It is the MetaBAT2 2.18 problem again (section (a)), and needs the same fix in COMEBIN if 1.1.0 is pinned. Binette reads only the first word of each header, and the harness built DAS Tool's contig table the same way, so neither was affected.

**No FAILED.txt in the runs that worked.** Every `FAILED.txt` seen came from the two 1.1.0 breakages above, each on attempt 1, each with exit 0.

## (g) SemiBin2 in co-assembly mode

In co-assembly mode the module runs `generate_sequence_features_single`, then `train_self` (a model trained on this assembly) and `bin_short`. Measured here, pooled BAM, 8 CPUs:

| Run | Wall | CPU-h | Cores used | Peak anon | Bins (high / medium) |
|---|---|---|---|---|---|
| 2.1.0 `train_self` | 7 h 11 min | 18.3 | 2.55 | 20.5 GB | 76 (17 / 18) |
| 2.6.0 `train_self` | 8 h 32 min | 24.8 | 2.91 | 19.6 GB | 66 (24 / 11) |
| 2.6.0 `train_self`, per-sample BAMs | 8 h 26 min | 24.8 | 2.94 | 21.7 GB | 67 (23 / 12) |
| 2.6.0 pre-trained model (`single_easy_bin --environment global`) | 28 min | 0.65 | 1.40 | 1.13 GB | 64 (22 / 11) |

MEGAHIT took 16.9 CPU-hours for this co-assembly ([assembly_2026.md](assembly_2026.md)), so self-training costs more than the assembly.

**Where the time goes.** Phase 5b's 2.1.0 log: features 12.5 min, `train_self` 9 h 43 min of 10 h 14 min, binning and reclustering ~13 min. With 2.6.0, features took 16 min (0.32 CPU-h, 0.41 GB). The training draws up to 4 M cannot-link contig pairs per epoch (11,007 contigs in the training data), for 15 epochs in batches of 2,048: it is PyTorch training of a small network on many pairs. A profile of a reduced run (500 k pairs) at 1 thread: backward pass 81 s, dropout 44 s, linear layers 38 s, data loading ~25 s of 159 s. At 8 threads the same work took 949 CPU-s, mostly in dropout and batch norm.

**More CPUs do not help.** One-epoch runs of `train_self`, at the same time, load 90–160:

| Threads | Wall | CPU-s | Average CPU | Peak anon |
|---|---|---|---|---|
| 1 | 28.0 min | 1,109 | 66 % | 11.3 GB |
| 4 | 23.1 min | 2,822 | 204 % | 9.3 GB |
| 8 | 26.9 min | 5,512 | 342 % | 10.1 GB |
| 16 | 38.0 min | 13,887 | 609 % | 11.2 GB |

Four threads were the fastest, by 18 % over one, for 2.5× the CPU time. Eight and sixteen were slower than four and cost 2–5× more. The network is too small for PyTorch's threads to pay off, and on a loaded machine they make it worse.

**A pre-trained model instead.** `--environment` skips training: 28 min and 0.65 CPU-hours instead of 8.5 h and 24.8, and 1.1 GB instead of 20. It made 64 bins (22 high) against 66 (24 high). After Binette with MetaBAT2's bins, 36 bins in both (23 / 13 vs 25 / 11 high / medium; 120.5 vs 118.7 Mbp). The model used was `global`; the pipeline's default `human_gut` matches this set and was not tested.

## Recommendations

Approved by the owner on 2026-10-26 (see Decision, top); Phase 7b applies them.

**L1, refiner: Binette 1.2.1, replacing MetaWRAP.** On both sets it gave as many or more bins, MIMAG high bins and sequence in passing bins as MetaWRAP, and every GTDB species MetaWRAP found, for an eighth of the CPU time and a third of the memory. It also reads only the first word of each contig header, which protects it from the MetaBAT2 2.18 and COMEBin 1.1.0 header changes. Keep its thresholds at the current `--metawrap_completeness 50` and `--metawrap_contamination 10` (renamed or aliased in 7b). The new module must keep the refined-bin directory contract (CHECKM2_BATCH, GTDB_TK_BATCH, the BIN_* reports, Bakta, ARG bin clustering, `SKIPPED.txt`); GTDB_TK_BATCH finds refined bins by MetaWRAP's directory suffix today. DAS Tool is the cheaper alternative, but it lost bins and kept one 36.8 %-complete bin; if it were chosen, 7b would have to add a completeness and contamination filter.

**M7, CheckM2 on raw bins: take it from Binette.** With Binette, the separate CHECKM2_BATCH run is redundant whenever refinement runs (176 / 176 and 45 / 45 bins in the same class). The pipeline can publish Binette's input-bin and final-bin reports in CHECKM2_BATCH's place, saving 4.5 CPU-hours and a 12 GB task on the asm set. CHECKM2_BATCH stays for single-binner runs.

**U5, pins:**

- SemiBin2 2.6.0. 7b: use `model.pt` in co-assembly mode; pass `--random-seed` (a fixed default) so reruns give the same bins.
- MetaBAT2 2.18. 7b: cut bin headers to the first word in METABAT2.
- MMseqs2 18.8cc5c.
- GTDB-Tk 2.8.0, from one of two images that gave identical results here: the Seqera build with skani 0.3.2 (docker URI for every engine, as its sif lacks `find`) or the biocontainer with skani 0.4.0 (a depot sif exists but was not checked here). The owner chose the biocontainer; 7b checks its depot sif for every command GTDB_TK_BATCH calls and runs it once under apptainer before pinning it.
- COMEBin 1.1.0 only if COMEBin stays, and only with all of: a `realpath` that supports `-e` (or a patched `run_comebin.sh`), `-d cpu`, and bin headers cut to the first word. Its PyTorch is the CPU build, so `-d cpu` loses nothing. Otherwise keep 1.0.4.

**COMEBin: keep as an option, not a default (it is not one today), with a cost warning.** After Binette it added a genome and upgraded a bin on sample C, and gave the only passing bin on sample D. It costs 1–1.6 CPU-hours per Mbp of contigs (38 CPU-hours for sample C, ~370× SemiBin2) and used at most 6.4–6.9 of 16 CPUs on average (here and in Phase 5): request 8 CPUs, and say in README and `docs/manual.md` what it costs. Removing it is the alternative if that cost is judged too high for one or two bins per sample. Either way, 7b fixes its failure handling: fail the task (so Nextflow retries with the smaller batch settings) and write `FAILED.txt` only after the last attempt.

**M10, pplacer threads 4 (decided after Phase 5).** GTDB_TK_BATCH passes `--pplacer_cpus 1` today; 7b adds a parameter for it with default 4. Every GTDB-Tk run here used 4. Keep GTDB_TK_BATCH's 256 GB request (160 GB would cover every case measured, with little margin): one genome in a class-level tree took 60 GB, a few 124–131 GB, fifty 142 GB.

**SemiBin2 in co-assembly mode: use a pre-trained model by default** (`single_easy_bin --environment` on the pooled BAM). It saves ~24 CPU-hours per co-assembly, more than MEGAHIT costs, for no loss after Binette when MetaBAT2 gets per-sample depths. `train_self` stays available behind a parameter. If self-training is kept anywhere, give it 4 CPUs, not 8: more threads made it slower.

**M1, co-assembly coverage: per-sample depths for MetaBAT2, the pooled BAM for SemiBin2.** In 7c, CALCULATE_DEPTH takes the per-sample BAMs (made anyway for issue #28) and gives MetaBAT2 one depth column per sample; SemiBin2's pre-trained model keeps the pooled BAM (7c's merge of the per-sample BAMs). Measured: MetaBAT2 high-quality bins 11 → 23, refined 36 bins with 25 high at SemiBin2's pre-trained cost. With `train_self`, SemiBin2 could take the per-sample BAMs too, but that gained nothing here.

## Caveats

- Two small sets: four per-sample assemblies and one co-assembly. COMEBin ran on two samples. Each variant ran once; the seeded SemiBin2 runs are reproducible, but MetaBAT2 2.15 vs 2.18 and every refiner were single runs.
- Bin classes are CheckM2's. Binette selects bins with CheckM2, so it is judged here by its own criterion; GTDB-Tk species give a check that does not depend on CheckM2. Bins near the 90 % line move between high and medium between CheckM2 runs (Binette's own CheckM2 report and the separate run differ by up to 2.6 points per bin), so differences of one or two high-quality bins are within noise.
- SemiBin2's pre-trained model was `global`, not the pipeline's default `human_gut`.
- The machine was loaded by other users throughout (load 60–180 on 80 threads); wall times are inflated and CPU-hours carry some contention.
- The co-assembly is the full-depth co-assembly of a five-sample set, not a large one; SemiBin2's training cost and the per-sample coverage gain may scale differently with more samples.
