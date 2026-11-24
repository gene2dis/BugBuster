# YAML Parameters Guide: Writing a `-params-file`

A params file holds the pipeline's parameters in YAML (or JSON), so a run can be repeated
and shared without a long command line.

## Basic Usage

```yaml
# params.yaml: keys are the parameter names without the leading --
input: "samplesheet.csv"
output: "results"
assembly_mode: "none"        # read-level results only
nonpareil: true
rgi_prediction: true
max_memory: "64.GB"
```

```bash
nextflow run main.nf -params-file params.yaml -profile docker
```

A container profile (`-profile docker`, `singularity`, `apptainer` or `podman`, alone or
with a site profile such as `slurm`) is needed in every run; without one, every task runs
on the host and fails at the first tool that is not installed.

## Types and Values

- **On/off options** are YAML booleans, unquoted: `include_binning: true`,
  `nonpareil: false`.
- **Numbers** unquoted: `min_read_sample: 20000000`, `bbmap_length: 1500`.
- **Text and paths** in quotes, safest for anything with `:`, `#`, `,` or spaces in it:
  `binners: "semibin,metabat2"`, `max_memory: "128.GB"`, `max_time: "48.h"`.
- **Comma-separated lists** are one quoted string, not a YAML list:
  `binners: "semibin,metabat2"`, `custom_host_fasta: "/g/human.fa,/g/mouse.fa"`.
  The one parameter that takes a YAML list is `custom_sourmash_db` (the k-mer database
  and its lineages file: `custom_sourmash_db: ["/db/kmers.zip", "/db/lineages.csv"]`).
- **Paths** are absolute, or relative to the directory you launch Nextflow in.
- **Unknown keys stop the run** at launch, as unknown `--options` do: the file is
  validated against the pipeline's parameter schema
  ([`parameters.md`](parameters.md) lists every parameter).

The file keeps these types on every Nextflow version. On the command line, Nextflow 26.04
and later pass every value as text (`--nonpareil false` is the string `"false"`); the
pipeline reads its options so that both forms work, but a params file is the one place
where `false` and numbers are typed exactly as written.

## Combining With the Command Line and Config Files

Checked on Nextflow 25.10 and 26.04: a `--parameter` on the command line overrides the
same key in the params file, and the params file overrides `params` set in a config file
passed with `-c` (and the pipeline's defaults). Profiles are chosen on the command line
only (`-profile`), not in the params file.

```bash
# The file's settings, with another output folder and more memory for this run
nextflow run main.nf -params-file params.yaml --output results_b --max_memory 256.GB -profile docker

# Resume a failed run: the same command plus -resume
nextflow run main.nf -params-file params.yaml -profile docker -resume

# Check the parameters and the workflow without running any task
nextflow run main.nf -params-file params.yaml -profile docker -preview
```

Resource and executor settings that are not parameters (a process's memory, the SLURM
queue size) go in a config file passed with `-c`; see the manual,
[Advanced Configuration](manual.md#11-advanced-configuration).

## Checking a Params File

```bash
# YAML syntax
python3 -c "import yaml; print(yaml.safe_load(open('params.yaml')))"

# Parameter names and types, without running any task
nextflow run main.nf -params-file params.yaml -profile docker -preview

# Also the input files: a stub run reads the samplesheet and runs no tool
nextflow run main.nf -params-file params.yaml -profile docker -stub --output stub_check
```

## Examples

- Several host genomes for decontamination, and a pre-built combined index:
  [`DECONTAMINATION_QUICK_REFERENCE.md`](DECONTAMINATION_QUICK_REFERENCE.md#multiple-host-genomes-in-a-params-file),
  and the example file `examples/params_multiple_genomes.yaml`.
- Choosing binners (`binners: "semibin,metabat2"`): the manual,
  [Binning](manual.md#57-binning).
- A full run with binning as a params file: the manual, Example 3b
  ([Usage Examples](manual.md#10-usage-examples)).

Problems with a params file (a file not found, an option that seems ignored, a memory
limit): [`troubleshooting.md`](troubleshooting.md#inputsamplesheet-errors).
