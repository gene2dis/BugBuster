# Contributing to BugBuster

Thank you for your interest in contributing to BugBuster! This document provides guidelines for contributing to the project.

## Table of Contents

- [Code of Conduct](#code-of-conduct)
- [Getting Started](#getting-started)
- [Development Setup](#development-setup)
- [Making Changes](#making-changes)
- [Submitting Changes](#submitting-changes)
- [Coding Standards](#coding-standards)
  - [Nextflow/Groovy](#nextflowgroovy)
  - [Adding a parameter](#adding-a-parameter)
  - [Software versions](#software-versions)
  - [Databases](#databases)
  - [Process Structure](#process-structure)
  - [Module meta.yml](#module-metayml)
  - [Results summary](#results-summary)
  - [Workflow diagrams](#workflow-diagrams)
- [Testing](#testing)
  - [Local Testing](#local-testing)
  - [Adding a test](#adding-a-test)
  - [What CI runs](#what-ci-runs)
  - [Linting](#linting)

## Code of Conduct

Please be respectful and constructive in all interactions. We welcome contributions from everyone regardless of experience level.

## Getting Started

1. Fork the repository on GitHub
2. Clone your fork locally:
   ```bash
   git clone https://github.com/YOUR_USERNAME/BugBuster.git
   cd BugBuster
   ```
3. Add the upstream repository:
   ```bash
   git remote add upstream https://github.com/gene2dis/BugBuster.git
   ```
4. Start your work from `dev`, the development branch every pull request goes into
   (`main` holds releases):
   ```bash
   git fetch upstream && git switch -c feature/your-feature upstream/dev
   ```

## Development Setup

### Requirements

- Nextflow >= 25.10.0
- Docker or Singularity
- Git

### Running Tests

```bash
# Run with test profile
nextflow run main.nf -profile test,docker

# Dry run (stub mode)
nextflow run main.nf -profile test,docker -stub

# nf-test suite
nf-test test
```

`-profile test` reads `assets/test_samplesheet.csv`, which points to nf-core test reads
on GitHub (downloaded at launch). CI passes `--input` with a samplesheet of the committed
reads in `tests/data/reads/` instead; see [What CI runs](#what-ci-runs).

> **Note**: The pipeline is strict-syntax clean (the parser that is default
> since Nextflow 26.04). Keep it that way: run `nextflow lint .` on your changes
> and fix any errors it reports before opening a PR.

## Making Changes

### Branch Naming

Use descriptive branch names with a type prefix, as the repository's branches do:
- `feature/add-new-tool` - New features
- `fix/semibin-no-long-contigs` - Bug fixes
- `docs/post-roadmap-sync` - Documentation
- `refactor/module-cleanup` - Code refactoring
- `roadmap/p19-report-coverage` - Maintainers' roadmap phases

### Commit Messages

Follow conventional commit format:

```
type(scope): subject

body (optional)

footer (optional)
```

Types:
- `feat`: New feature
- `fix`: Bug fix
- `docs`: Documentation
- `refactor`: Code refactoring
- `test`: Adding tests
- `chore`: Maintenance

Example:
```
feat(modules): add DIAMOND module for protein alignment

- Added DIAMOND blastp process
- Added meta.yml descriptor
- Added stub block for testing
```

## Submitting Changes

1. Create a new branch from `dev`:
   ```bash
   git switch -c feature/your-feature upstream/dev
   ```

2. Make your changes and commit:
   ```bash
   git add .
   git commit -m "feat: your feature description"
   ```

3. Push to your fork:
   ```bash
   git push origin feature/your-feature
   ```

4. Open a Pull Request on GitHub **into `dev`** (the repository's default branch is
   `main`, so change the base branch in the PR form). `main` only receives releases.

### Pull Request Guidelines

- Provide a clear description of changes
- Reference any related issues
- Ensure CI tests pass
- Request review from maintainers
- Work through the checklist of the PR template (`.github/PULL_REQUEST_TEMPLATE.md`):
  module conventions, new parameters in `nextflow_schema.json` and `docs/parameters.md`
  ([Adding a parameter](#adding-a-parameter)), and a `CHANGELOG.md` line under
  `## [Unreleased]` for every change a user can see (behaviour, outputs, parameters,
  messages)

## Coding Standards

### Nextflow/Groovy

- Use DSL2 syntax
- **One process per module file**, `modules/local/<name>/main.nf` with its `meta.yml`
  beside it, except for closely related processes that share code
  (`modules/local/format_db/main.nf` holds every `FORMAT_*` / `DOWNLOAD_*` database
  process, `modules/local/bwamem3_samtools/main.nf` holds `BWAMEM3_INDEX` and
  `BWAMEM3_SAMTOOLS`)
- Include `meta.yml` for new modules
- Add `versions.yml` output for software tracking
- Add `stub` blocks for dry-run testing
- **A module's tool flags and publish rules go in `config/modules.config`**, never in the
  process: a `withName: 'TOOL_NAME'` block sets `ext.args` (built from params there) and
  `publishDir`, and the process reads `task.ext.args` (see the template under
  [Process Structure](#process-structure)). A process without a block publishes through
  the default rule of `conf/base.config` (`<output>/<first word of the process name>/`)
- **nf-core modules are installed, not written or edited**: `modules/nf-core/` and
  `modules.json` are managed with nf-core/tools 3.5.1 (`.nf-core.yml`): `nf-core modules
  install` and `update`, and a local change is kept as a `.diff` patch
  (`nf-core modules patch`), as for `pyrodigal` and `bakta/bakta`. `nf-core modules list
  local` lists what is installed (checked with the `nfcore/tools:3.5.1` image).
  `nf-core modules lint <module>` runs but is not a CI gate: on `fastp` it reports 32
  passed, 3 warnings and 1 failure in the upstream `meta.yml` (author citation style), and
  prints a prettier pre-commit warning. `nf-core pipelines lint` is not used: this
  pipeline does not follow the nf-core template, and CI runs no nf-core command. Import an nf-core module with an alias,
  `include { FASTP as NFCORE_FASTP } from '../../modules/nf-core/fastp/main'`, so a
  `withName:` selector and the trace show it is vendored; the two taxonomy modules keep
  their tool names (`KRAKEN2_KRAKEN2 as KRAKEN2`, `BRACKEN_BRACKEN as BRACKEN`)
- **Read booleans through `flagOn()` and numbers with `as Integer`.** Read boolean
  parameters through `flagOn(params.x)` (include it from
  `subworkflows/local/utils_params.nf`) in `main.nf` and subworkflows, or inline
  `params.x.toString().toBoolean()` in modules and config files: on Nextflow >= 26.04 a
  command-line value arrives as a String, so `--x false` would otherwise be the truthy
  `"false"`. Convert numeric parameters with `as Integer` before comparing them or doing
  arithmetic; interpolating them into a command line needs no conversion
- **Give a task that consumes many tasks' outputs a fixed input order**: use
  `.collect(sort: true)` for lists of files and `.toSortedList { a, b -> a[0].id <=> b[0].id }`
  for `[meta, files]` tuples, and sort the grouped lists after `groupTuple`. Plain
  `.collect()`, `.toList()` and `groupTuple` follow task completion order, which changes the
  task's cache key between runs and makes `-resume` re-run it. (`sort: true` orders by work
  directory path, which is stable while the upstream tasks are cached, not by sample name.)
- **Every `errorStrategy` starts with the terminate-on-unschedulable prefix.** In
  `nextflow.config`, `config/modules.config` and `conf/` it must start
  with `task.exitStatus == Integer.MAX_VALUE && task.executor == 'local' ? 'terminate' : ...`
  (then your own rule, e.g. `'retry'`). A local task with no exit status asked for more CPUs
  or memory than the local executor's pool and can never be scheduled (a task killed while
  running reports its signal status, e.g. 137, instead); with `'finish'`, `'retry'` or
  `'ignore'` instead, the run can hang forever waiting for it. A rule that retries listed
  exit codes also retries the cloud preemption codes,
  `task.exitStatus == ['google-batch': 14, 'azurebatch': 50001][task.executor]` (as in
  `nextflow.config`'s default rule): a `withLabel` or `withName` rule replaces the default
  one under every profile, so a profile cannot add them. The slurm and cloud profiles set
  only `maxRetries` and leave the rule to `nextflow.config`
- **Settings a profile may override go before the `profiles {}` block** of
  `nextflow.config` (the default `process` rules, `conf/base.config`): with the strict
  config parser a later block wins, and placed after the profiles they replaced the slurm
  and cloud settings. The resource labels stay after it (under the legacy parser a profile's
  `withLabel` block placed later would replace the whole label block), and so do
  `config/databases.config` and `config/modules.config`. `tests/bin/test_profile_config.sh`
  checks the resolved config of every profile
- **Check a separate singularity image by hand before pinning it.** CI and nf-test run
  under docker only, so nothing automated exercises the image a process
  gets under singularity/apptainer. A separately built singularity image (a Seqera
  `oras://` variant, a depot.galaxyproject.org sif) can lack tools its docker twin has:
  the Seqera singularity variants of the wget, SUPER-FOCUS and HUMAnN images have no
  `gzip` or `tar`, and those processes now pin the docker URI for every engine (the engine
  converts it at pull time). Before pinning a separate singularity image, check that it has
  every command the script calls, and run the process once under apptainer. Also remember
  that a singularity/apptainer image is read-only: a tool must not write inside its own
  install directory
- **Test the engine with `in ['singularity', 'apptainer']`.** A container directive or
  config setting that picks the singularity image tests
  `workflow.containerEngine in ['singularity', 'apptainer']`, never `== 'singularity'`:
  `-profile apptainer` sets the engine to `apptainer`, and an equality test silently sends
  it the docker URI instead (`tests/bin/test_module_meta.sh` checks modules and configs)
- **A tool that names temp files by PID gets `ext.tmpdir_in_task = true`.** Under
  singularity/apptainer every task has its own PID namespace (PIDs repeat from task
  to task) but shares the host `/tmp`. A tool that names its temp files by PID (tRNAscan-SE's
  `/tmp/tscan<PID>`, inside Bakta) can collide with a task running beside it, so give such a
  process a `TMPDIR` in its task directory: `ext.tmpdir_in_task = true` in its
  `config/modules.config` block (see NFCORE_BAKTA_BAKTA), which `conf/base.config`'s `containerOptions`
  turns into `--env TMPDIR=$PWD`. Not a `beforeScript`, which would replace the `module load` an
  HPC config puts there, and not a `withName` `containerOptions`, which would replace base.config's
- **Upstream bugs are worked around in `bin/`, never by editing an image:**
  `bin/run_microbe_census_py3fix.py` and `bin/emapper_ogs_fix.py` replace a broken (or,
  for MicrobeCensus' base count, slow) function and then run the stock CLI with the same
  arguments, and `bin/ps` stands in for
  the `ps` the eggNOG-mapper `.sif` lacks. A patching wrapper must check that the installed
  code is the exact version it was written for and stop otherwise, so an upstream update
  surfaces as an error instead of a silent double fix. Each one names its removal
  condition (a fixed upstream release); remove it with the re-pin

### Adding a parameter

Validation is strict (`validation.unrecognisedParams = 'error'` in `nextflow.config`): a
parameter the schema does not declare stops every run at launch with an "unrecognised
parameter" error, the contributor's own included. A new parameter needs, in the same PR:

1. Its default in the `params {}` block of `nextflow.config` (`null` for an unset path).
2. An entry in `nextflow_schema.json`, in the group where it belongs, with the same
   default, a type, and an `enum` or `pattern` where the values are fixed.
   `tests/bin/test_schema_params.sh` (CI's lint job) fails when the schema and the config
   do not hold the same parameters or their defaults differ.
3. An entry in `docs/parameters.md` (type, default, description, example), and in the
   `--help` text of `main.nf` if users set it often.
4. Reads that follow the rule above: `flagOn()` for booleans, `as Integer` for numbers.
5. A `CHANGELOG.md` line.

A new database registry map in `config/databases.config` (the `*_ref_db` maps) is a
parameter too, but not one users set: add its name to `validation.defaultIgnoreParams`
in `nextflow.config` instead of the schema, or every run stops at launch.
`test_schema_params.sh` checks that every map-valued param is listed there.

### Software versions

`pipeline_info/software_versions.yml` (read by `RESULTS_SUMMARY`) is built by
`subworkflows/local/software_versions.nf` from two sources, written in one format,
`"PROCESS":` followed by indented `tool: version` lines:

- **Local modules** write a `versions.yml` file (template below) and emit it as
  `versions`; `main.nf` (or the stage's subworkflow) mixes `.out.versions` into
  `ch_versions`. Use `.first()` for per-sample processes.
- **nf-core modules** at current upstream emit one tuple per tool on the `versions`
  topic, `tuple val("${task.process}"), val('tool'), eval('tool --version'), topic: versions`,
  and have no `versions.yml`. `main.nf` reads them with `channel.topic('versions')`;
  do not mix such a module's `.out` into `ch_versions` (it has no file to read). When
  `nf-core modules update` moves a module to the topic, drop its `.out.versions` mixes.
  A module kept on `versions.yml` by a local `.diff` patch (`pyrodigal`, `bakta`) stays
  in `ch_versions`.

The tuples of one process become one block with its tools sorted, identical blocks
from per-sample tasks collapse to one, and the file keeps its format.
`tests/subworkflows/software_versions*.nf.test` cover both sources.

Every tool name a module writes there needs a row in `assets/report/tool_references.tsv`
(display name, kind `tool` / `database` / `helper`, and the citekeys of
`assets/report/references.yml`, each checked in PubMed or Crossref): the report's Methods
name and cite tools from it. `tests/bin/test_tool_references.sh` (run in CI) fails on a key
without a row (it reads versions.yml heredocs, `printf` and Python `f.write` lines, and
nf-core topic tuples); a key built at run time (`${aligner}:`) is listed in that script with
the names it can take.

### Databases

A process that downloads or builds a reference database from the registry
(`config/databases.config`) is stored with `storeDir`, not published. Each database
lives in `<databases_dir>/<database>/<registry key>/`, and a run that finds the
process's declared outputs there skips the process. Set this in `config/modules.config`,
in the same form as the existing blocks:

```groovy
withName: 'FORMAT_TOOL_DB' {
    storeDir = { params.databases_dir != null ? "${params.databases_dir}/${workflow.stubRun ? 'stub/' : ''}tool/${params.tool_db}" : null }
    publishDir = [enabled: false]
}
```

- **Switch `publishDir` off.** Without it, the default rule in `conf/base.config`
  also copies the stored database into `<output>/<first word of the process name>/`
  (`format/`, `download/`, ...), a second full copy in the results directory.

- **The closure returns `null` without `databases_dir`**, which turns storeDir off
  (module-level nf-tests run without `--output`). The `stub/` prefix keeps `-stub`
  placeholders away from the folders real runs use.
- **The key is the registry key**, plus anything else that changes the files: the
  aligner for SUPER-FOCUS, `_dense` for the bwa-mem3 index, `_wildcard_k<kmer>` for RGI.
  A database with no upstream version is keyed by the tool version that downloads it
  (`deeparg/1.0.4`); change the key together with the image.
- **Changing a registry `url` needs a new key.** Under the same key, every existing
  `databases_dir` keeps the old download.
- **Declare outputs with fixed names** (`path("k2_ref_db")`, not `path("*")`): storeDir
  only checks that the declared outputs exist, so a glob would match anything left in
  the folder. Make the stub create the same names.
- **When a process starts producing new files** (for example an index built once),
  declare them as their own output or change the key. Otherwise a store filled before
  the change still passes the check and is reused without them.
- **A `versions.yml` output is stored with the database** and re-emitted when the process
  is skipped, so `software_versions.yml` reports the versions that built it.
- **A build from the user's own files** (`--custom_host_fasta`, `--custom_rgi_wildcard`,
  a `--custom_sourmash_db` lineages file) has no version to key on. It runs every time
  and keeps a `publishDir` copy in `databases_dir`, for reuse through the matching
  `--custom_*` option. Its publish path takes the same `stub/` prefix as storeDir.
- **Every database has a `--custom_*` parameter** that bypasses the process and uses the
  path in place.
- **A new registry map goes in `validation.defaultIgnoreParams`** (`nextflow.config`), and
  its `--<tool>_db` key and `--custom_<tool>_db` path in the schema
  ([Adding a parameter](#adding-a-parameter)).

`tests/main.nf.test` ("Stub run skips a database already stored in databases_dir")
pre-fills a store and checks that the process is skipped.

### Process Structure

A module is `modules/local/<name>/main.nf` (the process, in upper case, named after the
folder) plus `modules/local/<name>/meta.yml`. The template follows the rules above:
container pinned to a version, with a separate singularity image only when one was checked
(otherwise one docker URI for every engine), a resource `label`, `task.ext.args` for
optional flags, `versions.yml`, and a `stub:` block that writes the same file names.

```groovy
process TOOL_NAME {
    tag "${meta.id}"
    label 'process_medium'

    // The galaxyproject sif of the same build under singularity/apptainer, the
    // quay image otherwise. Without a checked sif, pin the docker URI alone:
    // container 'quay.io/biocontainers/tool:1.2.3--h1234567_0'
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tool:1.2.3--h1234567_0' :
        'quay.io/biocontainers/tool:1.2.3--h1234567_0' }"

    input:
    tuple val(meta), path(input_files)

    output:
    tuple val(meta), path("*.out"), emit: results
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    tool_command \\
        ${args} \\
        --input ${input_files} \\
        --output ${prefix}.out \\
        --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tool: \$(tool --version 2>&1 | sed 's/.*version //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tool: 1.2.3
    END_VERSIONS
    """
}
```

Its flags and publish rule go in `config/modules.config`, in a `withName:` block next to
those of the same stage (output folders follow the numbered layout, `01_quality_control/`
to `07_functional_annotation/`):

```groovy
withName: 'TOOL_NAME' {
    ext.args = { [
        "--min-score ${params.tool_min_score}",
        params.tool_fast.toString().toBoolean() ? '--fast' : ''
    ].join(' ').trim() }
    publishDir = [
        path: { "${params.output}/05_arg_prediction/tool/${meta.id}" },
        mode: params.publish_dir_mode,
        pattern: '*.out'
    ]
}
```

Wire the process into `main.nf` or its stage's subworkflow under the parameter that gates
it, mix `.out.versions` into `ch_versions`, and add it to `docs/WORKFLOW_DIAGRAM.md`
([Workflow diagrams](#workflow-diagrams)).

### Module meta.yml

Include for all new modules:

```yaml
name: tool_name
description: Brief description
keywords:
  - keyword1
  - keyword2
tools:
  - toolname:
      description: Tool description
      homepage: https://tool.url
      documentation: https://docs.url
      doi: "10.xxxx/xxxxx"
      licence: ["MIT"]

input:
  - meta:
      type: map
      description: Sample metadata
  - input_file:
      type: file
      description: Input file description
      pattern: "*.ext"

output:
  - results:
      type: file
      description: Output description
      pattern: "*.out"
  - versions:
      type: file
      description: Software versions
      pattern: "versions.yml"

authors:
  - "@your_github_handle"
```

`input:` lists the process inputs in order (a `meta` map plus one entry per `val`/`path`
variable, using its name) and `output:` lists every `emit:` name (a `meta` entry for tuple
outputs is optional). When one `main.nf` defines several processes (`format_db`,
`bwamem3_samtools`), keep a single `meta.yml` with `tools` once at the top and a
`processes:` list instead of top-level `input`/`output`, each item with `name` (the process
name), `description`, `input` and `output`. `tests/bin/test_module_meta.sh` (run in CI)
checks every local module against these rules and that every process has a `tag` and a
`label`.

### Results summary

`RESULTS_SUMMARY` (`bin/results_summary.py` and its `bin/summary_*.py` modules, output
`00_summary/`) documents every table it writes:

- **Every column the summary emits has an entry in `TABLE_DOCS`** (`bin/summary_docs.py`):
  the README, the LLM guide and the report's data dictionary are written from it, and
  `tests/bin/test_results_summary.sh` fails on an undocumented column.
- **Every stage is listed in `bin/summary_catalog.py`** (its params, processes, and what
  the summary takes from it).
- A new stage's run-level table reaches the summary through the module's `stage_reports`
  input: a `case` for its file name in `modules/local/results_summary/main.nf` and a
  `mix()` into `ch_summary_reports` in `main.nf`.

### Workflow diagrams

A process or subworkflow you include in `main.nf` or a subworkflow must also appear in
`docs/WORKFLOW_DIAGRAM.md`: as a box in its stage's diagram (with the parameter that gates
it as the edge label) and as a row of the module reference. When you remove one, remove it
there too. `tests/bin/test_workflow_diagram.sh` (run in CI) fails on an included name the
file does not show, and on a diagram box or table row whose process the code no longer
includes. Render the Mermaid blocks to check they parse, e.g.
`docker run --rm -u $(id -u):$(id -g) -v $PWD/docs:/data minlag/mermaid-cli:11.17.1 -i /data/WORKFLOW_DIAGRAM.md -o /tmp/out.md`.
If the change alters the overview figure (`docs/images/bugbuster_overview.svg`, hand-written),
edit it and run `docs/images/render_overview.sh`.

## Testing

### Local Testing

```bash
# Full test
nextflow run main.nf -profile test,docker

# Stub/dry-run test
nextflow run main.nf -profile test,docker -stub

# Resume failed run
nextflow run main.nf -profile test,docker -resume

# nf-test suite (configured in nf-test.config: testsDir tests/, docker profile)
nf-test test

# Single nf-test file
nf-test test tests/modules/fastp.nf.test

# One test of a file, by the hash nf-test list prints
nf-test list tests/main.nf.test
nf-test test tests/main.nf.test@<hash>

# A shell test suite
tests/bin/test_report_scripts.sh
```

Local runs (nextflow and nf-test alike) leave working files in the repo root: `work/`,
`.nextflow/`, `.nextflow.log*`, `.nf-test/`, `.nf-test.log*` and `test_results/`.
A `-profile test` run also writes the decontamination index it builds from the mock host
into the databases directory, which defaults to `./databases/` beside the output:
`databases/stub/bwamem3_index/` for a `-stub` run, `databases/bwamem3_index/` for a full
run. All of these are gitignored, but they accumulate; clean them with:

```bash
tests/bin/clean_test_artifacts.sh
```

It removes the files above and `databases/stub/`, which only `-stub` runs write. It leaves
real run outputs alone (`results/`, `results_stub/`) and the rest of `databases/`, which may
hold databases you downloaded; the full test run's `databases/bwamem3_index/` is kept too,
because a real run with `--custom_host_fasta` writes the same folder. Delete it by hand if
it only holds the test's index.

Parameter validation uses the `nf-schema` plugin, which Nextflow downloads
automatically on first run. For offline environments, pre-install it once:

```bash
nextflow plugin install nf-schema@2.7.3
```

### Adding a test

New code comes with a test. No test downloads data: every input is committed under
`tests/data/` or written by the test itself.

- **A new module**: a module nf-test, `tests/modules/<name>.nf.test`, that runs the real
  script on a small committed fixture where the tool can run without a large database (a
  stub test otherwise), and checks the output content, not only success.
- **A change to the wiring** (a new gate, stage or output path): a pipeline-level test in
  `tests/main.nf.test` (or the topic files beside it, such as `tests/coassembly_groups.nf.test`),
  usually a `-stub` run that checks which tasks ran and which paths were published, or a
  `-preview` launch for a launch check. Build the input with
  `TestData.localSamplesheet(outputDir, projectDir)` (or `localSamplesheet3` / `N`,
  `samplesheetWithIds`) from `tests/lib/TestData.groovy`, which uses the committed reads in
  `tests/data/reads/`. A helper called from a test's `params {}` block must be a method of
  a `tests/lib` class: a function defined at the top of a `.nf.test` file is not visible
  there. A function of `subworkflows/local/utils_*.nf` can be tested on its own with a
  `nextflow_function` test (`tests/functions/`).
- **A new or changed `bin/` script**: a case in the `tests/bin/test_*.sh` suite of that
  script, run in the image the module pins. A new suite needs a step in
  `.github/workflows/ci.yml`: `tests/bin/test_ci_wiring.sh` fails CI on a
  `tests/bin/test_*.sh` that no step runs, or runs where its failure cannot fail the job.
  For an expected failure, check the message (`expect_err` in the suites), not only the
  exit status.
- **A new fixture folder**: a row in `tests/data/README.md` (generated by a
  `tests/bin/make_*` script, hand-written, or captured from a real run), and a generator
  script for generated data. Nothing there is real sample data.
- **A change to what a stub run publishes**: five test files compare the published paths
  with the lists in `tests/data/stub_paths/` (`TestData.publishedPathsDiff`). Regenerate
  them with `BUGBUSTER_UPDATE_STUB_PATHS=1 nf-test test tests/nonpareil.nf.test
  tests/coassembly_groups.nf.test`, review the diff, then run the other list tests
  (`mag_virus_quality`, `dereplicate`, `multisample`).
- **A new parameter**: see [Adding a parameter](#adding-a-parameter);
  `tests/bin/test_schema_params.sh` checks the schema against the config.
- **A result that depends on the Nextflow version**: nf-test runs on Nextflow 25.10.7 (the
  minimum) and on the newest release; branch on `TestData.nextflowBefore(outputDir, 26, 4)`
  rather than skipping the test.

Show that a new test fails when the code it covers is broken (one line changed and
reverted) before relying on it.

### What CI runs

`.github/workflows/ci.yml`, on every pull request into `dev` or `main` and every push to them:

- **Pipeline tests** on Nextflow 25.10.7 and the newest release: `-profile test,docker`
  in `-stub` and full mode, both with `--input` pointing to a samplesheet of the committed
  reads (so CI downloads no test reads), the resource-limit override check
  (`test_resource_limits.sh`) and each profile's resolved config
  (`test_profile_config.sh`).
- **Lint, config and schema checks**: the commands under [Linting](#linting).
- **nf-test suite** on Nextflow 25.10.7 and the newest release.
- **The `tests/bin` suites** (36 today), in two jobs: the database reformat scripts, and
  the taxonomy, report, summary and wrapper scripts. These include three repository
  consistency checks that a docs-only or metadata-only change can fail:
  `test_module_meta.sh` (every module's `meta.yml` against its process, a `tag` and a
  `label` per process, the engine test), `test_workflow_diagram.sh` (every included
  process in `docs/WORKFLOW_DIAGRAM.md`) and `test_tool_references.sh` (a Methods
  reference for every tool version key).

The K-MARVEL and eggNOG-mapper image workflows (`kmarvel-image.yml`,
`eggnog-image.yml`) build and test their image when a pull request touches its folder
under `containers/` or the workflow file; publishing an image is a manual run.

### Linting

The lint job of CI runs these; run them before opening a PR:

```bash
# Strict-syntax check of the scripts and config (errors fail CI; fix them)
nextflow lint .

# The config parses
nextflow config -flat . > /dev/null

# Database-key enums of the schema against config/databases.config
tests/bin/test_schema_db_enums.sh

# The schema is valid JSON Schema, declares the same params as the config, with
# the same defaults (needs python3 with jsonschema)
tests/bin/test_schema_params.sh

# Every tests/bin/test_*.sh is run by ci.yml where it can fail CI (needs PyYAML)
tests/bin/test_ci_wiring.sh

# Optional: resolve the workflow without running any task
nextflow run main.nf -profile test,docker -preview
```

## Questions?

- Open an issue on GitHub
- Contact: ffuentessantander@gmail.com

Thank you for contributing!
