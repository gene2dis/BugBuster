# Contributing to BugBuster

Thank you for your interest in contributing to BugBuster! This document provides guidelines for contributing to the project.

## Table of Contents

- [Code of Conduct](#code-of-conduct)
- [Getting Started](#getting-started)
- [Development Setup](#development-setup)
- [Making Changes](#making-changes)
- [Submitting Changes](#submitting-changes)
- [Coding Standards](#coding-standards)
- [Testing](#testing)

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

## Development Setup

### Requirements

- Nextflow >= 23.04.0
- Docker or Singularity
- Git

### Running Tests

```bash
# Run with test profile
nextflow run main.nf -profile test,docker

# Dry run (stub mode)
nextflow run main.nf -profile test,docker -stub
```

## Making Changes

### Branch Naming

Use descriptive branch names:
- `feature/add-new-tool` - New features
- `fix/memory-issue` - Bug fixes
- `docs/update-readme` - Documentation
- `refactor/module-cleanup` - Code refactoring

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

1. Create a new branch:
   ```bash
   git checkout -b feature/your-feature
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

4. Open a Pull Request on GitHub

### Pull Request Guidelines

- Provide a clear description of changes
- Reference any related issues
- Ensure CI tests pass
- Request review from maintainers

## Coding Standards

### Nextflow/Groovy

- Use DSL2 syntax
- One process per module file
- Include `meta.yml` for new modules
- Add `versions.yml` output for software tracking
- Add `stub` blocks for dry-run testing

### Process Structure

```groovy
process TOOL_NAME {
    tag "${meta.id}"
    container 'quay.io/biocontainers/tool:version'
    
    label 'process_medium'
    
    input:
    tuple val(meta), path(input_files)
    
    output:
    tuple val(meta), path("*.out"), emit: results
    path "versions.yml"           , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = meta.id
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
    def prefix = meta.id
    """
    touch ${prefix}.out
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tool: 1.0.0
    END_VERSIONS
    """
}
```

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

## Testing

### Local Testing

```bash
# Full test
nextflow run main.nf -profile test,docker

# Stub/dry-run test
nextflow run main.nf -profile test,docker -stub

# Resume failed run
nextflow run main.nf -profile test,docker -resume
```

### Linting

```bash
# Check Nextflow syntax
nextflow run main.nf -profile test,docker -preview
```

## Questions?

- Open an issue on GitHub
- Contact: ffuentessantander@gmail.com

Thank you for contributing!
