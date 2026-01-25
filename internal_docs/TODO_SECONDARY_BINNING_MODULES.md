# TODO: Secondary Binning Modules Unification

**Status**: ⏸️ DEFERRED  
**Priority**: Low (modules not currently used in pipeline)  
**Estimated Effort**: 2-3 hours  
**Date Created**: January 22, 2026

---

## Overview

AUTOMETA and VAMB are secondary binning modules that exist in the codebase but are **not currently integrated** into the main pipeline workflow. These modules follow the same duplication pattern as the core binning modules (separate per-sample and co-assembly processes) but were deferred during Phase 2-3 refactoring because:

1. They are not called in `main.nf` or any subworkflow
2. They are not required for current pipeline functionality
3. AUTOMETA has complex multi-step workflows with extensive bash scripting that requires careful handling

---

## Current Status

### AUTOMETA
- **Location**: `modules/local/autometa/main.nf`
- **Processes**: 
  - `AUTOMETA` (lines 1-125)
  - `AUTOMETA_COASSEMBLY` (lines 127-255)
- **Duplication**: ~99% identical (130 lines duplicated)
- **Complexity**: HIGH
  - 15+ autometa command steps
  - Complex bash variable substitutions
  - File manipulation and cleanup logic
  - Diamond BLAST integration
  - Taxonomy classification workflow
- **Used By**: Only referenced in `modules/local/basalt/main.nf` (BASALT also not used)

### VAMB
- **Location**: `modules/local/vamb/main.nf`
- **Processes**:
  - `VAMB` (lines 1-23)
  - `VAMB_COASSEMBLY` (lines 25-44)
- **Duplication**: ~95% identical (15 lines duplicated)
- **Complexity**: LOW (simple single-command process)
- **Used By**: Not referenced anywhere in pipeline

### BASALT
- **Location**: `modules/local/basalt/main.nf`
- **Status**: Also not integrated into pipeline
- **Note**: Uses AUTOMETA bins as input, so should be unified together with AUTOMETA

---

## Unification Plan

### Phase A: VAMB Unification (Easy - 30 minutes)

VAMB is straightforward and should be unified first as a warm-up.

**Current Structure**:
```groovy
process VAMB {
    input:
    tuple val(meta), path(contigs), path(depth)
    
    output:
    tuple val(meta), path("*_vamb_bins"), emit: vamb
    
    script:
    def prefix = "${meta.id}"
    """
    vamb bin default \\
        --fasta ${contigs} \\
        --jgi ${depth} \\
        --outdir ${prefix}_vamb_bins \\
        --minfasta ${params.vamb_min_fasta}
    """
}

process VAMB_COASSEMBLY {
    input:
    path(contigs_and_depth)
    
    output:
    path("*_vamb_bins"), emit: vamb
    
    script:
    def prefix = "co_assembly"
    # Similar logic
}
```

**Unified Structure**:
```groovy
process VAMB {
    tag "${meta.id}"
    label 'process_medium'
    
    container 'quay.io/biocontainers/vamb:3.0.2--py36h91eb985_2'
    
    input:
    tuple val(meta), path(contigs), path(depth)
    
    output:
    tuple val(meta), path("${meta.id}_vamb_bins"), emit: bins
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    prefix = meta.id
    def min_fasta = params.vamb_min_fasta ?: 200000
    
    // Input validation
    if (!meta || !meta.id) {
        error "VAMB: meta.id is required"
    }
    
    """
    set -euo pipefail
    
    vamb bin default \\
        --fasta ${contigs} \\
        --jgi ${depth} \\
        --outdir ${prefix}_vamb_bins \\
        --minfasta ${min_fasta}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vamb: \$(vamb --version 2>&1 | sed 's/vamb //')
    END_VERSIONS
    """
    
    stub:
    """
    mkdir -p ${meta.id}_vamb_bins
    touch ${meta.id}_vamb_bins/bin_1.fa
    touch ${meta.id}_vamb_bins/bin_2.fa
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vamb: 3.0.2
    END_VERSIONS
    """
}
```

**Changes Needed**:
1. Unify VAMB and VAMB_COASSEMBLY into single process
2. Add input validation
3. Add version tracking
4. Add stub section
5. Use `${meta.id}` directly in output paths
6. Update `modules.config` with conditional publishDir

---

### Phase B: AUTOMETA Unification (Complex - 2-3 hours)

AUTOMETA requires careful handling due to complex bash scripting.

**Key Differences**:
- **Per-sample**: Uses single BAM file directly
- **Co-assembly**: Merges multiple BAM files first with `samtools merge`

**Unification Strategy**:

```groovy
process AUTOMETA {
    tag "${meta.id}"
    label 'process_high'
    
    container 'quay.io/biocontainers/autometa:2.2.0--pyh7cba7a3_0'
    
    input:
    tuple val(meta), path(contigs), path(bam)
    path(ncbi_db)
    
    output:
    tuple val(meta), path("${meta.id}_autometa_bins"), emit: bins
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    prefix = meta.id
    def is_coassembly = (meta.id == 'coassembly')
    
    // Input validation
    if (!meta || !meta.id) {
        error "AUTOMETA: meta.id is required"
    }
    
    """
    set -euo pipefail
    
    # Prepare BAM file (merge if co-assembly)
    if [ "${is_coassembly}" = "true" ]; then
        samtools merge ${prefix}_merged.bam ${bam}
        bam_file="${prefix}_merged.bam"
    else
        bam_file="${bam}"
    fi
    
    # Prepare contigs (handle potential list for co-assembly)
    if [ "${is_coassembly}" = "true" ]; then
        contigs_input=\$(ls *.fa | tr '\\n' ',' | sed 's/,\$//')
    else
        contigs_input="${contigs}"
    fi
    
    # Run autometa workflow (15+ steps)
    autometa-length-filter \\
        --assembly \${contigs_input} \\
        --cutoff 1000 \\
        --output-fasta ${prefix}_contigs.filtered.fna \\
        --output-stats ${prefix}_contigs.stats.tsv \\
        --output-gc-content ${prefix}_contigs.gc_content.tsv
    
    autometa-coverage \\
        --assembly \${contigs_input} \\
        --bam \${bam_file} \\
        --out ./${prefix}_contigs_coverages.tsv \\
        --cpus ${task.cpus}
    
    # ... (continue with all autometa steps)
    
    # Rename output directory
    mv autometa_bins ${prefix}_autometa_bins
    
    # Generate versions
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        autometa: \$(autometa --version 2>&1 | head -n1 | sed 's/autometa //')
        diamond: \$(diamond version 2>&1 | head -n1 | sed 's/diamond version //')
        samtools: \$(samtools --version 2>&1 | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """
    
    stub:
    """
    mkdir -p ${meta.id}_autometa_bins
    touch ${meta.id}_autometa_bins/bin_1.fa
    touch ${meta.id}_autometa_bins/bin_2.fa
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        autometa: 2.2.0
        diamond: 2.0.15
        samtools: 1.17
    END_VERSIONS
    """
}
```

**Challenges**:
1. **Bash variable escaping**: Many `$` and backticks need careful handling
2. **File list handling**: Co-assembly mode uses bash command substitution to create file lists
3. **Conditional BAM merging**: Need to merge BAMs only for co-assembly
4. **Cleanup logic**: Complex file cleanup at end needs to be preserved
5. **Testing**: Requires NCBI database and significant compute resources

---

### Phase C: BASALT Unification (Optional)

BASALT is a bin refinement tool that uses AUTOMETA bins. Should be unified if AUTOMETA is integrated.

**Current Status**: Not used in pipeline
**Complexity**: Medium
**Dependencies**: Requires unified AUTOMETA

---

## Integration Plan (When Needed)

If these modules are needed in the future:

### 1. Unify Modules
- [ ] Unify VAMB (30 min)
- [ ] Unify AUTOMETA (2-3 hours)
- [ ] Unify BASALT if needed (1 hour)

### 2. Update Subworkflow
- [ ] Add to `subworkflows/local/binning.nf`
- [ ] Wire up channels correctly
- [ ] Handle optional execution (params)

### 3. Update Configuration
- [ ] Add conditional publishDir in `config/modules.config`
- [ ] Add parameters in `nextflow.config`
- [ ] Update `nextflow_schema.json`

### 4. Update Documentation
- [ ] Update pipeline README
- [ ] Add to usage documentation
- [ ] Update output documentation

### 5. Testing
- [ ] Stub-run tests
- [ ] Small dataset tests
- [ ] Compare with original outputs
- [ ] Performance benchmarking

---

## Technical Notes

### AUTOMETA Bash Escaping Issues

The AUTOMETA module contains complex bash that caused JSON parsing errors during initial unification attempts:

```bash
# Examples of problematic patterns:
contigs=`ls | grep -E '.+?filtered_contigs.+' | tr '\\n' ',' | sed 's/.\$//'`
bam_list=`ls | grep -E '.+?all_reads.bam' | tr '\\n' ' ' | sed 's/.\$//'`
file_name=`echo \$file | sed 's/fna/fa/g'`
```

**Solution**: When unifying, use heredoc or create bash script file to avoid JSON escaping issues:

```groovy
script:
"""
#!/bin/bash
set -euo pipefail

# Complex bash logic here without JSON escaping issues
"""
```

### VAMB Input Requirements

VAMB requires:
- Contigs in FASTA format
- Depth file from `jgi_summarize_bam_contig_depths` (MetaBAT2 format)

Currently, CALCULATE_DEPTH generates this format, so VAMB can use the same input as METABAT2.

---

## Files to Modify (When Implementing)

### Modules
- `modules/local/vamb/main.nf` - Unify VAMB processes
- `modules/local/autometa/main.nf` - Unify AUTOMETA processes
- `modules/local/basalt/main.nf` - Unify BASALT processes (optional)

### Subworkflows
- `subworkflows/local/binning.nf` - Add VAMB/AUTOMETA calls

### Configuration
- `config/modules.config` - Add publishDir configs
- `nextflow.config` - Add parameters
- `nextflow_schema.json` - Add parameter schemas

### Documentation
- `README.md` - Update pipeline description
- `docs/usage.md` - Add usage instructions
- `docs/output.md` - Document outputs

---

## Estimated Metrics

### VAMB Unification
- **Processes reduced**: 2 → 1 (-50%)
- **Lines duplicated**: ~15 lines
- **Complexity**: Low
- **Time**: 30 minutes

### AUTOMETA Unification
- **Processes reduced**: 2 → 1 (-50%)
- **Lines duplicated**: ~130 lines
- **Complexity**: High
- **Time**: 2-3 hours

### Total Impact (if implemented)
- **Processes reduced**: 4 → 2 (-50%)
- **Code duplication eliminated**: ~145 lines
- **Maintenance burden**: Reduced significantly

---

## Decision Criteria for Future Implementation

Implement VAMB/AUTOMETA unification when:

1. **User Request**: Users request these binning tools
2. **Comparative Study**: Need to compare multiple binning algorithms
3. **Specific Use Case**: Certain sample types benefit from these tools
4. **Pipeline Completeness**: Goal to offer comprehensive binning options

**Current Recommendation**: ✅ Keep deferred until needed

---

## References

- VAMB: https://github.com/RasmussenLab/vamb
- AUTOMETA: https://github.com/KwanLab/Autometa
- BASALT: Custom bin refinement wrapper

---

## Completion Checklist (For Future)

When implementing, check off:

- [ ] VAMB module unified
- [ ] VAMB tests passing
- [ ] AUTOMETA module unified
- [ ] AUTOMETA tests passing
- [ ] Subworkflow updated
- [ ] Configuration updated
- [ ] Documentation updated
- [ ] Stub-run tests passing
- [ ] Real data tests passing
- [ ] Performance benchmarked
- [ ] Code reviewed
- [ ] User documentation complete

---

**Last Updated**: January 22, 2026  
**Next Review**: When user requests these features or pipeline expansion is planned
