# BugBuster Pipeline Refactoring Plan

**Date**: January 22, 2026  
**Objective**: Eliminate code duplication and improve pipeline maintainability

---

## Overview

The BugBuster pipeline contains **15 duplicate process pairs** (30 total processes) with ~85% code duplication. This plan outlines a systematic refactoring to unify these processes while maintaining functionality and improving performance.

---

## Phase Structure

### Phase 1: Assembly Modules ✅ COMPLETED

**Modules**: MEGAHIT, BBMAP, BOWTIE2_SAMTOOLS  
**Status**: Already unified  
**Outcome**: 50% code reduction, 20-30% performance improvement

**What was done**:
- Unified duplicate processes with conditional logic
- Added error handling and validation
- Optimized file operations (SAM→BAM piping)
- Updated assembly.nf subworkflow
- Updated config/modules.config with conditional publishDir

---

### Phase 2: Core Binning Modules (HIGH PRIORITY)

**Duration**: 12-16 hours  
**Dependencies**: Phase 1 complete

**Modules to Unify** (7 pairs):
1. METABAT2 / METABAT2_COASSEMBLY
2. SEMIBIN / SEMIBIN_COASSEMBLY
3. COMEBIN / COMEBIN_COASSEMBLY
4. METAWRAP / METAWRAP_COASSEMBLY
5. CHECKM2 / CHECKM2_COASSEMBLY
6. GTDB_TK / GTDB_TK_COASSEMBLY
7. CALCULATE_DEPTH / CALCULATE_DEPTH_COASSEMBLY

**Tasks**:
- [ ] Unify each module pair
- [ ] Update binning.nf subworkflow to use unified processes
- [ ] Update config/modules.config with conditional publishDir
- [ ] Add error handling and validation to all modules
- [ ] Add version tracking
- [ ] Create stub sections for testing

**Deliverables**:
- 7 unified module files
- Updated binning.nf subworkflow
- Updated config for binning modules
- Passing stub-run tests for both modes

**Testing**:
- Stub-run for per-sample mode
- Stub-run for co-assembly mode
- Verify channel structures
- Validate output paths match existing structure

---

### Phase 3: Secondary Binning Modules (MEDIUM PRIORITY)

**Duration**: 8-10 hours  
**Dependencies**: Phase 2 complete

**Modules to Unify** (2 pairs):
1. AUTOMETA / AUTOMETA_COASSEMBLY (complex, 257 lines)
2. VAMB / VAMB_COASSEMBLY (incomplete implementation)

**Tasks**:
- [ ] Unify AUTOMETA (preserve complex multi-step workflow)
- [ ] Unify VAMB (complete implementation if needed)
- [ ] Update any subworkflows using these modules
- [ ] Update config

**Deliverables**:
- 2 unified module files
- Updated configurations
- Passing tests

---

### Phase 4: Annotation Modules (LOW PRIORITY)

**Duration**: 4-6 hours  
**Dependencies**: Phase 3 complete

**Modules to Unify** (3 pairs):
1. DEEPARG_BINS / DEEPARG_CONTIGS
2. PRODIGAL_BINS / PRODIGAL_CONTIGS
3. METACERBERUS_CONTIGS / METACERBERUS_BINS

**Tasks**:
- [ ] Unify each annotation module pair
- [ ] Update relevant subworkflows
- [ ] Update config

**Deliverables**:
- 3 unified module files
- Updated configurations
- Passing tests

---

### Phase 5: Integration Testing (CRITICAL)

**Duration**: 10-12 hours  
**Dependencies**: Phases 2, 3, 4 complete

**Testing Components**:

**5.1 Unit Testing**:
- Test each unified process with stub-run
- Verify per-sample mode works
- Verify co-assembly mode works
- Validate error handling
- Check version outputs

**5.2 Integration Testing**:
- Full per-sample assembly + binning workflow
- Full co-assembly + binning workflow
- Assembly only (no binning)
- Error condition handling

**5.3 Regression Testing**:
- Compare outputs with original pipeline
- Verify file structures match
- Check all expected files generated
- Validate file contents are equivalent

**5.4 Performance Benchmarking**:
- Measure execution time
- Track memory usage
- Monitor disk I/O
- Compare with baseline

**Deliverables**:
- Test results documentation
- Performance benchmarks
- Regression test report
- Issue log and resolutions

---

### Phase 6: Documentation & Finalization (CRITICAL)

**Duration**: 4-6 hours  
**Dependencies**: Phase 5 complete

**6.1 Code Documentation**:
- [ ] Update module headers with concise descriptions
- [ ] Document all parameters
- [ ] Add usage examples where helpful
- [ ] Update inline comments

**6.2 User Documentation**:
- [ ] Update README.md
- [ ] Create migration guide (if breaking changes)
- [ ] Document any parameter changes
- [ ] Update usage examples

**6.3 Technical Documentation**:
- [ ] Architecture overview
- [ ] Process flow diagrams
- [ ] Channel structure documentation
- [ ] Configuration guide

**6.4 Final Validation**:
- [ ] Code review
- [ ] Documentation review
- [ ] Final full pipeline test
- [ ] Create release notes

**Deliverables**:
- Complete documentation set
- Migration guide
- Release notes
- Production-ready pipeline

---

## Expected Outcomes

| Metric | Before | After | Improvement |
|--------|--------|-------|-------------|
| Total processes | 60+ | ~45 | -25% |
| Duplicate code | ~8,000 lines | ~4,000 lines | -50% |
| Config blocks | ~45 | ~30 | -33% |
| Error handling coverage | 40% | 100% | +60% |
| Version tracking | 60% | 100% | +40% |
| Performance | Baseline | +20-30% | Faster |

---

## Unified Process Pattern

All unified processes follow this pattern:

```groovy
process TOOL_NAME {
    tag "${meta.id}"
    label "${ meta.id == 'coassembly' ? 'process_high' : 'process_medium' }"
    
    container "..."
    
    input:
    tuple val(meta), path(inputs)
    
    output:
    tuple val(meta), path(outputs), emit: main
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    prefix = meta.id
    def args = task.ext.args ?: ''
    def is_coassembly = (meta.id == 'coassembly')
    
    """
    set -euo pipefail
    
    # Input validation
    if [ ! -f "input_file" ]; then
        echo "ERROR: Input not found" >&2
        exit 1
    fi
    
    # Conditional logic for mode-specific operations
    if [ "${is_coassembly}" == "true" ]; then
        # Co-assembly specific logic
    else
        # Per-sample specific logic
    fi
    
    # Version tracking
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tool: \$(tool --version)
    END_VERSIONS
    """
    
    stub:
    """
    # Create realistic test outputs
    mkdir -p output_dir
    echo "stub" > output_dir/test_file
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tool: "stub"
    END_VERSIONS
    """
}
```

**Configuration Pattern**:

```groovy
withName: 'TOOL_NAME' {
    ext.args = "tool-specific arguments"
    publishDir = [
        path: { meta.id == 'coassembly' ? 
            "${params.output}/category/coassembly" : 
            "${params.output}/category/per_sample/${meta.id}" },
        mode: params.publish_dir_mode,
        pattern: 'output_pattern'
    ]
}
```

---

## Success Criteria

**Must Have**:
- [ ] All duplicate processes unified
- [ ] All tests passing (unit, integration, regression)
- [ ] Output structure unchanged (backward compatible)
- [ ] Performance maintained or improved
- [ ] Documentation complete

**Should Have**:
- [ ] 50% code reduction achieved
- [ ] 100% error handling coverage
- [ ] All processes have version tracking
- [ ] Stub sections for all processes
- [ ] Migration guide created

**Nice to Have**:
- [ ] Performance improvements >20%
- [ ] Automated testing CI/CD
- [ ] Code coverage reports
- [ ] Performance benchmarks documented

---

## Risk Mitigation

**High Risks**:
- Breaking existing workflows → Comprehensive testing, gradual rollout
- Channel structure mismatches → Careful validation, stub testing
- Config conflicts → Systematic config review

**Rollback Plan**:
```bash
# Tag before starting
git tag pre-refactor-stable

# Work in feature branch
git checkout -b refactor-unified-processes

# If rollback needed
git checkout main
git reset --hard pre-refactor-stable
```

---

## Next Steps

1. **Review this plan** and confirm approach
2. **Begin Phase 2**: Start with METABAT2 unification
3. **Test incrementally**: Each module before moving to next
4. **Follow pattern**: Use Phase 1 (assembly) as template
5. **Maintain consistency**: Same approach for all modules

---

## Detailed Implementation

For detailed code examples and step-by-step implementation guides for each module, see:
**`docs/COMPLETE_REFACTORING_PLAN.md`**
