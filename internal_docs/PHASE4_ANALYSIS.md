# Phase 4: Performance Optimization Analysis

## Objective
Optimize CPU allocation and parallelization for slow-running processes, particularly mapping and annotation steps.

## Current Resource Allocation

### Process Labels (from `nextflow.config`)

| Label | CPUs | Memory | Time | Use Case |
|-------|------|--------|------|----------|
| `process_single` | 1 | 4 GB | 4h | Simple tasks, reporting |
| `process_low` | 4 | 12 GB | 8h | Light computation |
| `process_medium` | 8 | 36 GB | 24h | Moderate computation |
| `process_high` | 16 | 72 GB | 72h | Heavy computation |
| `process_high_memory` | - | 256 GB | - | Memory-intensive |

**Max resources**: 16 CPUs, 128 GB RAM, 240h

## Performance Bottlenecks Identified

### 1. Mapping Steps (Bowtie2)

#### BOWTIE2 (Host/PhiX decontamination)
- **Current**: `process_medium` (8 CPUs)
- **Usage**: Aligns all samples against host/PhiX databases
- **Issue**: Runs sequentially per sample, not utilizing full CPU capacity
- **Recommendation**: Increase to `process_high` (16 CPUs) - alignment is highly parallelizable

#### BOWTIE2_SAMTOOLS (Contig alignment)
- **Current**: `process_medium` (8 CPUs)
- **Usage**: Aligns reads back to assembled contigs for coverage calculation
- **Issue**: Critical bottleneck - runs for every sample AND co-assembly
- **Recommendation**: Increase to `process_high` (16 CPUs)
- **Rationale**: 
  - Bowtie2 alignment scales well with threads
  - Samtools sort also benefits from multiple threads
  - This is a critical path process affecting binning quality

#### BOWTIE2_SAMTOOLS_DEPTH (Bin alignment)
- **Current**: `process_medium` (8 CPUs)
- **Usage**: Aligns reads to individual bins
- **Issue**: Loops through bins sequentially, not parallelized
- **Recommendation**: Keep at `process_medium` but optimize script to parallelize bin processing

### 2. Annotation Modules

#### PRODIGAL_BINS
- **Current**: `process_medium` (8 CPUs)
- **Usage**: Gene prediction on bins
- **Current implementation**: Manual parallelization with bash loops
- **Issue**: Inefficient parallelization, waits every N iterations
- **Recommendation**: 
  - Keep `process_medium` (appropriate for gene prediction)
  - Improve parallelization logic (use GNU parallel or better bash implementation)

#### DEEPARG_BINS
- **Current**: `process_low` (4 CPUs)
- **Usage**: ARG prediction on bin proteins
- **Issue**: Loops through bins sequentially, no parallelization
- **Recommendation**: 
  - Increase to `process_medium` (8 CPUs)
  - Add parallelization to process multiple bins simultaneously

#### DEEPARG_CONTIGS
- **Current**: `process_low` (4 CPUs)
- **Usage**: ARG prediction on contig proteins
- **Recommendation**: Increase to `process_medium` (8 CPUs)

#### METACERBERUS_CONTIGS & METACERBERUS_BINS
- **Current**: `process_high` (16 CPUs)
- **Usage**: Functional annotation
- **Status**: Already optimized, no changes needed

### 3. Binning Modules

#### COMEBIN
- **Current**: `process_high` (16 CPUs)
- **Status**: Appropriate, no changes needed

#### METAWRAP (Bin refinement)
- **Current**: `process_high` (16 CPUs)
- **Status**: Appropriate, no changes needed

#### SEMIBIN
- **Current**: `process_medium` (8 CPUs)
- **Recommendation**: Consider increasing to `process_high` (16 CPUs) for faster binning

#### GTDB-TK (Taxonomy)
- **Current**: `process_high` + `process_high_memory` (16 CPUs, 256 GB)
- **maxForks**: 2 (limits parallelism to avoid memory issues)
- **Status**: Appropriate, no changes needed

### 4. Assembly

#### MEGAHIT
- **Current**: `process_medium` (8 CPUs)
- **Recommendation**: Increase to `process_high` (16 CPUs)
- **Rationale**: Assembly is CPU-intensive and benefits from more threads

### 5. Other High-Resource Processes

#### AUTOMETA
- **Current**: `process_high` (16 CPUs)
- **Status**: Appropriate, no changes needed

#### VAMB_COASSEMBLY
- **Current**: `process_high` (16 CPUs)
- **Status**: Appropriate, no changes needed

#### CLUSTERING (MMseqs2)
- **Current**: `process_high` (16 CPUs)
- **Status**: Appropriate, no changes needed

## Optimization Strategy

### Priority 1: Critical Path Optimizations (High Impact)

1. **BOWTIE2_SAMTOOLS**: `process_medium` → `process_high`
   - Impact: HIGH (critical for binning quality)
   - Effort: LOW (label change only)

2. **BOWTIE2**: `process_medium` → `process_high`
   - Impact: HIGH (affects all samples)
   - Effort: LOW (label change only)

3. **MEGAHIT**: `process_medium` → `process_high`
   - Impact: HIGH (assembly is often the slowest step)
   - Effort: LOW (label change only)

### Priority 2: Annotation Optimizations (Medium Impact)

4. **DEEPARG_BINS**: Add parallelization + increase to `process_medium`
   - Impact: MEDIUM (runs on all bins)
   - Effort: MEDIUM (script modification needed)

5. **DEEPARG_CONTIGS**: `process_low` → `process_medium`
   - Impact: MEDIUM (runs per sample)
   - Effort: LOW (label change only)

6. **PRODIGAL_BINS**: Improve parallelization logic
   - Impact: MEDIUM (gene prediction is moderately fast)
   - Effort: MEDIUM (script optimization)

### Priority 3: Additional Optimizations (Lower Impact)

7. **SEMIBIN**: `process_medium` → `process_high`
   - Impact: LOW-MEDIUM (binning step)
   - Effort: LOW (label change only)

8. **BOWTIE2_SAMTOOLS_DEPTH**: Optimize bin loop parallelization
   - Impact: LOW (only if using bin-level depth)
   - Effort: MEDIUM (script modification)

## Expected Performance Improvements

### Time Savings Estimates (per sample)

| Process | Current Time | Optimized Time | Savings |
|---------|-------------|----------------|---------|
| BOWTIE2 | ~2-4h | ~1-2h | 50% |
| BOWTIE2_SAMTOOLS | ~4-8h | ~2-4h | 50% |
| MEGAHIT | ~6-12h | ~3-6h | 50% |
| DEEPARG_BINS | ~2-4h | ~1-2h | 50% |
| DEEPARG_CONTIGS | ~1-2h | ~30m-1h | 50% |

**Total potential savings**: 8-15 hours per sample (30-40% reduction in total runtime)

## Implementation Plan

1. Update process labels for simple optimizations (Priority 1 & some Priority 2)
2. Optimize DEEPARG_BINS parallelization
3. Improve PRODIGAL_BINS parallelization
4. Test with small dataset to verify improvements
5. Document changes and performance metrics

## Notes

- All optimizations respect the `max_cpus = 16` limit
- Memory allocations remain appropriate for each process
- No changes to core algorithms or outputs
- Optimizations focus on better CPU utilization
