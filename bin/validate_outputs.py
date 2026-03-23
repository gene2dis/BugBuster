#!/usr/bin/env python3
"""
BugBuster Output Validation Script

Validates pipeline outputs to ensure all expected files are present
and contain valid data.

Usage:
    python validate_outputs.py --output /path/to/results [--mode full|qc|taxonomy]
"""

import argparse
import json
import os
import sys
from pathlib import Path


def check_file_exists(filepath, description):
    """Check if a file exists and is not empty."""
    path = Path(filepath)
    if not path.exists():
        return False, f"Missing: {description} ({filepath})"
    if path.stat().st_size == 0:
        return False, f"Empty file: {description} ({filepath})"
    return True, f"OK: {description}"


def check_directory_exists(dirpath, description):
    """Check if a directory exists and is not empty."""
    path = Path(dirpath)
    if not path.exists():
        return False, f"Missing directory: {description} ({dirpath})"
    if not path.is_dir():
        return False, f"Not a directory: {description} ({dirpath})"
    if not any(path.iterdir()):
        return False, f"Empty directory: {description} ({dirpath})"
    return True, f"OK: {description}"


def validate_qc_outputs(output_dir):
    """Validate QC-related outputs."""
    results = []
    
    # Check for pipeline_info directory
    pipeline_info = output_dir / "pipeline_info"
    results.append(check_directory_exists(pipeline_info, "Pipeline info directory"))
    
    # Check for execution reports
    if pipeline_info.exists():
        timeline_files = list(pipeline_info.glob("execution_timeline_*.html"))
        results.append((len(timeline_files) > 0, f"Timeline report: {len(timeline_files)} found"))
        
        report_files = list(pipeline_info.glob("execution_report_*.html"))
        results.append((len(report_files) > 0, f"Execution report: {len(report_files)} found"))
        
        trace_files = list(pipeline_info.glob("execution_trace_*.txt"))
        results.append((len(trace_files) > 0, f"Trace file: {len(trace_files)} found"))
    
    # Check for workflow directory
    workflow_dir = output_dir / "workflow"
    results.append(check_directory_exists(workflow_dir, "Workflow output directory"))
    
    return results


def validate_taxonomy_outputs(output_dir):
    """Validate taxonomy-related outputs."""
    results = []
    
    # Check for taxonomy reports
    workflow_dir = output_dir / "workflow"
    if workflow_dir.exists():
        # Look for kraken or sourmash outputs
        for sample_dir in workflow_dir.iterdir():
            if sample_dir.is_dir():
                kraken_dir = sample_dir / "kraken2"
                sourmash_dir = sample_dir / "sourmash"
                
                if kraken_dir.exists():
                    results.append(check_directory_exists(kraken_dir, f"Kraken2 output for {sample_dir.name}"))
                if sourmash_dir.exists():
                    results.append(check_directory_exists(sourmash_dir, f"Sourmash output for {sample_dir.name}"))
    
    return results


def validate_assembly_outputs(output_dir):
    """Validate assembly-related outputs."""
    results = []
    
    workflow_dir = output_dir / "workflow"
    if workflow_dir.exists():
        for sample_dir in workflow_dir.iterdir():
            if sample_dir.is_dir():
                assembly_dir = sample_dir / "assembly"
                if assembly_dir.exists():
                    contig_files = list(assembly_dir.glob("*_contigs.fa"))
                    results.append((len(contig_files) > 0, f"Contigs for {sample_dir.name}: {len(contig_files)} found"))
    
    return results


def validate_binning_outputs(output_dir):
    """Validate binning-related outputs."""
    results = []
    
    workflow_dir = output_dir / "workflow"
    if workflow_dir.exists():
        for sample_dir in workflow_dir.iterdir():
            if sample_dir.is_dir():
                bins_dir = sample_dir / "bins"
                if bins_dir.exists():
                    results.append(check_directory_exists(bins_dir, f"Bins for {sample_dir.name}"))
                    
                    # Check for quality reports
                    checkm_dir = sample_dir / "checkm2"
                    if checkm_dir.exists():
                        results.append(check_directory_exists(checkm_dir, f"CheckM2 for {sample_dir.name}"))
    
    return results


def main():
    parser = argparse.ArgumentParser(description="Validate BugBuster pipeline outputs")
    parser.add_argument("--output", "-o", required=True, help="Path to output directory")
    parser.add_argument("--mode", "-m", default="full", 
                        choices=["full", "qc", "taxonomy", "assembly", "binning"],
                        help="Validation mode")
    parser.add_argument("--json", "-j", action="store_true", help="Output as JSON")
    
    args = parser.parse_args()
    output_dir = Path(args.output)
    
    if not output_dir.exists():
        print(f"Error: Output directory does not exist: {output_dir}")
        sys.exit(1)
    
    all_results = []
    
    # Always check QC outputs
    all_results.extend(validate_qc_outputs(output_dir))
    
    if args.mode in ["full", "taxonomy"]:
        all_results.extend(validate_taxonomy_outputs(output_dir))
    
    if args.mode in ["full", "assembly"]:
        all_results.extend(validate_assembly_outputs(output_dir))
    
    if args.mode in ["full", "binning"]:
        all_results.extend(validate_binning_outputs(output_dir))
    
    # Output results
    passed = sum(1 for ok, _ in all_results if ok)
    failed = sum(1 for ok, _ in all_results if not ok)
    
    if args.json:
        output = {
            "passed": passed,
            "failed": failed,
            "total": len(all_results),
            "results": [{"passed": ok, "message": msg} for ok, msg in all_results]
        }
        print(json.dumps(output, indent=2))
    else:
        print("=" * 60)
        print("BugBuster Output Validation Report")
        print("=" * 60)
        print(f"Output directory: {output_dir}")
        print(f"Validation mode: {args.mode}")
        print("-" * 60)
        
        for ok, msg in all_results:
            status = "✓" if ok else "✗"
            print(f"  {status} {msg}")
        
        print("-" * 60)
        print(f"Total: {len(all_results)} | Passed: {passed} | Failed: {failed}")
        print("=" * 60)
    
    # Exit with error code if any checks failed
    sys.exit(0 if failed == 0 else 1)


if __name__ == "__main__":
    main()
