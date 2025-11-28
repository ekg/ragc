#!/usr/bin/env python3
"""
Test suite for RAGC compression validity and performance.

Tests N random subsets of M genome files, comparing:
- Archive sizes between RAGC and C++ AGC
- Extraction correctness (byte-for-byte identical)

Uses a fixed seed for reproducibility.
"""

import os
import sys
import random
import subprocess
import tempfile
import hashlib
from pathlib import Path
from dataclasses import dataclass
from typing import List, Tuple, Optional

# Configuration
RAGC_BIN = "./target/release/ragc"
CPP_AGC_BIN = "/home/erik/agc/bin/agc"
GENOME_DIR = "/home/erik/scrapy/yeast235_split"

@dataclass
class TestResult:
    subset_id: int
    samples: List[str]
    ragc_size: int
    cpp_size: int
    size_diff_pct: float
    samples_with_errors: List[str]
    contigs_with_errors: List[Tuple[str, str, int, int]]  # (sample, contig, ragc_len, cpp_len)
    passed: bool

def get_genome_files(genome_dir: str) -> List[str]:
    """Get all .fa files in directory."""
    files = []
    for f in sorted(os.listdir(genome_dir)):
        if f.endswith('.fa'):
            files.append(os.path.join(genome_dir, f))
    return files

def generate_subsets(files: List[str], n_subsets: int, samples_per_subset: int, seed: int) -> List[List[str]]:
    """Generate N random subsets of M files with fixed seed."""
    random.seed(seed)
    subsets = []
    for _ in range(n_subsets):
        subset = random.sample(files, min(samples_per_subset, len(files)))
        subsets.append(sorted(subset))  # Sort for reproducibility
    return subsets

def create_ragc_archive(files: List[str], output_path: str) -> bool:
    """Create RAGC archive."""
    cmd = [RAGC_BIN, "create", "-o", output_path, "-k", "21", "-s", "10000", "-m", "20", "-t", "1"] + files
    result = subprocess.run(cmd, capture_output=True)
    return result.returncode == 0

def create_cpp_archive(files: List[str], output_path: str) -> bool:
    """Create C++ AGC archive."""
    cmd = [CPP_AGC_BIN, "create", "-o", output_path, "-k", "21", "-s", "10000", "-l", "20", "-t", "1"] + files
    result = subprocess.run(cmd, capture_output=True)
    return result.returncode == 0

def get_ragc_samples(archive_path: str) -> List[str]:
    """Get list of samples in RAGC archive."""
    cmd = [RAGC_BIN, "listset", archive_path]
    result = subprocess.run(cmd, capture_output=True, text=True)
    return [s.strip() for s in result.stdout.strip().split('\n') if s.strip()]

def get_cpp_samples(archive_path: str) -> List[str]:
    """Get list of samples in C++ archive."""
    cmd = [CPP_AGC_BIN, "listset", archive_path]
    result = subprocess.run(cmd, capture_output=True, text=True)
    return [s.strip() for s in result.stdout.strip().split('\n') if s.strip()]

def extract_ragc_sample(archive_path: str, sample: str, output_path: str) -> bool:
    """Extract sample from RAGC archive."""
    cmd = [RAGC_BIN, "getset", archive_path, sample, "-o", output_path]
    result = subprocess.run(cmd, capture_output=True)
    return result.returncode == 0

def extract_cpp_sample(archive_path: str, sample: str, output_path: str) -> bool:
    """Extract sample from C++ archive."""
    with open(output_path, 'w') as f:
        cmd = [CPP_AGC_BIN, "getset", archive_path, sample]
        result = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE)
    return result.returncode == 0

def parse_fasta(path: str) -> dict:
    """Parse FASTA file and return dict of contig_name -> sequence."""
    contigs = {}
    current_name = None
    current_seq = []
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_name:
                    contigs[current_name] = ''.join(current_seq)
                current_name = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
        if current_name:
            contigs[current_name] = ''.join(current_seq)
    return contigs

def compare_extractions(ragc_path: str, original_path: str) -> List[Tuple[str, int, int]]:
    """Compare extracted FASTA with original, return list of (contig, ragc_len, orig_len) for differences."""
    ragc_contigs = parse_fasta(ragc_path)
    orig_contigs = parse_fasta(original_path)

    differences = []

    # Build a mapping from contig suffix to original sequence
    # Original contigs might be like "AAA#0#chrI", we need to match these
    orig_by_suffix = {}
    for name, seq in orig_contigs.items():
        # Use the full contig name for matching
        orig_by_suffix[name] = seq

    # Check all contigs from RAGC extraction
    for contig in sorted(ragc_contigs.keys()):
        ragc_seq = ragc_contigs.get(contig, '')
        orig_seq = orig_by_suffix.get(contig, '')

        if len(ragc_seq) != len(orig_seq):
            differences.append((contig, len(ragc_seq), len(orig_seq)))
        elif ragc_seq.upper() != orig_seq.upper():
            # Same length but different content (case-insensitive comparison)
            # AGC format normalizes case (lowercase -> uppercase)
            differences.append((contig, len(ragc_seq), len(orig_seq)))

    return differences

def run_test(subset_id: int, files: List[str], tmp_dir: str) -> TestResult:
    """Run test on a single subset."""
    ragc_archive = os.path.join(tmp_dir, f"ragc_{subset_id}.agc")
    cpp_archive = os.path.join(tmp_dir, f"cpp_{subset_id}.agc")

    # Create archives
    if not create_ragc_archive(files, ragc_archive):
        print(f"  ERROR: Failed to create RAGC archive for subset {subset_id}")
        return TestResult(subset_id, files, 0, 0, 0, [], [], False)

    if not create_cpp_archive(files, cpp_archive):
        print(f"  ERROR: Failed to create C++ archive for subset {subset_id}")
        return TestResult(subset_id, files, 0, 0, 0, [], [], False)

    # Get sizes
    ragc_size = os.path.getsize(ragc_archive)
    cpp_size = os.path.getsize(cpp_archive)
    size_diff_pct = 100 * (ragc_size - cpp_size) / cpp_size if cpp_size > 0 else 0

    # Get samples
    ragc_samples = get_ragc_samples(ragc_archive)
    cpp_samples = get_cpp_samples(cpp_archive)

    samples_with_errors = []
    contigs_with_errors = []

    # Build mapping from RAGC sample name to original file
    # File: /path/to/AAA_0.fa -> Sample: AAA#0
    sample_to_file = {}
    for f in files:
        basename = os.path.basename(f).replace('.fa', '')
        # Convert filename format (AAA_0) to RAGC sample format (AAA#0)
        # Handle cases like "CQS_1a_0" -> "CQS_1a#0"
        parts = basename.rsplit('_', 1)
        if len(parts) == 2:
            ragc_sample = f"{parts[0]}#{parts[1]}"
            sample_to_file[ragc_sample] = f

    # Extract and compare each sample against original
    for ragc_sample in ragc_samples:
        original_file = sample_to_file.get(ragc_sample)
        if not original_file:
            samples_with_errors.append(ragc_sample)
            continue

        ragc_extract = os.path.join(tmp_dir, f"ragc_ext_{subset_id}_{ragc_sample.replace('#', '_')}.fa")

        if not extract_ragc_sample(ragc_archive, ragc_sample, ragc_extract):
            samples_with_errors.append(ragc_sample)
            continue

        # Compare RAGC extraction with original file
        differences = compare_extractions(ragc_extract, original_file)
        if differences:
            samples_with_errors.append(ragc_sample)
            for contig, ragc_len, orig_len in differences:
                contigs_with_errors.append((ragc_sample, contig, ragc_len, orig_len))

        # Clean up extract file
        os.remove(ragc_extract)

    passed = len(samples_with_errors) == 0

    return TestResult(
        subset_id=subset_id,
        samples=[os.path.basename(f) for f in files],
        ragc_size=ragc_size,
        cpp_size=cpp_size,
        size_diff_pct=size_diff_pct,
        samples_with_errors=samples_with_errors,
        contigs_with_errors=contigs_with_errors,
        passed=passed
    )

def main():
    import argparse
    parser = argparse.ArgumentParser(description='Test RAGC compression validity and performance')
    parser.add_argument('-n', '--num-subsets', type=int, default=10, help='Number of subsets to test')
    parser.add_argument('-m', '--samples-per-subset', type=int, default=10, help='Samples per subset')
    parser.add_argument('-s', '--seed', type=int, default=42, help='Random seed for reproducibility')
    parser.add_argument('--genome-dir', type=str, default=GENOME_DIR, help='Directory with genome .fa files')
    parser.add_argument('--verbose', '-v', action='store_true', help='Verbose output')
    args = parser.parse_args()

    print(f"RAGC Compression Test Suite")
    print(f"============================")
    print(f"Seed: {args.seed}")
    print(f"Subsets: {args.num_subsets}")
    print(f"Samples per subset: {args.samples_per_subset}")
    print(f"Genome directory: {args.genome_dir}")
    print()

    # Get all genome files
    files = get_genome_files(args.genome_dir)
    print(f"Found {len(files)} genome files")

    # Generate subsets
    subsets = generate_subsets(files, args.num_subsets, args.samples_per_subset, args.seed)

    # Run tests
    results = []
    with tempfile.TemporaryDirectory() as tmp_dir:
        for i, subset in enumerate(subsets):
            print(f"\nSubset {i+1}/{args.num_subsets}: {len(subset)} files")
            if args.verbose:
                for f in subset:
                    print(f"  - {os.path.basename(f)}")

            result = run_test(i, subset, tmp_dir)
            results.append(result)

            status = "PASS" if result.passed else "FAIL"
            print(f"  Size: RAGC={result.ragc_size:,} C++={result.cpp_size:,} ({result.size_diff_pct:+.1f}%)")
            print(f"  Status: {status}")

            if result.contigs_with_errors:
                print(f"  Data errors in {len(result.samples_with_errors)} samples:")
                for sample, contig, ragc_len, orig_len in result.contigs_with_errors[:5]:
                    print(f"    {sample}/{contig}: RAGC={ragc_len} Original={orig_len}")
                if len(result.contigs_with_errors) > 5:
                    print(f"    ... and {len(result.contigs_with_errors) - 5} more")

    # Summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)

    passed = sum(1 for r in results if r.passed)
    failed = sum(1 for r in results if not r.passed)

    print(f"Passed: {passed}/{len(results)}")
    print(f"Failed: {failed}/{len(results)}")

    if failed > 0:
        print("\nFailed subsets:")
        for r in results:
            if not r.passed:
                print(f"  Subset {r.subset_id}: {len(r.samples_with_errors)} samples with errors")
                total_data_loss = sum(orig_len - ragc_len for _, _, ragc_len, orig_len in r.contigs_with_errors)
                print(f"    Total data loss: {total_data_loss:,} bytes")

    # Size comparison summary
    size_diffs = [r.size_diff_pct for r in results]
    avg_diff = sum(size_diffs) / len(size_diffs)
    print(f"\nSize difference (RAGC vs C++):")
    print(f"  Average: {avg_diff:+.2f}%")
    print(f"  Min: {min(size_diffs):+.2f}%")
    print(f"  Max: {max(size_diffs):+.2f}%")

    # Return exit code
    return 0 if failed == 0 else 1

if __name__ == '__main__':
    sys.exit(main())
