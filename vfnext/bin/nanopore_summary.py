#!/usr/bin/env python3

import argparse
import gzip
import statistics


def parse_args():
    parser = argparse.ArgumentParser(
        description="Summarize the current ViralFlow Nanopore threshold behavior."
    )
    parser.add_argument("--raw-vcf", required=True)
    parser.add_argument("--filtered-vcf", required=True)
    parser.add_argument("--consensus", required=True)
    parser.add_argument("--coverage", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--clair3-qual", required=True)
    parser.add_argument("--mapping-quality", required=True)
    parser.add_argument("--af-threshold", required=True)
    parser.add_argument("--min-depth", required=True, type=int)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def count_fasta(path):
    length = 0
    n_bases = 0
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                continue
            sequence = line.strip().upper()
            length += len(sequence)
            n_bases += sequence.count("N")
    return length, n_bases


def count_variants(path, min_depth=None):
    total = 0
    depth_counts = {"below": 0, "equal": 0, "above": 0, "missing": 0}
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            total += 1
            if min_depth is None:
                continue
            fields = line.rstrip().split("\t")
            format_keys = fields[8].split(":")
            sample_values = fields[9].split(":")
            sample = dict(zip(format_keys, sample_values))
            depth_value = sample.get("DP", ".")
            if depth_value in ("", "."):
                depth_counts["missing"] += 1
                continue
            depth = int(depth_value)
            if depth < min_depth:
                depth_counts["below"] += 1
            elif depth == min_depth:
                depth_counts["equal"] += 1
            else:
                depth_counts["above"] += 1
    return total, depth_counts


def read_depths(path):
    depths = []
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            if len(fields) >= 3:
                depths.append(int(fields[2]))
    return depths


def main():
    args = parse_args()
    depths = read_depths(args.coverage)
    raw_variants, _ = count_variants(args.raw_vcf)
    filtered_variants, variant_depths = count_variants(
        args.filtered_vcf, min_depth=args.min_depth
    )
    reference_length, _ = count_fasta(args.reference)
    consensus_length, consensus_n_bases = count_fasta(args.consensus)
    masked_bases = sum(depth <= args.min_depth for depth in depths)
    zero_depth_bases = sum(depth == 0 for depth in depths)
    callable_bases = consensus_length - consensus_n_bases
    callable_percent = (
        100.0 * callable_bases / consensus_length if consensus_length else 0.0
    )

    metrics = [
        ("clair3_qual", args.clair3_qual),
        ("mapping_quality", args.mapping_quality),
        ("af_threshold", args.af_threshold),
        ("consensus_min_depth", args.min_depth),
        ("consensus_mask_rule", "depth <= consensus_min_depth"),
        ("raw_variant_count", raw_variants),
        ("af_filtered_variant_count", filtered_variants),
        ("filtered_variants_depth_below_min", variant_depths["below"]),
        ("filtered_variants_depth_equal_min", variant_depths["equal"]),
        ("filtered_variants_depth_above_min", variant_depths["above"]),
        ("filtered_variants_depth_missing", variant_depths["missing"]),
        ("reference_length", reference_length),
        ("mean_depth", f"{statistics.mean(depths):.6f}" if depths else "0.000000"),
        (
            "median_depth",
            f"{statistics.median(depths):.6f}" if depths else "0.000000",
        ),
        ("zero_depth_bases", zero_depth_bases),
        ("masked_bases", masked_bases),
        ("consensus_length", consensus_length),
        ("consensus_n_bases", consensus_n_bases),
        ("callable_bases", callable_bases),
        ("callable_percent", f"{callable_percent:.6f}"),
    ]

    with open(args.output, "w") as handle:
        handle.write("metric\tvalue\n")
        for metric, value in metrics:
            handle.write(f"{metric}\t{value}\n")


if __name__ == "__main__":
    main()
