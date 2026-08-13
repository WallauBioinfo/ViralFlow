#!/usr/bin/env python3

import gzip
import hashlib
from pathlib import Path


SEED = 20260609
READ_LENGTH = 1000
START_STEP = 25
REPLICATES = 3
BOUNDARY_REPLICATES = 25
GAP_START = 15001
GAP_END = 15500
EXPECTED_MASKED_BASES = 588
VARIANTS = [
    ("NC_045512.2", 241, "C", "T"),
    ("NC_045512.2", 3037, "C", "T"),
    ("NC_045512.2", 14408, "C", "T"),
    ("NC_045512.2", 21990, "TTTA", "T"),
    ("NC_045512.2", 23403, "A", "G"),
]


def read_fasta(path):
    lines = path.read_text().splitlines()
    name = lines[0][1:].split()[0]
    sequence = "".join(line.strip() for line in lines[1:])
    return name, sequence


def write_fasta(path, name, sequence):
    with path.open("w") as handle:
        handle.write(f">{name}\n")
        for index in range(0, len(sequence), 70):
            handle.write(sequence[index : index + 70] + "\n")


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def apply_variants(reference):
    truth = reference
    for _chrom, position, ref, alt in sorted(
        VARIANTS, key=lambda variant: variant[1], reverse=True
    ):
        start = position - 1
        observed = truth[start : start + len(ref)]
        if observed != ref:
            raise ValueError(
                f"Reference mismatch at {position}: expected {ref}, found {observed}"
            )
        truth = truth[:start] + alt + truth[start + len(ref) :]
    return truth


def truth_offset(reference_position):
    deleted_bases_before = sum(
        len(ref) - len(alt)
        for _chrom, position, ref, alt in VARIANTS
        if position < reference_position
    )
    return reference_position - 1 - deleted_bases_before


def segment_starts(start, end):
    last_start = end - READ_LENGTH
    starts = list(range(start, last_start + 1, START_STEP))
    if starts[-1] != last_start:
        starts.append(last_start)
    return starts


def write_reads(path, truth):
    segments = [
        (0, GAP_START - 1),
        (truth_offset(GAP_END + 1), len(truth)),
    ]
    read_number = 0
    with path.open("wb") as raw_handle:
        with gzip.GzipFile(
            filename="", mode="wb", fileobj=raw_handle, mtime=0
        ) as gzip_handle:
            for segment_number, (start, end) in enumerate(segments, start=1):
                starts = segment_starts(start, end)
                for start_index, read_start in enumerate(starts):
                    copies = REPLICATES
                    if start_index in (0, len(starts) - 1):
                        copies = BOUNDARY_REPLICATES
                    sequence = truth[read_start : read_start + READ_LENGTH]
                    for copy_number in range(copies):
                        read_number += 1
                        strand = "+"
                        if (read_number + SEED) % 2:
                            sequence_out = reverse_complement(sequence)
                            strand = "-"
                        else:
                            sequence_out = sequence
                        header = (
                            f"@truth_{read_number:05d} "
                            f"segment={segment_number} start={read_start + 1} "
                            f"strand={strand} copy={copy_number + 1}"
                        )
                        record = (
                            f"{header}\n{sequence_out}\n+\n{'I' * len(sequence_out)}\n"
                        )
                        gzip_handle.write(record.encode())
    return read_number


def write_vcf(path, contig, reference_length):
    with path.open("w") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write(f"##contig=<ID={contig},length={reference_length}>\n")
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for chrom, position, ref, alt in VARIANTS:
            handle.write(f"{chrom}\t{position}\t.\t{ref}\t{alt}\t60\tPASS\t.\n")


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main():
    fixture_dir = Path(__file__).resolve().parent
    project_dir = fixture_dir.parents[3]
    reference_path = project_dir.parent / "test_files/nanopore/sars-cov-2.fna"
    contig, reference = read_fasta(reference_path)
    truth = apply_variants(reference)

    reference_fasta = fixture_dir / "reference.fasta"
    truth_fasta = fixture_dir / "truth.fasta"
    truth_vcf = fixture_dir / "truth.vcf"
    reads_fastq = fixture_dir / "reads.fastq.gz"
    write_fasta(reference_fasta, contig, reference)
    write_fasta(truth_fasta, "synthetic_truth", truth)
    write_vcf(truth_vcf, contig, len(reference))
    read_count = write_reads(reads_fastq, truth)

    expected_metrics = fixture_dir / "expected_metrics.tsv"
    callable_bases = len(truth) - EXPECTED_MASKED_BASES
    callable_percent = 100.0 * callable_bases / len(truth)
    expected_metrics.write_text(
        "metric\tvalue\n"
        f"reference_length\t{len(reference)}\n"
        f"consensus_length\t{len(truth)}\n"
        f"masked_bases\t{EXPECTED_MASKED_BASES}\n"
        f"consensus_n_bases\t{EXPECTED_MASKED_BASES}\n"
        f"callable_bases\t{callable_bases}\n"
        f"callable_percent\t{callable_percent:.6f}\n"
        f"raw_variant_count\t{len(VARIANTS)}\n"
        f"af_filtered_variant_count\t{len(VARIANTS)}\n"
        f"read_count\t{read_count}\n"
    )

    checksum_paths = [
        reference_fasta,
        truth_fasta,
        truth_vcf,
        reads_fastq,
        expected_metrics,
    ]
    with (fixture_dir / "SHA256SUMS").open("w") as handle:
        for path in checksum_paths:
            handle.write(f"{sha256(path)}  {path.name}\n")


if __name__ == "__main__":
    main()
