# Nanopore

This document provides guidelines and instructions on how to run nanopore data on ViralFlow.

> this is a provisory documentation for development purpose, the oficial documentation should be included on the regular ViralFlow documentation.

---

## Build container

```
cd /.../ViralFlow/vfnext/containers/
singularity build baseContainer.sif Nanopore_baseContainer.sing
# or
apptainer build baseContainer.sif Nanopore_baseContainer.sing
```
### Docker

Where Singularity and Apptainer are not available — CI runners, and developer
machines such as Apple Silicon Macs — build the equivalent Docker image and use
the `docker` profile:

```bash
cd /.../ViralFlow/vfnext/containers/
docker build -f nanopore_base.Dockerfile -t viralflow/nanopore-base:2.0.0a1 .
```

```bash
nextflow run /../ViralFlow/vfnext/main.nf \
        -profile docker \
        --mode NANOPORE \
        --inDir /path/to/np_input_dir/ \
        --referenceGenome /path/to/reference.fna \
        -resume
```

`nanopore_base.Dockerfile` mirrors `Nanopore_baseContainer.sing` and pins the
same tool versions; `tests/test_container_recipes.py` fails if the two drift
apart. The profile sets `params.base_container` to the image above, so the
image tag must track the pipeline version.

Clair3's published image is amd64 only, so on Apple Silicon it runs under
emulation. The workflow still completes end to end, the truth test included,
but expect Clair3 to take a couple of minutes where it would take seconds
natively.

### Apptainer setup

Unfortunately, apptainer does not support `library://` protocol. To make it work on this protocol run the following commands:

```bash
apptainer remote add --no-login SylabsCloud cloud.sycloud.io
apptainer remote use SylabsCloud
apptainer remote list
```

---

## Run pipeline

To run the nanopore mode of the viral pipeline simply run:

```bash
nextflow run /../ViralFlow/vfnext/main.nf \
        --mode NANOPORE \
        --inDir /path/to/np_input_dir/ \
        --referenceGenome /path/to/reference.fna \
        -resume
```

to run it using apptainer, just add `-profile apptainer` to your nextflow command.

## Primer clipping (optional, off by default)

Primer clipping does not run unless you ask for it. Supplying `--primersBED`
turns it on, exactly as in ILLUMINA mode:

```bash
nextflow run /../ViralFlow/vfnext/main.nf \
        --mode NANOPORE \
        --inDir /path/to/np_input_dir/ \
        --referenceGenome /path/to/reference.fna \
        --primersBED /path/to/primers.bed \
        -resume
```

When enabled, `samtools ampliconclip --strand --hard-clip` runs between the
Minimap2 alignment and variant calling, and the clipped BAM replaces the raw
alignment for **everything** downstream — Clair3, the depth used for masking and
the consensus. Clipping only some of those would let the consensus disagree with
the variants it was built from.

Each sample directory gains `<sample>.primer_clip.bam` (plus its index) and
`<sample>.ampliconclip.txt`, samtools' report of what was clipped.

Without `--primersBED` the pipeline logs a warning and leaves the alignment
untouched, so primer-derived bases remain in the consensus. That is the default
because it is the right behaviour for non-amplicon data; if you are working with
amplicon protocols, supply the BED.

## Read-end trimming (optional, off by default)

`--trimLen` trims a fixed number of bases from both ends of every read. It is 0
by default, which disables the step. ILLUMINA applies the same parameter inside
fastp, before alignment; NANOPORE has no equivalent FASTQ step, so it trims the
aligned BAM with `bam trimBam` instead. When both are enabled, trimming runs
after primer clipping.

```bash
nextflow run /../ViralFlow/vfnext/main.nf \
        --mode NANOPORE \
        --inDir /path/to/np_input_dir/ \
        --referenceGenome /path/to/reference.fna \
        --trimLen 5 \
        -resume
```

**How bamUtil trims here.** The step runs `bam trimBam --clip`, which *soft
clips* the trimmed bases: the alignment start moves past them and the CIGAR
gains leading and trailing `S` operations, so the bases stay in the record but
no longer occupy reference positions. They therefore drop out of the depth that
`np_min_depth` masking is computed from, which is what you want - a position
covered only by trimmed bases is correctly seen as uncovered.

`--clip` is not bamUtil's default. Left to itself, `trimBam` *masks* the bases
instead, setting them to `N` with quality `!` while the alignment keeps spanning
the same reference positions. That inflated the consensus depth with bases
supporting no allele, and Clair3's pileup aborted on such reads outright
(`munmap_chunk(): invalid pointer`), so trimming could not be used at all.

A trim wider than a read leaves that read unclipped and marked unmapped, rather
than producing a degenerate alignment. Each sample directory gains
`<sample>.trim.sorted.bam` and its index.

## Current threshold behavior

The current threshold logic is as follows:

- Clair3 receives `--qual` from `clair3_qual` and `--min_mq` from
  `mapping_quality`. Clair3 does not drop a call below `--qual`; it keeps it
  with `FILTER=LowQual`.
- BCFtools retains variants when `FILTER=PASS` and
  `FORMAT/AF >= af_threshold`, so `LowQual` calls stay in Clair3's
  `merge_output.vcf.gz` but are left out of the filtered VCF and the
  consensus. No variant depth condition is applied.
- Consensus coverage is calculated with `samtools depth -J -aa` without
  additional mapping-quality or base-quality filters.
- Consensus positions with depth less than or equal to `np_min_depth` are
  masked.
- Positions no read reached count as depth 0 and are masked like any other,
  including whole contigs with no aligned reads. A sample where nothing aligns,
  such as a negative control, therefore produces an all-`N` consensus with
  `callable_percent` 0, not a copy of the reference.

Consequently, a low-depth variant can remain in the filtered VCF while the same
position is masked in the consensus. Each sample directory contains a
`<sample>.nanopore_summary.tsv` file that reports the configured thresholds,
variant counts, depth summary, masked bases, and callable consensus percentage.
This file is descriptive and does not affect pipeline success or filtering.

### Reading `masked_bases` against `consensus_n_bases`

These two look interchangeable and are equal on an N-free reference, but they
count different things:

| Metric | Counts | Coordinate space |
|---|---|---|
| `masked_bases` | positions ViralFlow masked because `depth <= np_min_depth` | reference |
| `consensus_n_bases` | every `N` in the consensus, whatever its origin | consensus |

`consensus_n_bases - masked_bases` is therefore the N that did **not** come from
your coverage threshold — in practice, ambiguity codes the reference already
carried. With SARS-CoV-2 (`NC_045512.2`) the difference is zero, because that
reference contains no N; with a reference that does, they diverge by exactly
that count.

`callable_bases` and `callable_percent` are derived from `consensus_n_bases`,
since a reference N is no more callable than a low-coverage one.

One related behaviour worth knowing when interpreting these numbers: `bcftools
consensus` does not apply a variant that falls inside a masked region, so a
deletion there leaves the consensus length unchanged.

## Reproducibility metadata

Every run writes reproducibility records under `RUN_METADATA` inside the output
directory:

- `run_manifest.json`: pipeline revision, parameters, runtime context, and final status.
- `input_checksums.tsv`: SHA-256 checksums for reads and reference inputs.
- `software_versions.tsv`: versions of the core tools used by the selected mode.
- `container_manifest.tsv`: container identities and local SIF checksums.
- `execution_trace.tsv`: per-task status and resource usage.
- `execution_report.html` and `execution_timeline.html`: Nextflow execution reports.

Metadata collection is part of the workflow. Missing tools, unreadable inputs, or
an invalid local container path cause the run to fail rather than recording
incomplete provenance.

## Truth integration test

The deterministic FASTQ-to-consensus truth test is intentionally separate from
the regular unit-test suite because it runs Porechop, Minimap2, Clair3,
BCFtools, consensus generation, and the Nanopore summary:

```bash
cd vfnext
NXF_VER=26.04.6 nf-test test integration_tests/nanopore-truth.nf.test
```

`integration_tests/` also holds a multi-sample fan-out test and a negative
control (`nanopore-no-reads.nf.test`): reads that align nowhere, which must come
out as an all-`N` consensus rather than a copy of the reference. Run the whole
directory, as CI does, with `nf-test test integration_tests/`.
