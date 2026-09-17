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

Note that Clair3's published image is amd64-only. On Apple Silicon it runs only
under emulation, so the NANOPORE workflow is not usable end to end there; the
profile is still useful for the modules that run in the base container.

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

## Current threshold behavior

The current threshold logic is as follows:

- Clair3 receives `--qual` from `clair3_qual` and `--min_mq` from
  `mapping_quality`.
- BCFtools retains variants when `FORMAT/AF >= af_threshold`. No additional
  variant depth or `FILTER=PASS` condition is applied.
- Consensus coverage is calculated with `samtools depth -J -a` without
  additional mapping-quality or base-quality filters.
- Consensus positions with depth less than or equal to `np_min_depth` are
  masked.

Consequently, a low-depth variant can remain in the filtered VCF while the same
position is masked in the consensus. Each sample directory contains a
`<sample>.nanopore_qc.tsv` file that reports the configured thresholds, variant
counts, depth summary, masked bases, and callable consensus percentage. This
file is descriptive and does not affect pipeline success or filtering.

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
BCFtools, consensus generation, and Nanopore QC:

```bash
cd vfnext
NXF_VER=26.04.6 nf-test test integration_tests/nanopore-truth.nf.test
```
