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
### Aptainner setup

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