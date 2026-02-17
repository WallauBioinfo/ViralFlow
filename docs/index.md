# ViralFlow Documentation

```{admonition} 🌐 Language / Idioma / Lengua
:class: tip

Esta documentação também está disponível em **Português**: [Clique aqui](https://viralflow.readthedocs.io/pt-br/latest/)

Esta documentación también está disponible en **Español**: [Haga clic aquí](https://viralflow.readthedocs.io/es/latest/)
```

**ViralFlow** is a workflow developed for viral genomic surveillance that performs several genomic analyses based on reference genome assembly.

The current version of the pipeline was written in the Nextflow workflow language, and can be applied to different viruses. Only a single command line is needed to run the entire workflow. To make using ViralFlow more accessible, an installing "wrapper" was written to make the pipeline easier for people unfamiliar with Nextflow.

Currently, the code has been tested only for Illumina paired-end and single-end data. The ViralFlow development team can not provide support for errors/problems of applying ViralFlow on sequencing reads from other platforms.

## Workflow Overview

![ViralFlow Workflow](_static/images/flowchart_en.png)

## Intrahost Analysis

ViralFlow presents its own algorithm to detect iSNV (intrahost Single Nucleotide Variant) regions.

For a given multiallelic site to have a lower frequency allele considered as iSNV, three conditions must be met, these conditions are schematized in the upper part of the figure below.

In this logic, the bottom part of the figure below represents 3 multi-allelic sites and exemplifies the logic for considering which of them should be considered as iSNVs.

![iSNV Detection Logic](_static/images/viralflow_snv.png)

## Documentation Contents

```{toctree}
:maxdepth: 2

installation
dependencies
quickstart
parameters
outputs
```

## Publications

- **Version 0.1**: Dezordi, F. Z., et al. (2022). ViralFlow: A Versatile Automated Workflow for SARS-CoV-2 Genome Assembly, Lineage Assignment, Mutations and Intrahost Variant Detection. *Viruses*, 14(2), 217. [https://www.mdpi.com/1999-4915/14/2/217](https://www.mdpi.com/1999-4915/14/2/217)

- **Version 1.0**: da Silva, A. F., et al. (2024). ViralFlow v1.0: characterization and annotation of viral genomes. *NAR Genomics and Bioinformatics*, 6(2), lqae056. [https://academic.oup.com/nargab/article/6/2/lqae056/7682253](https://academic.oup.com/nargab/article/6/2/lqae056/7682253)

## Quick Links

- [GitHub Repository](https://github.com/WallauBioinfo/ViralFlow)
- [Tool Versions](https://github.com/WallauBioinfo/ViralFlow/tree/main/versions)
- [Issues/Bugs](https://github.com/WallauBioinfo/ViralFlow/issues)
- **License:** Not specified
