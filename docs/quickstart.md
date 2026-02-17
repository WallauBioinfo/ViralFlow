# Quick Start Guide

ViralFlow provides 2 usage modes: **sars-cov2** and **custom**. Regardless of the mode, the user must provide the absolute paths (the entire path to the directory or file) for each input file in the command line, otherwise the pipeline will be interrupted during execution.

```{important}
Always use **absolute paths** when specifying directories and files, e.g., `/home/user/data/` instead of `./data/`
```

## Customizing the snpEff Database

By default, the snpEff tool in ViralFlow is configured with the NC_045512.2 genome of the SARS-CoV-2 virus only. If you want to include the snpEff analysis for other viruses, you must update the snpEff database with the following command (example with Dengue):

```bash
viralflow -add_entry_to_snpeff --org_name Dengue --genome_code NC_001474.2
```

## SARS-CoV-2 Mode

In this model, the analysis is performed based on the reference genome NC_045512.2, and has, as additional analysis, the signature of strains with the Pangolin tool, and signing clades and mutations with the Nextclade tool.

For this, the user needs to build a file with the analysis parameters. An example can be [seen here](https://github.com/WallauBioinfo/ViralFlow/blob/main/test_files/sars-cov-2.params).

```bash
viralflow -run --params_file test_files/sars-cov-2.params
```

## Custom Mode

In this model, the analysis is performed based on the files for the virus that the user wants to analyze. In this mode, the user is responsible for providing each of the files necessary for the analysis. If the user wants to perform the snpEff analysis, he must pass the refseq code viral genome.

An example parameter file can be [seen here](https://github.com/WallauBioinfo/ViralFlow/blob/main/test_files/denv.params).

```bash
viralflow -run --params_file test_files/denv.params
```

## Pangolin Update

Periodically the pangolin tool updates the lineage database, as well as the usher classification phylogeny, the scorpion mutation constellations, and the pangoLearn trained model. To update the tool and/or its databases, just run ViralFlow with one of the commands:

```bash
# Update the tool and databases
viralflow -update_pangolin

# Update only the databases
viralflow -update_pangolin_data
```
