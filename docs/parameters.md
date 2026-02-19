# Parameters

This page describes all command-line arguments and parameters available in ViralFlow.

## Parameter File

ViralFlow requires a parameter file that contains all configuration options. Example parameter files can be found in the [test_files directory](https://github.com/WallauBioinfo/ViralFlow/tree/main/test_files).

## Parameter Reference

| Argument | Default Value | Description |
|----------|---------------|-------------|
| `virus` | sars-cov2 | Analysis type (sars-cov2 or custom) |
| `primersBED` | null | Absolute path to bed file with primers information used in genomic amplification (optional) |
| `outDir` | launchDir/output/ | Absolute path to the directory where the results will be stored |
| `inDir` | launchDir/input/ | Absolute path to the directory with the input data (directory with the FASTQ files) |
| `runSnpEff` | true | Needed to run the snpEff tool (true or false) |
| `writeMappedReads` | true | Needed to generate the FASTQ files containing the sequencing reads that mapped to the reference genome |
| `minLen` | 75 | Minimum size the reads must have. Reads below this threshold will be eliminated by FastP |
| `depth` | 5 | Minimum coverage depth to call consensus bases. Positions with lower coverage depth will not be called and a "-" will be added to the respective consensus genomic position |
| `mapping_quality` | 30 | Mapping quality threshold used to variant calling |
| `base_quality` | 30 | Base quality threshold used to variant calling |
| `minDpIntrahost` | 100 | Minimum coverage depth per genomic site to be considered in the intrahost analysis |
| `trimLen` | 0 | Number of bases that should be trimmed at both ends of the reads |
| `refGenomeCode` | null | Code of the genome to be used in the custom analysis |
| `referenceGFF` | null | GFF genome file to be used in custom analysis |
| `referenceGenome` | null | Fasta genome file to be used in custom analysis |
| `nextflowSimCalls` | null | Number of simultaneous calls that nextflow can perform |
| `fastp_threads` | 1 | Number of threads to be used in the fastp read filtering step |
| `bwa_threads` | 1 | Number of threads to be used in the bwa mapping step |
| `mafft_threads` | 1 | Number of threads to be used in the mafft alignment step |
| `dedup` | false | This argument enable dedup mode of fastp. To activate it change value to true on params test file |
| `ndedup` | 3 | When dedup mode is active you can use accuracy levels (1 - 6). You can change this value, but we recommend the standard. How much higher, more RAM and time are consumed. To activate it change value from 1 to 6 on params test file |
