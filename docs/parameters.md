# Parameters

This page describes all command-line arguments and parameters available in ViralFlow.

## Parameter File

ViralFlow requires a parameter file that contains all configuration options. Example parameter files can be found in the [test_files directory](https://github.com/WallauBioinfo/ViralFlow/tree/main/test_files).

## Parameter Reference

`--outDir` is the single output location for published results and `RUN_METADATA`. Do not combine it with Nextflow's `-output-dir`; ViralFlow rejects conflicting values because execution observer paths are resolved from `--outDir` during configuration.

| Argument | Default Value | Description |
|----------|---------------|-------------|
| `mode` | ILLUMINA | Sequencing technology of the input data (ILLUMINA or NANOPORE) |
| `virus` | sars-cov2 | Analysis type (sars-cov2 or custom) |
| `primersBED` | null | Absolute path to bed file with primers information used in genomic amplification (optional). Supplying it enables primer clipping in both ILLUMINA and NANOPORE modes; without it, no clipping is performed |
| `outDir` | launchDir/output/ | ViralFlow output directory where results and run metadata will be stored |
| `samplesheet` | null | CSV file with `sample_id`, `fastq_1` and `fastq_2` columns. Repeat a sample ID to provide ordered chunks or lanes. Cannot be combined with `inDir` |
| `inDir` | null | Deprecated, and removed in v3; use `samplesheet`. Absolute path to the directory with the input data (directory with the FASTQ files). When unset, ViralFlow falls back to `launchDir/input/` |
| `runSnpEff` | true | Needed to run the snpEff tool (true or false) |
| `writeMappedReads` | true | Needed to generate the FASTQ files containing the sequencing reads that mapped to the reference genome |
| `minLen` | 75 | Minimum size the reads must have. Reads below this threshold will be eliminated by FastP |
| `depth` | 25 | Minimum coverage depth to call consensus bases. Positions with lower coverage depth will not be called and a "-" will be added to the respective consensus genomic position |
| `mapping_quality` | 30 | Mapping quality threshold used to variant calling |
| `base_quality` | 30 | Base quality threshold used to variant calling |
| `minDpIntrahost` | 100 | Minimum coverage depth per genomic site to be considered in the intrahost analysis |
| `trimLen` | 0 | Number of bases trimmed from both ends of the reads; 0 disables it. ILLUMINA applies it in fastp before alignment; NANOPORE masks the bases in the aligned BAM with bamUtil |
| `refGenomeCode` | null | Code of the genome to be used in the custom analysis |
| `referenceGFF` | null | GFF genome file to be used in custom analysis |
| `referenceGenome` | null | Fasta genome file to be used in custom analysis |
| `nextflowSimCalls` | 6 | Number of simultaneous calls that nextflow can perform |
| `fastp_threads` | 1 | Number of threads to be used in the fastp read filtering step |
| `bwa_threads` | 1 | Number of threads to be used in the bwa mapping step |
| `mafft_threads` | 1 | Number of threads to be used in the mafft alignment step |
| `dedup` | false | This argument enable dedup mode of fastp. To activate it change value to true on params test file |
| `ndedup` | 3 | When dedup mode is active you can use accuracy levels (1 - 6). You can change this value, but we recommend the standard. How much higher, more RAM and time are consumed. To activate it change value from 1 to 6 on params test file |

## NANOPORE Parameters

These apply only when `mode` is `NANOPORE` and are ignored otherwise. `virus`,
`refGenomeCode`, `referenceGFF`, `runSnpEff` and the fastp/bwa/mafft settings
above belong to ILLUMINA mode; `base_quality` is likewise ILLUMINA only, while
`mapping_quality` is used by both.

| Argument | Default Value | Description |
|----------|---------------|-------------|
| `base_container` | projectDir/containers/baseContainer.sif | Container providing Porechop_ABI, Minimap2, Samtools, BCFtools and bamUtil. A local `.sif` path under Singularity/Apptainer; the `docker` profile overrides it with an image reference |
| `clair3_container` | docker://hkubal/clair3@sha256:1430f7b5… | Clair3 image, pinned by digest so a rebuild cannot change the variant caller. Corresponds to release v1.2.0 |
| `clair3_model` | r941_prom_sup_g5014 | Basecalling model Clair3 uses, passed as `--model_path`. Must name a directory present under `/opt/models` inside the Clair3 image, and should match the basecaller and chemistry that produced the reads |
| `clair3_qual` | 10 | Minimum variant quality for a call to reach the consensus, passed to Clair3 as `--qual`. Clair3 still reports calls below it, marked `FILTER=LowQual` in `merge_output.vcf.gz`; the filtered VCF and the consensus keep only `PASS` calls |
| `clair3_chunk_size` | 10000 | Size in bases of the chunks Clair3 splits the reference into for parallel calling, passed as `--chunk_size`. Affects runtime and memory, not results |
| `mapping_quality` | 30 | Minimum mapping quality for a read to be used in variant calling, passed to Clair3 as `--min_mq` |
| `af_threshold` | 0.51 | Allele-frequency cutoff applied to Clair3's output: BCFtools keeps a variant when `FORMAT/AF >= af_threshold` and Clair3 marked it `FILTER=PASS` (see `clair3_qual`). No variant-depth condition is applied alongside it. The default above 0.5 keeps the majority allele at each site |
| `np_min_depth` | 20 | Consensus masking threshold. Coverage comes from `samtools depth -J -aa`, and every position whose depth is **less than or equal to** this value is written as `N`. At the default, a position needs at least 21 reads to be called |
| `porechop_cpus` | 4 | CPUs for the Porechop_ABI adapter-removal step |
| `porechop_memory` | 4.GB | Memory for the Porechop_ABI step |
| `minimap_cpus` | 4 | CPUs for the Minimap2 alignment step |
| `minimap_memory` | 4.GB | Memory for the Minimap2 step |
| `clair3_cpus` | 4 | CPUs for Clair3, passed as `--threads` |
| `clair3_memory` | 4.GB | Memory for Clair3 |

`af_threshold` and `np_min_depth` are applied independently, so a low-depth
variant can survive the allele-frequency filter while the same position is
masked in the consensus. Each sample directory contains a
`<sample>.nanopore_summary.tsv` reporting the configured thresholds, variant counts
and masked-base totals so this can be inspected per run. Note that `masked_bases`
and `consensus_n_bases` count different things — see `NANOPORE.md` for why they
differ on a reference containing ambiguity codes.
