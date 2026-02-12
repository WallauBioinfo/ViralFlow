#!/usr/bin/env python3
import os
from importlib.metadata import version

import click
from . import (
    build_containers as _build_containers,
    update_pangolin as _update_pangolin,
    update_pangolin_data as _update_pangolin_data,
    add_entries_to_DB as _add_entries_to_DB,
    run_vfnext as _run_vfnext,
    concat_fastqs as _concat_fastqs
)

__version__ = version("ViralFlow")

# Get root paths
script_file = os.path.realpath(__file__)
VF_ROOT_PATH = "/".join(script_file.split("/")[0:-2]) + "/"


@click.group(invoke_without_command=True)
@click.version_option(version=__version__, prog_name="ViralFlow")
@click.pass_context
def cli(ctx):
    """
    ViralFlow — a computational workflow for streamlining viral genomic surveillance.

    This is just a wrapper for a nextflow pipeline.
    For users familiar with nextflow, the directory vfnext holds the pipeline
    and usage directly via nextflow is strongly recommended.

    This wrapper was designed to make vfnext accessible for non-technical users.
    """
    if ctx.invoked_subcommand is None:
        click.echo(f"ViralFlow v{__version__} — a computational workflow for streamlining viral genomic surveillance.")
        click.echo("\nUse 'viralflow --help' to see available commands.")


# =============================================================================
# Container Management Commands
# =============================================================================

@cli.command("build-containers")
@click.option(
    "--arch",
    default="amd64",
    type=click.Choice(["amd64", "arm64"], case_sensitive=False),
    help="Architecture to build containers"
)
def build_containers(arch):
    """Build containers for vfnext."""
    _build_containers(VF_ROOT_PATH, arch)


@cli.command("update-pangolin")
def update_pangolin():
    """Update pangolin container to the latest pangolin version."""
    _update_pangolin(VF_ROOT_PATH)


@cli.command("update-pangolin-data")
def update_pangolin_data():
    """Update pangolin container with the latest pangolin version databases."""
    _update_pangolin_data(VF_ROOT_PATH)


# =============================================================================
# Database Management Commands
# =============================================================================

@cli.command("add-entry-to-snpeff")
@click.option(
    "--org-name",
    required=True,
    type=str,
    help="Organism name"
)
@click.option(
    "--genome-code",
    required=True,
    type=str,
    help="Organism reference genome code"
)
@click.option(
    "--arch",
    default="amd64",
    type=click.Choice(["amd64", "arm64"], case_sensitive=False),
    help="System architecture"
)
def add_entry_to_snpeff(org_name, genome_code, arch):
    """Add a new entry to the SnpEff database."""
    _add_entries_to_DB(VF_ROOT_PATH, org_name, genome_code, arch)


# =============================================================================
# Pipeline Execution Commands
# =============================================================================

@cli.command("run")
@click.option("--params-file", type=click.Path(exists=True), default=None,
              help="Path to a parameters file. File values override defaults.")
@click.option("--profile", type=str, default=None,
              help="Nextflow profile (e.g. apptainer, fiocruz_default, fiocruz_pbs)")
@click.option("--mode", type=click.Choice(["ILLUMINA", "NANOPORE"], case_sensitive=True),
              default="ILLUMINA", show_default=True, help="Sequencing technology used")
@click.option("--virus", type=click.Choice(["sars-cov2", "custom"], case_sensitive=True),
              default="sars-cov2", show_default=True, help="Virus preset")
@click.option("--in-dir", type=click.Path(), default="./input/",
              show_default=True, help="Input directory with fastq files")
@click.option("--out-dir", type=click.Path(), default="./output/",
              show_default=True, help="Output directory")
@click.option("--primers-bed", type=click.Path(), default=None,
              help="BED file with primer positions")
@click.option("--run-snpeff", is_flag=True, default=False,
              show_default=True, help="Run SnpEff annotation")
@click.option("--write-mapped-reads", is_flag=True, default=False,
              show_default=True, help="Write mapped/unmapped reads")
@click.option("--min-len", type=int, default=75,
              show_default=True, help="Minimum read length")
@click.option("--depth", type=int, default=25,
              show_default=True, help="Minimum depth for consensus")
@click.option("--min-dp-intrahost", type=int, default=100,
              show_default=True, help="Minimum depth for intrahost variants")
@click.option("--trim-len", type=int, default=0,
              show_default=True, help="Trim length")
@click.option("--ref-genome-code", type=str, default=None,
              help="Reference genome code (e.g. NC_045512.2)")
@click.option("--reference-gff", type=click.Path(), default=None,
              help="Reference GFF file (custom virus)")
@click.option("--reference-genome", type=click.Path(), default=None,
              help="Reference genome fasta (custom virus)")
@click.option("--nextflow-sim-calls", type=int, default=6,
              show_default=True, help="Max simultaneous Nextflow calls")
@click.option("--fastp-threads", type=int, default=1,
              show_default=True, help="Number of threads for fastp")
@click.option("--bwa-threads", type=int, default=1,
              show_default=True, help="Number of threads for bwa")
@click.option("--mafft-threads", type=int, default=1,
              show_default=True, help="Number of threads for mafft")
@click.option("--mapping-quality", type=int, default=30,
              show_default=True, help="Minimum mapping quality")
@click.option("--base-quality", type=int, default=30,
              show_default=True, help="Minimum base quality")
@click.option("--dedup/--no-dedup", default=False,
              show_default=True, help="Enable deduplication")
@click.option("--ndedup", type=int, default=3,
              show_default=True, help="Number of allowed duplicates")
def run(params_file, profile, mode, virus, in_dir, out_dir, primers_bed, run_snpeff,
        write_mapped_reads, min_len, depth, min_dp_intrahost, trim_len,
        ref_genome_code, reference_gff, reference_genome, nextflow_sim_calls,
        fastp_threads, bwa_threads, mafft_threads, mapping_quality,
        base_quality, dedup, ndedup):
    """Run the ViralFlow pipeline.

    All parameters have sensible defaults from nextflow.config.
    Optionally pass --params-file to load from a file (CLI options still override).
    """
    cli_to_nf = {
        "virus": virus,
        "inDir": in_dir,
        "outDir": out_dir,
        "primersBED": primers_bed,
        "runSnpEff": run_snpeff,
        "writeMappedReads": write_mapped_reads,
        "minLen": min_len,
        "depth": depth,
        "minDpIntrahost": min_dp_intrahost,
        "trimLen": trim_len,
        "refGenomeCode": ref_genome_code,
        "referenceGFF": reference_gff,
        "referenceGenome": reference_genome,
        "nextflowSimCalls": nextflow_sim_calls,
        "fastp_threads": fastp_threads,
        "bwa_threads": bwa_threads,
        "mafft_threads": mafft_threads,
        "mapping_quality": mapping_quality,
        "base_quality": base_quality,
        "dedup": dedup,
        "ndedup": ndedup,
    }
    cli_params = {k: v for k, v in cli_to_nf.items() if v is not None}

    # Validate paths only when no params file is provided
    if not params_file:
        if not os.path.exists(in_dir):
            raise click.BadParameter(f"Path '{in_dir}' does not exist.", param_hint=f"'--in-dir'")

    for k, v in cli_params.items():
        if isinstance(v, bool):
            cli_params[k] = str(v).lower()

    click.echo(f"ViralFlow v{__version__}")
    _run_vfnext(VF_ROOT_PATH, params_file, mode, cli_params, profile)

# =============================================================================
#  Concat fastq commands
# =============================================================================

@cli.command("concat-fastq")
@click.option(
    "--path",
    required=True,
    type=click.Path(exists=True),
    help="Path to an input fastq files"
)
@click.option(
    "--prefix",
    type=str,
    default="barcode",
    help="Prefix of the fastq files, default: barcode"
)
@click.option(
    "--extension",
    type=str,
    default=".fastq.gz",
    help="Extension of the fastq files, default: .fastq.gz"
)
@click.option(
    "--min-len",
    type=int,
    default=200,
    help="Minimum length of the reads to write on output fastq files, default: 200"
)
@click.option(
    "--max-len",
    type=int,
    default=500,
    help="Maximum length of the reads to write on output fastq files, default: 500"
)
def concat_fastqs(path, prefix, extension, min_len, max_len):
    click.echo(f"Concat fastq files on path {path} with prefix {prefix} and extension {extension} with min length {min_len} and max length {max_len}")
    _concat_fastqs(path, prefix, extension, min_len, max_len)

if __name__ == "__main__":
    cli()
