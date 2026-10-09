#!/usr/bin/env python3
import os
import re
import subprocess
from importlib.metadata import version
from typing import Any, NamedTuple

import click
from . import (
    build_containers as _build_containers,
    update_pangolin as _update_pangolin,
    update_pangolin_data as _update_pangolin_data,
    add_entries_to_DB as _add_entries_to_DB,
    run_vfnext as _run_vfnext,
    concat_fastqs as _concat_fastqs,
    container_reference as _container_reference,
)

__version__ = version("ViralFlow")

# Get root paths
script_file = os.path.realpath(__file__)
VF_ROOT_PATH = "/".join(script_file.split("/")[0:-2]) + "/"


def _memory(_ctx, _param, value):
    """Accept what Nextflow's memory directive accepts: 8.GB, 8GB, 8 GB, 1.5 GB."""
    if value is None or re.fullmatch(
        r"\d+(\.\d+)?\s*\.?\s*[KMGTP]?B", value.strip(), flags=re.IGNORECASE
    ):
        return value
    raise click.BadParameter(f"'{value}' is not a memory size such as 8.GB")


def _container(_ctx, _param, value):
    return None if value is None else _container_reference(value)


class NanoporeOption(NamedTuple):
    flag: str
    param: str  # the nextflow.config name, as wrapper.NANOPORE_PARAMS lists it
    type: Any
    default: Any  # shown in --help only; nextflow.config applies it
    help: str
    callback: Any = None


# One option per NANOPORE parameter. Each defaults to None and is forwarded only
# when given, so nextflow.config stays the single source of these defaults; the
# `default` column is what --help shows, and tests/test_wrapper_nanopore.py fails
# if it stops matching nextflow.config.
NANOPORE_OPTIONS = (
    NanoporeOption(
        "--clair3-model",
        "clair3_model",
        str,
        "r941_prom_sup_g5014",
        "Clair3 model; must match the flowcell and basecaller (the default is "
        "for R9.4.1)",
    ),
    NanoporeOption(
        "--clair3-qual",
        "clair3_qual",
        click.IntRange(min=0),
        10,
        "Clair3 --qual: calls at or below it are labelled LowQual",
    ),
    NanoporeOption(
        "--clair3-chunk-size",
        "clair3_chunk_size",
        click.IntRange(min=1),
        10000,
        "Clair3 --chunk_size; affects runtime and memory, not results",
    ),
    NanoporeOption(
        "--af-threshold",
        "af_threshold",
        click.FloatRange(0, 1),
        0.51,
        "Keep a variant when FORMAT/AF is at least this",
    ),
    NanoporeOption(
        "--np-min-depth",
        "np_min_depth",
        click.IntRange(min=0),
        20,
        "Mask consensus positions with depth at or below this",
    ),
    NanoporeOption(
        "--base-container",
        "base_container",
        str,
        None,
        "Base image: a local SIF, or an image reference under -profile docker "
        "[default: set by nextflow.config or the profile]",
        _container,
    ),
    NanoporeOption(
        "--clair3-container",
        "clair3_container",
        str,
        None,
        "Clair3 image [default: the digest pinned in nextflow.config]",
        _container,
    ),
    NanoporeOption(
        "--porechop-cpus", "porechop_cpus", click.IntRange(min=1), 4, "Porechop CPUs"
    ),
    NanoporeOption(
        "--porechop-memory", "porechop_memory", str, "4.GB", "Porechop memory", _memory
    ),
    NanoporeOption(
        "--minimap-cpus", "minimap_cpus", click.IntRange(min=1), 4, "Minimap2 CPUs"
    ),
    NanoporeOption(
        "--minimap-memory", "minimap_memory", str, "4.GB", "Minimap2 memory", _memory
    ),
    NanoporeOption(
        "--clair3-cpus", "clair3_cpus", click.IntRange(min=1), 4, "Clair3 CPUs"
    ),
    NanoporeOption(
        "--clair3-memory", "clair3_memory", str, "4.GB", "Clair3 memory", _memory
    ),
)


def _nanopore_options(function):
    # Applied in reverse so --help lists them in table order.
    for option in reversed(NANOPORE_OPTIONS):
        shown = "" if option.default is None else f" [default: {option.default}]"
        function = click.option(
            option.flag,
            option.param,
            type=option.type,
            default=None,
            callback=option.callback,
            help=f"NANOPORE only. {option.help}{shown}",
        )(function)
    return function


def _call_helper(function, *args):
    try:
        return function(*args)
    except subprocess.CalledProcessError as error:
        raise click.ClickException(
            f"Command failed with exit status {error.returncode}: {error.cmd}"
        ) from error
    except (OSError, RuntimeError, ValueError) as error:
        raise click.ClickException(str(error)) from error


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
        click.echo(
            f"ViralFlow v{__version__} — a computational workflow for streamlining viral genomic surveillance."
        )
        click.echo("\nUse 'viralflow --help' to see available commands.")


# =============================================================================
# Container Management Commands
# =============================================================================


@cli.command("build-containers")
@click.option(
    "--arch",
    default="amd64",
    type=click.Choice(["amd64", "arm64"], case_sensitive=False),
    help="Architecture to build containers",
)
def build_containers(arch):
    """Build containers for vfnext."""
    _call_helper(_build_containers, VF_ROOT_PATH, arch)


@cli.command("update-pangolin")
def update_pangolin():
    """Update pangolin container to the latest pangolin version."""
    _call_helper(_update_pangolin, VF_ROOT_PATH)


@cli.command("update-pangolin-data")
def update_pangolin_data():
    """Update pangolin container with the latest pangolin version databases."""
    _call_helper(_update_pangolin_data, VF_ROOT_PATH)


# =============================================================================
# Database Management Commands
# =============================================================================


@cli.command("add-entry-to-snpeff")
@click.option("--org-name", required=True, type=str, help="Organism name")
@click.option(
    "--genome-code", required=True, type=str, help="Organism reference genome code"
)
@click.option(
    "--arch",
    default="amd64",
    type=click.Choice(["amd64", "arm64"], case_sensitive=False),
    help="System architecture",
)
def add_entry_to_snpeff(org_name, genome_code, arch):
    """Add a new entry to the SnpEff database."""
    _call_helper(_add_entries_to_DB, VF_ROOT_PATH, org_name, genome_code, arch)


# =============================================================================
# Pipeline Execution Commands
# =============================================================================


@cli.command("run")
@click.option(
    "--params-file",
    type=click.Path(exists=True),
    default=None,
    help="Path to a parameters file. File values override defaults.",
)
@click.option(
    "--profile",
    type=str,
    default=None,
    help="Nextflow profile (e.g. apptainer, fiocruz_default, fiocruz_pbs)",
)
@click.option(
    "--mode",
    type=click.Choice(["ILLUMINA", "NANOPORE"], case_sensitive=True),
    default=None,
    help="Sequencing technology used (defaults to ILLUMINA without --params-file)",
)
@click.option(
    "--virus",
    type=click.Choice(["sars-cov2", "custom"], case_sensitive=True),
    default="sars-cov2",
    show_default=True,
    help="Virus preset",
)
@click.option(
    "--in-dir",
    type=click.Path(),
    default=None,
    help="Deprecated input directory with FASTQ files (removed in v3)",
)
@click.option(
    "--samplesheet",
    type=click.Path(exists=True, dir_okay=False),
    default=None,
    help="CSV sample sheet with sample_id, fastq_1, and fastq_2 columns",
)
@click.option(
    "--out-dir",
    type=click.Path(),
    default="./output/",
    show_default=True,
    help="Output directory",
)
@click.option(
    "--primers-bed",
    type=click.Path(),
    default=None,
    help="BED file with primer positions",
)
@click.option(
    "--run-snpeff",
    is_flag=True,
    default=False,
    show_default=True,
    help="Run SnpEff annotation",
)
@click.option(
    "--write-mapped-reads",
    is_flag=True,
    default=False,
    show_default=True,
    help="Write mapped/unmapped reads",
)
@click.option(
    "--min-len", type=int, default=75, show_default=True, help="Minimum read length"
)
@click.option(
    "--depth",
    type=int,
    default=25,
    show_default=True,
    help="Minimum depth for consensus",
)
@click.option(
    "--min-dp-intrahost",
    type=int,
    default=100,
    show_default=True,
    help="Minimum depth for intrahost variants",
)
@click.option("--trim-len", type=int, default=0, show_default=True, help="Trim length")
@click.option(
    "--ref-genome-code",
    type=str,
    default=None,
    help="Reference genome code (e.g. NC_045512.2)",
)
@click.option(
    "--reference-gff",
    type=click.Path(),
    default=None,
    help="Reference GFF file (custom virus)",
)
@click.option(
    "--reference-genome",
    type=click.Path(),
    default=None,
    help="Reference genome fasta (custom virus)",
)
@click.option(
    "--nextflow-sim-calls",
    type=int,
    default=6,
    show_default=True,
    help="Max simultaneous Nextflow calls",
)
@click.option(
    "--fastp-threads",
    type=int,
    default=1,
    show_default=True,
    help="Number of threads for fastp",
)
@click.option(
    "--bwa-threads",
    type=int,
    default=1,
    show_default=True,
    help="Number of threads for bwa",
)
@click.option(
    "--mafft-threads",
    type=int,
    default=1,
    show_default=True,
    help="Number of threads for mafft",
)
@click.option(
    "--mapping-quality",
    type=int,
    default=30,
    show_default=True,
    help="Minimum mapping quality, in both modes (NANOPORE passes it to Clair3 "
    "as --min_mq)",
)
@click.option(
    "--base-quality",
    type=int,
    default=30,
    show_default=True,
    help="Minimum base quality",
)
@click.option(
    "--dedup/--no-dedup", default=False, show_default=True, help="Enable deduplication"
)
@click.option(
    "--ndedup",
    type=int,
    default=3,
    show_default=True,
    help="Number of allowed duplicates",
)
@_nanopore_options
def run(
    params_file,
    profile,
    mode,
    virus,
    in_dir,
    samplesheet,
    out_dir,
    primers_bed,
    run_snpeff,
    write_mapped_reads,
    min_len,
    depth,
    min_dp_intrahost,
    trim_len,
    ref_genome_code,
    reference_gff,
    reference_genome,
    nextflow_sim_calls,
    fastp_threads,
    bwa_threads,
    mafft_threads,
    mapping_quality,
    base_quality,
    dedup,
    ndedup,
    **nanopore_options,
):
    """Run the ViralFlow pipeline.

    All parameters have sensible defaults from nextflow.config.
    A parameter file is authoritative and cannot be combined with --mode.
    """
    # Only the options the user gave: an unset one leaves nextflow.config's
    # default in force rather than repeating it here.
    nanopore_given = {
        param: value for param, value in nanopore_options.items() if value is not None
    }
    if nanopore_given:
        flags = ", ".join(
            option.flag for option in NANOPORE_OPTIONS if option.param in nanopore_given
        )
        # Before the --mode check below: a params file sets its own mode.
        if params_file:
            raise click.UsageError(
                f"{flags} cannot be used with --params-file; set them in the "
                "parameter file"
            )
        if mode != "NANOPORE":
            raise click.UsageError(f"{flags} only apply to --mode NANOPORE")

    cli_to_nf = {
        "virus": virus,
        "inDir": in_dir,
        "samplesheet": samplesheet,
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
        **nanopore_given,
    }
    cli_params = {k: v for k, v in cli_to_nf.items() if v is not None}

    if params_file and mode is not None:
        raise click.UsageError(
            "--mode cannot be used with --params-file; set mode in the parameter file"
        )

    if in_dir and samplesheet:
        raise click.UsageError("--samplesheet and --in-dir cannot be used together")

    # Validate paths only when no params file is provided
    if not params_file:
        if in_dir and not os.path.exists(in_dir):
            raise click.BadParameter(
                f"Path '{in_dir}' does not exist.", param_hint="'--in-dir'"
            )

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
    help="Path to an input fastq files",
)
@click.option(
    "--prefix",
    type=str,
    default="barcode",
    help="Prefix of the fastq files, default: barcode",
)
@click.option(
    "--extension",
    type=str,
    default=".fastq.gz",
    help="Extension of the fastq files, default: .fastq.gz",
)
@click.option(
    "--min-len",
    type=int,
    default=200,
    help="Minimum length of the reads to write on output fastq files, default: 200",
)
# No maximum by default. It used to be 500, which kept ~400 bp amplicon runs
# but dropped every read of a ~1200 bp amplicon scheme and nearly every read
# of a whole-genome run, while reporting success.
@click.option(
    "--max-len",
    type=int,
    default=None,
    help="Maximum length of the reads to write on output fastq files, default: no maximum",
)
def concat_fastqs(path, prefix, extension, min_len, max_len):
    max_len_text = "no maximum length" if max_len is None else f"max length {max_len}"
    click.echo(
        f"Concat fastq files on path {path} with prefix {prefix} and extension {extension} with min length {min_len} and {max_len_text}"
    )
    _call_helper(_concat_fastqs, path, prefix, extension, min_len, max_len)


if __name__ == "__main__":
    cli()
