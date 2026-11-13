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
@click.option(
    "--params-file",
    required=True,
    type=click.Path(exists=True),
    help="Path to an input parameters file"
)
@click.option(
    "--dedup",
    is_flag=True,
    default=False,
    help="Enable dedup mode on fastp"
)
@click.option(
    "--ndedup",
    type=int,
    default=3,
    help="Standard accuracy level for dedup mode on fastp"
)
def run(params_file, dedup, ndedup):
    click.echo(f"ViralFlow v{__version__}")
    _run_vfnext(VF_ROOT_PATH, params_file)


if __name__ == "__main__":
    cli()
