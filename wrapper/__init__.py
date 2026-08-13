import os
import subprocess
import shlex
import tempfile
import re
import sys
from pathlib import Path


def add_entries_to_DB(root_path, org_name, refseq_code, arch):
    """
    add entries provided to snpeff database
    """
    if any(ord(character) < 32 or ord(character) == 127 for character in org_name):
        raise ValueError("Organism name cannot contain control characters")
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", refseq_code):
        raise ValueError(
            "Genome code must start with a letter or number and contain only letters, numbers, dots, underscores, and hyphens"
        )
    containers_dir = Path(root_path) / "vfnext" / "containers"
    command = [
        "bash",
        str(containers_dir / "add_entries_SnpeffDB.sh"),
        org_name,
        refseq_code,
        arch,
    ]
    print(shlex.join(command))
    subprocess.run(command, cwd=containers_dir, check=True)


def parse_csv(csv_flpath):
    with open(csv_flpath, "r") as csv_fl:
        first_line = True
        entries_lst = []
        for line in csv_fl:
            # skip header
            if first_line:
                first_line = False
                continue
            ln_data = line.split(",")
            entry = [ln_data[0], ln_data[1].replace("\n", "")]
            entries_lst.append(entry)
    return entries_lst


def build_containers(root_path, arch: str):
    """
    run script to build container for vfnext
    """
    containers_dir = Path(root_path) / "vfnext" / "containers"
    subprocess.run(
        [sys.executable, "pull_containers.py", arch], cwd=containers_dir, check=True
    )
    subprocess.run(
        [sys.executable, "build_containers.py", arch], cwd=containers_dir, check=True
    )


# input args file load
def parse_params(in_flpath):
    """
    Load a legacy parameter file as a subprocess argument list.

    Path parameters consume the entire remainder of their line so unquoted paths
    containing spaces remain a single argument. All other parameters are scalar.
    """
    valid_args = {
        "mode",
        "virus",
        "primersBED",
        "outDir",
        "inDir",
        "samplesheet",
        "runSnpEff",
        "writeMappedReads",
        "minLen",
        "depth",
        "minDpIntrahost",
        "trimLen",
        "runSnpEff",
        "refGenomeCode",
        "referenceGFF",
        "referenceGenome",
        "nextflowSimCalls",
        "fastp_threads",
        "bwa_threads",
        "mafft_threads",
        "nxtclade_jobs",
        "mapping_quality",
        "base_quality",
        "dedup",
        "ndedup",
    }
    path_params = {
        "inDir",
        "samplesheet",
        "outDir",
        "referenceGFF",
        "referenceGenome",
        "primersBED",
    }
    parsed = {}

    with open(in_flpath, "r") as in_file:
        for line_number, raw_line in enumerate(in_file, start=1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue

            fields = line.split(maxsplit=1)
            key = fields[0]
            if key not in valid_args:
                raise ValueError(f"Line {line_number}: {key} is not a valid argument")

            value = fields[1].strip() if len(fields) == 2 else ""
            if not value or value == "null":
                continue
            if key not in path_params and len(value.split()) > 1:
                raise ValueError(f"Line {line_number}: {key} accepts a single value")

            parsed[key] = os.path.abspath(value) if key in path_params else value

    args = []
    for key, value in parsed.items():
        option = f"--{key}"
        args.extend([option, value])
    args.append("-resume")
    return args


def update_pangolin(root_path):
    containers_dir = Path(root_path) / "vfnext" / "containers"
    subprocess.run(
        [
            "singularity",
            "exec",
            "--writable",
            "./pangolin:4.4.sif",
            "pangolin",
            "--update",
        ],
        cwd=containers_dir,
        check=True,
    )


def update_pangolin_data(root_path):
    containers_dir = Path(root_path) / "vfnext" / "containers"
    subprocess.run(
        [
            "singularity",
            "exec",
            "--writable",
            "./pangolin:4.4.sif",
            "pangolin",
            "--update-data",
        ],
        cwd=containers_dir,
        check=True,
    )


def run_vfnext(root_path, params_fl, mode, cli_params=None, profile=None):
    """
    Run the vfnext pipeline.

    If params_fl is provided, file parameters are authoritative, including mode.
    If no params_fl, CLI parameters are used and mode defaults to ILLUMINA.
    """
    path_params = [
        "inDir",
        "samplesheet",
        "outDir",
        "referenceGFF",
        "referenceGenome",
        "primersBED",
    ]

    if params_fl:
        if mode is not None:
            raise ValueError("mode cannot be provided with a parameter file")
        # Params file takes full precedence — do not append CLI defaults.
        args = parse_params(params_fl)
    else:
        # No file provided — use CLI params
        if cli_params:
            for k in path_params:
                if k in cli_params:
                    cli_params[k] = os.path.abspath(str(cli_params[k]))
            args = []
            for key, value in cli_params.items():
                option = f"--{key}"
                args.extend([option, str(value)])
        else:
            raise ValueError(
                "No parameters provided. Use --params-file or individual CLI options."
            )

    if "-resume" not in args:
        args.append("-resume")

    nxtflw_ver = os.environ.get("NXF_VER", "26.04.6")
    resolved_mode = None if params_fl else (mode or "ILLUMINA")
    run_env = os.environ.copy()
    run_env["NXF_VER"] = nxtflw_ver
    command = ["nextflow", "run", f"{root_path}/vfnext/main.nf"]
    command.extend(args)
    if resolved_mode:
        command.extend(["--mode", resolved_mode])
    if profile:
        command.extend(["-profile", profile])
    print(f"NXF_VER={shlex.quote(nxtflw_ver)} {shlex.join(command)}")
    subprocess.run(command, env=run_env, check=True)


def concat_fastqs(path, prefix, extension, min_len, max_len):
    """
    Concatenate and filter fastq files from barcode directories.

    For each directory matching {prefix}* (e.g. barcode01, barcode02, ...),
    concatenates all fastq files and filters reads by min/max length using seqkit.
    """
    read_dir = Path(path).resolve()
    output_dir = read_dir / "filtered"
    output_dir.mkdir(exist_ok=True)

    def matching_directories(parent):
        return sorted(
            entry
            for entry in parent.iterdir()
            if entry.is_dir() and entry.name.startswith(prefix)
        )

    barcode_dirs = matching_directories(read_dir)
    if not barcode_dirs:
        barcode_dirs = sorted(
            barcode
            for parent in read_dir.iterdir()
            if parent.is_dir()
            for barcode in matching_directories(parent)
        )
    if not barcode_dirs:
        raise FileNotFoundError(
            f"No files matching '{prefix}*/*{extension}' found in '{read_dir}' or its subdirectories."
        )

    failures = []
    for barcode_path in barcode_dirs:
        barcode_id = barcode_path.name
        fastq_files = sorted(
            entry
            for entry in barcode_path.iterdir()
            if entry.is_file() and entry.name.endswith(extension)
        )
        if not fastq_files:
            print(f"Skipping {barcode_id}: no {extension} files found")
            continue

        output_file = output_dir / f"{barcode_id}.concat.fastq.gz"
        temporary_path = None
        print(f"Processing {barcode_id}...")
        processes = []
        try:
            with tempfile.NamedTemporaryFile(
                dir=output_dir, prefix=f".{barcode_id}.", suffix=".tmp", delete=False
            ) as temporary:
                temporary_path = Path(temporary.name)
                reader = (
                    ["gzip", "-cd", "--", *map(str, fastq_files)]
                    if extension.endswith(".gz")
                    else ["cat", "--", *map(str, fastq_files)]
                )
                source = subprocess.Popen(reader, stdout=subprocess.PIPE)
                processes.append(source)
                seqkit = subprocess.Popen(
                    [
                        "seqkit",
                        "seq",
                        "-w",
                        "0",
                        "-g",
                        "--min-len",
                        str(min_len),
                        "--max-len",
                        str(max_len),
                    ],
                    stdin=source.stdout,
                    stdout=subprocess.PIPE,
                )
                processes.append(seqkit)
                source.stdout.close()
                compressor = subprocess.Popen(
                    ["gzip", "-c"], stdin=seqkit.stdout, stdout=temporary
                )
                processes.append(compressor)
                seqkit.stdout.close()
                compressor_status = compressor.wait()
                seqkit_status = seqkit.wait()
                source_status = source.wait()
            statuses = {
                "reader": source_status,
                "seqkit": seqkit_status,
                "gzip": compressor_status,
            }
            failed_stages = [name for name, status in statuses.items() if status != 0]
            if failed_stages:
                raise RuntimeError(
                    f"failed pipeline stages: {', '.join(failed_stages)}"
                )
            os.replace(temporary_path, output_file)
            print(f"   The reads were written in {output_file}")
        except (OSError, RuntimeError) as error:
            for process in processes:
                if process.poll() is None:
                    process.terminate()
            for process in processes:
                if process.poll() is None:
                    process.wait()
            if temporary_path is not None:
                temporary_path.unlink(missing_ok=True)
            failures.append(f"{barcode_id}: {error}")

    if failures:
        raise RuntimeError("FASTQ concatenation failed:\n - " + "\n - ".join(failures))
