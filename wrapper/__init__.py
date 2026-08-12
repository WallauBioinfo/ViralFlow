import os
import glob
import subprocess
import shlex


def add_entries_to_DB(root_path, org_name, refseq_code, arch):
    """
    add entries provided to snpeff database
    """
    run_bash = f"bash {root_path}/vfnext/containers/add_entries_SnpeffDB.sh"
    print(f"{run_bash} {org_name} {refseq_code} {arch}")
    os.system(f"{run_bash} {org_name} {refseq_code} {arch}")

def parse_csv(csv_flpath):
    with open(csv_flpath, "r") as csv_fl:
        first_line = True
        entries_lst = []
        for line in csv_fl:
            # skip header
            if first_line == True:
                first_line = False
                continue
            ln_data = line.split(",")
            entry = [ln_data[0], ln_data[1].replace("\n","")]
            entries_lst.append(entry)
    return entries_lst

def build_containers(root_path, arch: str):
    """
    run script to build container for vfnext
    """
    # build containers
    cd_to_dir= f"cd {root_path}/vfnext/containers/" 
    build_sandbox = f"python ./build_containers.py {arch}"
    pull_containers = f"python ./pull_containers.py {arch}"
    os.system(cd_to_dir+';'+pull_containers) 
    print(cd_to_dir+';'+build_sandbox)
    os.system(cd_to_dir+';'+build_sandbox)
    

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
        "ndedup"
    }
    path_params = {"inDir", "samplesheet", "outDir", "referenceGFF", "referenceGenome", "primersBED"}
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
        option = "-output-dir" if key == "outDir" else f"--{key}"
        args.extend([option, value])
    args.append("-resume")
    return args

def update_pangolin(root_path):
    cd_to_dir= f"cd {root_path}/vfnext/containers/" 
    run_update = "singularity exec --writable ./pangolin:4.4.sif pangolin --update"
    os.system(cd_to_dir+';'+run_update)

def update_pangolin_data(root_path):
    cd_to_dir= f"cd {root_path}/vfnext/containers/" 
    run_update_data = "singularity exec --writable ./pangolin:4.4.sif pangolin --update-data"
    os.system(cd_to_dir+';'+run_update_data)

def run_vfnext(root_path, params_fl, mode, cli_params=None, profile=None):
    """
    Run the vfnext pipeline.

    If params_fl is provided, file parameters are authoritative, including mode.
    If no params_fl, CLI parameters are used and mode defaults to ILLUMINA.
    """
    path_params = ["inDir", "samplesheet", "outDir", "referenceGFF", "referenceGenome", "primersBED"]

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
                option = "-output-dir" if key == "outDir" else f"--{key}"
                args.extend([option, str(value)])
        else:
            raise ValueError("No parameters provided. Use --params-file or individual CLI options.")

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
    read_dir = os.path.abspath(path)
    output_dir = os.path.join(read_dir, "filtered")
    os.makedirs(output_dir, exist_ok=True)

    # Search for barcode directories in read_dir
    barcode_dirs = sorted(glob.glob(os.path.join(read_dir, f"{prefix}*")))

    # If not found, search in subdirectories (read_dir/*/)
    if not barcode_dirs:
        barcode_dirs = sorted(glob.glob(os.path.join(read_dir, "*", f"{prefix}*")))

    # If still not found, raise an error
    if not barcode_dirs:
        raise FileNotFoundError(
            f"No files matching '{prefix}*/*{extension}' found in '{read_dir}' or its subdirectories."
        )

    for barcode_path in barcode_dirs:
        if not os.path.isdir(barcode_path):
            continue

        barcode_id = os.path.basename(barcode_path)

        fastq_files = glob.glob(os.path.join(barcode_path, f"*{extension}"))
        if not fastq_files:
            print(f"Skipping {barcode_id}: no {extension} files found")
            continue

        output_file = os.path.join(output_dir, f"{barcode_id}.concat.fastq.gz")
        fastq_pattern = os.path.join(barcode_path, f"*{extension}")

        cmd = (
            f"zcat {fastq_pattern} | "
            f"seqkit seq -w 0 -g --min-len {min_len} --max-len {max_len} | "
            f"gzip > {output_file}"
        )

        print(f"Processing {barcode_id}...")
        os.system(cmd)
        print(f"   The reads were written in {output_file}")
