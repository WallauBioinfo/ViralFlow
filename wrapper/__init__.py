import os
import glob


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
    load text file containing viralflow arguments
    """
    valid_args = [
        "mode",
        "virus",
        "primersBED",
        "outDir",
        "inDir",
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
    ]
    path_params = ["inDir", "outDir", "referenceGFF", "referenceGenome", "primersBED"]
    in_file = open(in_flpath, "r")
    dct = {}
    for l in in_file:
        # skip lines
        if (l in ["", " ", "\n"]) or l.startswith("#"):
            continue

        # get line data
        l_dt = l.replace("\n", "").split(" ")
        
        # get content
        key = l_dt[0]
        if (key not in valid_args):
            raise Exception(f"ERROR: {key} not a valid argument")
        # fill dict
        if key in valid_args:
        
            vls_1 = l_dt[1 : len(l_dt)]
            vls = []
        
            for v in vls_1:
                if v in [""]:
                    continue
                vls.append(v)
            # if single value
            if len(vls) == 1:
                # skip null values
                if vls[0] == "null":
                    continue
                # be sure paths are absolute
                if key in path_params:
                    dct[key] = os.path.abspath(vls[0])
                    continue
                dct[key] = vls[0]
            # if a list of values
            if len(vls) > 1:
                dct[key] = vls
            continue
    # get arguments for nextflow
    
    args_str = ""
    for key in dct:
        args_str += f"--{key} {dct[key]} "
    args_str += "-resume"
    return args_str

def update_pangolin(root_path):
    cd_to_dir= f"cd {root_path}/vfnext/containers/" 
    run_update = "singularity exec --writable ./pangolin:4.3.sif pangolin --update"
    os.system(cd_to_dir+';'+run_update)

def update_pangolin_data(root_path):
    cd_to_dir= f"cd {root_path}/vfnext/containers/" 
    run_update_data = "singularity exec --writable ./pangolin:4.3.sif pangolin --update-data"
    os.system(cd_to_dir+';'+run_update_data)

def run_vfnext(root_path, params_fl, mode, cli_params=None, profile=None):
    """
    Run the vfnext pipeline.

    If params_fl is provided, file parameters are the rule (CLI defaults are ignored).
    If no params_fl, CLI parameters are used.
    """
    path_params = ["inDir", "outDir", "referenceGFF", "referenceGenome", "primersBED"]

    if params_fl:
        # Params file takes full precedence — ignore CLI defaults
        args_str = parse_params(params_fl)
    else:
        # No file provided — use CLI params
        if cli_params:
            for k in path_params:
                if k in cli_params:
                    cli_params[k] = os.path.abspath(str(cli_params[k]))
            args_str = " ".join(f"--{k} {v}" for k, v in cli_params.items())
        else:
            raise ValueError("No parameters provided. Use --params-file or individual CLI options.")

    if "-resume" not in args_str:
        args_str += " -resume"

    nxtflw_ver="22.04.0"
    profile_str = f" -profile {profile}" if profile else ""
    run_nxtfl_cmd = f"NXF_VER={nxtflw_ver} nextflow run {root_path}/vfnext/main.nf {args_str} --mode {mode}{profile_str}"
    print(run_nxtfl_cmd)
    os.system(run_nxtfl_cmd)


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
