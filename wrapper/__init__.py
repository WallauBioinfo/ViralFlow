from logging import root
import os
import shlex
import subprocess
import sys
import tempfile


def add_entries_to_DB(root_path, org_name, refseq_code, arch):
    """
    add entries provided to snpeff database
    """
    run_bash = ["bash", f"{root_path}/vfnext/containers/add_entries_SnpeffDB.sh",
                org_name, refseq_code, arch]
    subprocess.check_call(run_bash)

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
    containers_dir = os.path.join(root_path, "vfnext", "containers")
    for script in ("pull_containers.py", "build_containers.py"):
        subprocess.check_call([sys.executable, script, arch], cwd=containers_dir)
    

# input args file load
def parse_params(in_flpath):
    """
    load text file containing viralflow arguments
    """
    valid_args = [
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

def _update_pangolin(root_path, option):
    # Upstream setup.py still imports pkg_resources; constrain only build tools.
    with tempfile.TemporaryDirectory(prefix="viralflow-pangolin-") as temporary:
        with open(os.path.join(temporary, "build-constraints.txt"), "w") as constraints:
            constraints.write("setuptools<81\n")
        command = ["singularity", "exec", "--writable", "--bind", f"{temporary}:/tmp",
                   "./pangolin:4.4.sif", "env",
                   "PIP_BUILD_CONSTRAINT=/tmp/build-constraints.txt"]
        containers_dir = os.path.join(root_path, "vfnext", "containers")
        subprocess.check_call(command + ["pangolin", option], cwd=containers_dir)
        if option == "--update":
            # Conda package metadata may lag behind the upstream Python requirements.
            requirements = (
                "from importlib.metadata import requires; import subprocess, sys; "
                "subprocess.check_call([sys.executable, '-m', 'pip', 'install', "
                "*(requires('pangolin') or []), *(requires('snakemake') or [])])"
            )
            subprocess.check_call(command + ["python", "-c", requirements], cwd=containers_dir)
            subprocess.check_call(command + ["python", "-m", "pip", "check"], cwd=containers_dir)

def update_pangolin(root_path):
    _update_pangolin(root_path, "--update")

def update_pangolin_data(root_path):
    _update_pangolin(root_path, "--update-data")

def run_vfnext(root_path, params_fl):
    # get nextflow arguments
    args_str = parse_params(params_fl)
    nxtflw_ver = os.environ.get("NXF_VER", "22.04.0")
    configured_nextflow = os.environ.get("VIRALFLOW_NEXTFLOW")
    nextflow = configured_nextflow or "nextflow"
    if configured_nextflow and (
            not os.path.isabs(nextflow) or not os.path.isfile(nextflow)
            or not os.access(nextflow, os.X_OK)):
        raise FileNotFoundError(f"Nextflow executable not found at the configured absolute path: {nextflow}")
    run_nxtfl_cmd = (
        f"NXF_VER={shlex.quote(nxtflw_ver)} {shlex.quote(nextflow)} run "
        f"{shlex.quote(os.path.join(root_path, 'vfnext', 'main.nf'))} {args_str}"
    )
    print(run_nxtfl_cmd)
    result = subprocess.call(run_nxtfl_cmd, shell=True)
    if result:
        raise SystemExit(result)
