# Installation

## Ubuntu 26.04

This branch uses Micromamba 1.5.7, Nextflow 22.04.0 and Java 17.
ViralFlow invokes `singularity` for container build, download, maintenance and
execution, as in the upstream branch. The GUI installs Apptainer and provides
the `singularity` compatibility command. Pangolin and snpEff remain
**sandbox directories**, despite their `.sif` suffix; the other images are SIF.

The ViralFlow GUI manages the installation in separate steps: Micromamba,
Nextflow, Apptainer, ViralFlow and containers. Its default paths are
`~/ViralFlowGUI/bin`, `~/ViralFlowGUI/micromamba` and `~/ViralFlowGUI/ViralFlow`.
This branch does not provide `install.sh`.

### Official Quay containers

This branch incorporates upstream commit `77ca0ec` from
`WallauBioinfo/change_containers_repository`. Both repository manifests are
identical to upstream and reference `quay.io/wallaulabs/amd64-<image>:<tag>` or
`quay.io/wallaulabs/arm64-<image>:<tag>`. Downloads use `docker://` and the host's
`singularity` command to convert the OCI images into SIF files. Local filenames
remain `<image>:<tag>.sif`, preserving the existing Nextflow paths and reusing
containers already present. Failed downloads propagate a nonzero exit status.
No Sylabs remote or registry login is required for the public AMD64 images.

The AMD64 Pangolin and snpEff recipes are byte-for-byte copies of `77ca0ec`.
On 2026-11-08, both built successfully on Ubuntu 26.04.1 AMD64 with Apptainer
1.5.5. Pangolin 4.4 classified the public MN908947 reference; snpEff 5.0e
downloaded `NC_045512.2`, produced its database catalog and annotated a test VCF,
including HTML and gene reports. Pangolin's `pip check` reported mismatches
for pandas, scikit-learn and packaging; the official recipe is retained as-is.
No scientific algorithms, biological parameters or reference datasets were changed.

Anonymous requests for all 11 AMD64 manifests succeeded; the 11 ARM64 manifests
returned HTTP 401. ARM64 availability requires confirmation with the upstream
maintainers. This validation did not include a complete pipeline analysis or ARM64.

For manual installation, provide Micromamba and a working `singularity` command
on the host (Singularity or Apptainer's compatibility command), put them in
`PATH`, then create the environment from the appropriate YAML:

```bash
git clone -b feat/develop/ubuntu26.04-quay.io https://github.com/camilodotto/ViralFlow.git
cd ViralFlow
micromamba env create -n viralflow -f envs/amd64.yml -y
micromamba run -n viralflow python -m pip install -e .
NXF_VER=22.04.0 micromamba run -n viralflow nextflow -version
micromamba run -n viralflow viralflow --version
```

Java 17 and Nextflow 22.04.0 come from the YAML, matching the upstream
environment. The GUI also installs its own Nextflow launcher in `bin` and
selects that absolute path. The container runtime is provided by the host;
the YAML does not install a second runtime in the Conda environment.
For an existing environment, use
`micromamba update -n viralflow -f envs/amd64.yml --prune -y`.
The GUI pins its own Nextflow binary and keeps its cache under the installation
root. See [the migration record](../UBUNTU26.04-CHANGES.md) for changes and
validation. The remaining instructions describe the older installation paths.

## MacOS Installation

Due to the limitation of using singularity on MacOS, to run ViralFlow on this type of system, we suggest using a Linux virtualization software called [Lima](https://github.com/lima-vm/lima).

### Installing Lima

Install Lima using the Homebrew package manager:

```bash
brew install lima
```

### Installing the Ubuntu instance

Follow the step by step installation of an Ubuntu instance:

```bash
limactl start
```

### Start Ubuntu

When the virtual machine is installed, use the command below to start Ubuntu and follow the ViralFlow installation steps as usual:

```bash
lima
```

## Ubuntu Installation

To install ViralFlow, four steps are necessary:
1. Install system dependencies
2. Install Conda
3. Install ViralFlow  
4. Build the containers for analyses

This process is performed only once.

ViralFlow was developed and tested for the following operational systems:
- Ubuntu 20.04 LTS
- Ubuntu 22.04 LTS

### Installing System Dependencies

If you don't have the pip dependency installer, the git version control system, and the uidmap package, install them with:

```bash
sudo apt update -y && \
  sudo apt upgrade -y && \
  sudo apt install curl git python3-pip uidmap -y
```

### Installing and Configuring Micromamba

We recommend managing micromamba environments due to their parallelization when downloading and installing dependencies. If you use another environment manager (conda, miniconda, mamba), you can skip this step, however, conflicts during installation may occur.

If you don't have micromamba installed:

```bash
cd $HOME
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/1.5.7 | tar -xvj bin/micromamba
./bin/micromamba shell init -s bash -p ~/micromamba
source ~/.bashrc
micromamba activate
```

### Installing ViralFlow

If you already have the aforementioned dependencies and conda installed, you can install ViralFlow with 5 lines of code:

```bash
git clone https://github.com/WallauBioinfo/ViralFlow
cd ViralFlow/
micromamba env create -f envs/env.yml
micromamba activate viralflow
pip install -e .
```

### Building the Containers

All steps of ViralFlow are performed in controlled environments. ViralFlow has its own method to carry out all this construction of environments, running just one line of code.

For the building of containers, ViralFlow requires that the tool "unsquashfs" be available in the directory `/usr/local/bin/`. To ensure this, create a symbolic link:

```bash
sudo ln -s /usr/bin/unsquashfs /usr/local/bin/unsquashfs
```

After ensuring that "unsquashfs" is in the appropriate location, run the command to build the containers:

```bash
viralflow build_containers
```

```{note}
This process will download approximately 4.4GB of container images. Ensure you have sufficient disk space and a stable internet connection.
```
