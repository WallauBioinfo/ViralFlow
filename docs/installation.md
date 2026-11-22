# Installation

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

ViralFlow provides one environment file per architecture: `envs/amd64.yml` (x86_64) and `envs/arm64.yml` (ARM64/aarch64). Check your architecture with `uname -m` (`x86_64` → AMD64, `aarch64` → ARM64).

#### AMD64

```bash
git clone https://github.com/WallauBioinfo/ViralFlow
cd ViralFlow/
micromamba env create -f envs/amd64.yml
micromamba activate viralflow
pip install -e .
```

#### ARM64

On ARM64, Singularity is not available through conda and must be compiled from source before installing ViralFlow (micromamba must also be downloaded from `linux-aarch64` instead of `linux-64` in the step above):

```bash
sudo apt update
sudo apt install -y build-essential git wget pkg-config \
    libseccomp-dev squashfs-tools cryptsetup \
    libglib2.0-dev uuid-dev libssl-dev libgpgme-dev \
    libarchive-dev runc golang

git clone --recursive https://github.com/sylabs/singularity.git
cd singularity
git checkout v3.11.4
git submodule update --init --recursive
./mconfig
make -C builddir
sudo make -C builddir install
cd ..

git clone https://github.com/WallauBioinfo/ViralFlow
cd ViralFlow/
micromamba env create -f envs/arm64.yml
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
# AMD64
viralflow build-containers

# ARM64
viralflow build-containers --arch arm64
```

```{note}
This process will download approximately 4.4GB of container images. Ensure you have sufficient disk space and a stable internet connection.
```
