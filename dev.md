# ViralFlow - Development Documentation

For computing environments with ARM64 architecture the user should inform the `--arch arm64` flag into some commands.

## Installation

### AMD64

```bash
git clone -b develop https://github.com/WallauBioinfo/ViralFlow.git
cd ViralFlow/
micromamba env create -f envs/amd64.yml
micromamba activate viralflow
pip install -e .
sudo ln -s /usr/bin/unsquashfs /usr/local/bin/unsquashfs
viralflow build-containers
```

### ARM64

```bash
# install singularity
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

# install viralflow wrapper
git clone -b develop https://github.com/WallauBioinfo/ViralFlow.git
cd ViralFlow/
micromamba env create -f envs/arm64.yml
micromamba activate viralflow
pip install -e .

# build containers
viralflow build-containers --arch arm64
```

## Development quality checks

Python 3.12 is the supported development version. For a lightweight setup, use
`uv` to create the Python environment and install ViralFlow and `pre-commit`:

```bash
uv venv --python 3.12
source .venv/bin/activate
uv pip install -e . "pre-commit==4.6.0"
```

The Python test hook uses `uv` to create and cache its own Python 3.12
environment. Nextflow linting and the pre-push tests still require `nextflow` and
`nf-test` on `PATH`. Alternatively, both architecture-specific development
environments include Python 3.12, `pre-commit`, Nextflow, and nf-test.

Install both Git hook stages from the repository root:

```bash
pre-commit install --hook-type pre-commit --hook-type pre-push
```

The pre-commit stage runs repository hygiene checks, Ruff linting and formatting,
Nextflow linting, and the Python unit tests. Hooks that modify files will fail the
first run so the updated files can be reviewed and staged again.

Run the complete commit-time suite manually with:

```bash
pre-commit run --all-files
```

The pre-push stage runs the 14 Nextflow tests that do not require local container
images. Run it manually with:

```bash
pre-commit run --all-files --hook-stage pre-push
```

The BCFtools, container metadata, and Nanopore truth tests remain manual because
they require Singularity and a locally built `vfnext/containers/baseContainer.sif`
(and the truth test also uses the Clair3 container):

```bash
cd vfnext
NXF_VER=26.04.6 nf-test test \
  tests/workflows/bcftools-fixture.nf.test \
  tests/workflows/metadata-fixture.nf.test \
  --ci
NXF_VER=26.04.6 nf-test test integration_tests/nanopore-truth.nf.test --ci
```

For an emergency-only bypass, use `git commit --no-verify` or
`git push --no-verify`, then run the skipped hook stage manually before opening
or updating a pull request.

### Customizing snpEff catalog

#### AMD64

```bash
viralflow add-entry-to-snpeff --org-name Dengue --genome-code NC_001474.2
```

#### ARM64

```bash
viralflow add-entry-to-snpeff --org-name Dengue --genome-code NC_001474.2 --arch arm64
```

### Updating pangolin

```bash
viralflow update-pangolin

viralflow update-pangolin-data
```

### Running ViralFlow

```bash
viralflow run --params-file test_files/sars-cov-2.params
```

```bash
viralflow run --params-file test_files/denv.params
```
