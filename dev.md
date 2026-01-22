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
viralflow build_containers
```

### ARM64

```bash
git clone -b develop https://github.com/WallauBioinfo/ViralFlow.git
cd ViralFlow/
micromamba env create -f envs/arm64.yml
micromamba activate viralflow
pip install -e .
# include apptainer / singularity config here
viralflow build_containers --arch arm64
```

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