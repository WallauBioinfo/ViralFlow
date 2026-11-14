#!/bin/bash

set -euo pipefail

# get bash script location
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
cd "$SCRIPT_DIR"

# user input
organism_name=${1:?Usage: add_entries_SnpeffDB.sh <organism_name> <refseq_code> <arch>}
organism_refseq_code=${2:?Usage: add_entries_SnpeffDB.sh <organism_name> <refseq_code> <arch>}
arch=${3:?Usage: add_entries_SnpeffDB.sh <organism_name> <refseq_code> <arch>}

# hardcoded paths
SNPEFF_CTNR="snpeff:5.0.sif"
EFETCH_CTNR="edirect:1.1.0.sif"

if [ "$arch" == "amd64" ]; then
    # Detect the correct snpEff path (varies between systems: 5.0-2 or 5.0-3)
    if [ -d "$SNPEFF_CTNR/opt/conda/share/snpeff-5.0-3" ]; then
        SNPEFF_PATH="/opt/conda/share/snpeff-5.0-3"
    elif [ -d "$SNPEFF_CTNR/opt/conda/share/snpeff-5.0-2" ]; then
        SNPEFF_PATH="/opt/conda/share/snpeff-5.0-2"
    else
        echo "ERROR: Could not find snpEff installation directory inside $SNPEFF_CTNR for amd64"
        exit 1
    fi
elif [ "$arch" == "arm64" ]; then
    if [ -d "$SNPEFF_CTNR/usr/local/bin/mm/share/snpeff-5.0-3" ]; then
        SNPEFF_PATH="/usr/local/bin/mm/share/snpeff-5.0-3"
    elif [ -d "$SNPEFF_CTNR/usr/local/bin/mm/share/snpeff-5.0-2" ]; then
        SNPEFF_PATH="/usr/local/bin/mm/share/snpeff-5.0-2"
    else
        echo "ERROR: Could not find snpEff installation directory inside $SNPEFF_CTNR for arm64"
        exit 1
    fi
else
    echo "ERROR: Unsupported architecture: $arch. Expected 'amd64' or 'arm64'."
    exit 1
fi

echo "Using snpEff path: $SNPEFF_PATH"

TMP_DIR=$(mktemp -d)
trap 'rm -rf "$TMP_DIR"' EXIT

# GenBank assembly records may omit the sequence unless gbwithparts is used.
echo "@ downloading GenBank..."
singularity exec --fakeroot "$EFETCH_CTNR" efetch -db nucleotide \
    -id "$organism_refseq_code" -format gbwithparts > "$TMP_DIR/genes.gbk"
if ! grep -q '^ORIGIN' "$TMP_DIR/genes.gbk"; then
    echo "ERROR: Downloaded GenBank record contains no sequence."
    exit 1
fi

echo "@ adding new entry..."
CONFIG_FILE="$SNPEFF_CTNR/$SNPEFF_PATH/snpEff.config"
if ! grep -Fq "$organism_refseq_code.genome:" "$CONFIG_FILE"; then
    {
        printf '# %s, version %s\n' "$organism_name" "$organism_refseq_code"
        printf '%s.genome: %s\n' "$organism_refseq_code" "$organism_name"
        printf '%s.has_cds: true\n' "$organism_refseq_code"
        printf '%s.codonTable: Standard\n' "$organism_refseq_code"
    } >> "$CONFIG_FILE"
fi

# build the directory at DB
DATA_DIR="$SNPEFF_CTNR/$SNPEFF_PATH/data/$organism_refseq_code"
mkdir -p "$DATA_DIR"
cp "$TMP_DIR/genes.gbk" "$DATA_DIR/genes.gbk"

# build database
echo "@ rebuild database"
singularity exec --fakeroot --writable "$SNPEFF_CTNR" snpEff build -genbank -v "$organism_refseq_code"
if [ ! -s "$DATA_DIR/snpEffectPredictor.bin" ]; then
    echo "ERROR: snpEff did not produce the custom database."
    exit 1
fi

# update catalog
echo "@ update snpeff database catalog..."
singularity exec "$SNPEFF_CTNR" snpEff databases > "$TMP_DIR/snpEff_DB.catalog"
mv "$TMP_DIR/snpEff_DB.catalog" snpEff_DB.catalog
