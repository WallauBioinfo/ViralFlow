#!/bin/bash

# get bash script location
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
cd $SCRIPT_DIR

# user input
organism_name=$1 # Dengue
organism_refseq_code=$2 # NC_001474.2

# hardcoded paths
SNPEFF_CTNR="snpeff:5.0.sif"
EFETCH_CTNR="edirect:1.1.0.sif"

# Detect the correct snpEff path (varies between systems: 5.0-2 or 5.0-3)
if [ -d "$SNPEFF_CTNR/opt/conda/share/snpeff-5.0-3" ]; then
    SNPEFF_PATH="/opt/conda/share/snpeff-5.0-3/"
elif [ -d "$SNPEFF_CTNR/opt/conda/share/snpeff-5.0-2" ]; then
    SNPEFF_PATH="/opt/conda/share/snpeff-5.0-2/"
else
    echo "ERROR: Could not find snpEff installation directory in container"
    exit 1
fi
echo "Using snpEff path: $SNPEFF_PATH"

echo "@ adding new entry..."
echo -e "# $organism_name, version $organism_refseq_code" >> $SNPEFF_CTNR/$SNPEFF_PATH/snpEff.config
echo -e "$organism_refseq_code.genome: $organism_name" >> $SNPEFF_CTNR/$SNPEFF_PATH/snpEff.config
echo -e "$organism_refseq_code.has_cds: true" >> $SNPEFF_CTNR/$SNPEFF_PATH/snpEff.config
echo -e "$organism_refseq_code.codonTable: Standard" >> $SNPEFF_CTNR/$SNPEFF_PATH/snpEff.config

# build the directory at DB
mkdir -p $SNPEFF_CTNR/$SNPEFF_PATH/data/$organism_refseq_code

# download fasta
echo "@ downloading fasta..."
#sudo singularity exec $EFETCH_CTNR efetch -db nucleotide -id $organism_refseq_code -format gb > $SNPEFF_CTNR/$SNPEFF_PATH/data/$organism_refseq_code/genes.gbk
singularity exec --fakeroot $EFETCH_CTNR efetch -db nucleotide -id $organism_refseq_code -format gb > $SNPEFF_CTNR/$SNPEFF_PATH/data/$organism_refseq_code/genes.gbk

# build database
echo "@ rebuild database"
#sudo singularity exec --writable $SNPEFF_CTNR snpEff build -genbank -v $organism_refseq_code
singularity exec --fakeroot --writable $SNPEFF_CTNR snpEff build -genbank -v $organism_refseq_code

# update catalog
echo "@ update snpeff database catalog..."
singularity exec $SNPEFF_CTNR snpEff databases > snpEff_DB.catalog
