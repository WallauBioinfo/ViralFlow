#!/bin/bash -ue
# IVAR STEP 1 ----------------------------------------------------------------
samtools mpileup -aa -d 50000 --reference NC_045512.2.fa -a -B ART3.sorted.bam |        ivar variants -p ART3 -q 30 -t 0.05

# IVAR STEP 2 ----------------------------------------------------------------
samtools mpileup -aa -d 50000 --reference NC_045512.2.fa -a -B ART3.sorted.bam |        ivar consensus -p ART3 -q 30 -t 0 -m 25 -n N -c 0.51

# IVAR STEP 3 ----------------------------------------------------------------
samtools mpileup -aa -d 50000 --reference NC_045512.2.fa -a -B ART3.sorted.bam |        ivar consensus -p ART3.ivar060 -q 30 -t 0.60 -n N -m 25 -c 0.51
# EDIT FILE NAMES
mv ART3.fa ART3.depth25.fa
mv ART3.ivar060.fa ART3.depth25.amb.fa
sed -i -e 's/>.*/>ART3/g' ART3.depth25.fa
sed -i -e 's/>.*/>ART3/g' ART3.depth25.amb.fa
