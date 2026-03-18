#!/bin/bash -ue
# RUN READ COUNT
bam-readcount -d 50000 -q 30 -w 0 -f NC_045512.2.fa ART3.sorted.bam > ART3.depth25.fa.bc
