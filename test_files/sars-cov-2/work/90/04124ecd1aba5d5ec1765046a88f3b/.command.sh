#!/bin/bash -ue
if [[ true  == true ]]; then
    fastp -i Cneg_R1.fq.gz -I Cneg_R2.fq.gz             --detect_adapter_for_pe             --thread 1             -o Cneg.R1.fq.gz             -O Cneg.R2.fq.gz             -h Cneg.fastp.html             -j Cneg.fastp.json             -l 75 -f 0 -t 0             -F 0 -T 0             --cut_front --cut_tail --qualified_quality_phred 20             --dont_eval_duplication
else
    fastp -i Cneg_R1.fq.gz             --thread 1             -o Cneg.SE.fq.gz             -h Cneg.fastp.html             -j Cneg.fastp.json             -l 75 -f 0 -t 0             -F 0 -T 0             --cut_front --cut_tail --qualified_quality_phred 20             --dont_eval_duplication
fi
