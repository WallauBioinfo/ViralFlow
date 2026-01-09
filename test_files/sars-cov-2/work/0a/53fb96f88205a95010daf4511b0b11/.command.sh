#!/bin/bash -ue
if [[ true  == true ]]; then
    fastp -i ART1_R1.fq.gz -I ART1_R2.fq.gz             --detect_adapter_for_pe             --thread 1             -o ART1.R1.fq.gz             -O ART1.R2.fq.gz             -h ART1.fastp.html             -j ART1.fastp.json             -l 75 -f 0 -t 0             -F 0 -T 0             --cut_front --cut_tail --qualified_quality_phred 20             --dont_eval_duplication
else
    fastp -i ART1_R1.fq.gz             --thread 1             -o ART1.SE.fq.gz             -h ART1.fastp.html             -j ART1.fastp.json             -l 75 -f 0 -t 0             -F 0 -T 0             --cut_front --cut_tail --qualified_quality_phred 20             --dont_eval_duplication
fi
