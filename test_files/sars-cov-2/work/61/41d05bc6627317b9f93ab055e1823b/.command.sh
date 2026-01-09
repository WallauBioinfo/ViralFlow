#!/bin/bash -ue
if [[ true  == true ]]; then
    fastp -i ART2_R1.fq.gz -I ART2_R2.fq.gz             --detect_adapter_for_pe             --thread 1             -o ART2.R1.fq.gz             -O ART2.R2.fq.gz             -h ART2.fastp.html             -j ART2.fastp.json             -l 75 -f 0 -t 0             -F 0 -T 0             --cut_front --cut_tail --qualified_quality_phred 20             --dont_eval_duplication
else
    fastp -i ART2_R1.fq.gz             --thread 1             -o ART2.SE.fq.gz             -h ART2.fastp.html             -j ART2.fastp.json             -l 75 -f 0 -t 0             -F 0 -T 0             --cut_front --cut_tail --qualified_quality_phred 20             --dont_eval_duplication
fi
