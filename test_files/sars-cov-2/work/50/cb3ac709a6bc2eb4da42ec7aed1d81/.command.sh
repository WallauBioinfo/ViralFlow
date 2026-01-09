#!/bin/bash -ue
if [[ false  == true ]]; then
    fastp -i ART3_R1.fq.gz -I null             --detect_adapter_for_pe             --thread 1             -o ART3.R1.fq.gz             -O ART3.R2.fq.gz             -h ART3.fastp.html             -j ART3.fastp.json             -l 75 -f 0 -t 0             -F 0 -T 0             --cut_front --cut_tail --qualified_quality_phred 20             --dont_eval_duplication
else
    fastp -i ART3_R1.fq.gz             --thread 1             -o ART3.SE.fq.gz             -h ART3.fastp.html             -j ART3.fastp.json             -l 75 -f 0 -t 0             -F 0 -T 0             --cut_front --cut_tail --qualified_quality_phred 20             --dont_eval_duplication
fi
