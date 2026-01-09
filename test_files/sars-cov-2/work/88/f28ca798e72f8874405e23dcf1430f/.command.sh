#!/bin/bash -ue
java -jar /app/picard.jar CollectWgsMetrics -I ART3.sorted.bam                                 -R NC_045512.2.fa                                -O wgs -CAP 99999                                 -Q 30                                 -MQ 30                                 --COUNT_UNPAIRED

java -jar /app/picard.jar CollectMultipleMetrics -I ART3.sorted.bam                                      -R NC_045512.2.fa                                     -O metrics
