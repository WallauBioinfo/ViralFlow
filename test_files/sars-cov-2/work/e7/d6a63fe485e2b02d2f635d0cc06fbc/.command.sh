#!/bin/bash -ue
java -jar /app/picard.jar CollectWgsMetrics -I Cneg.sorted.bam                                 -R NC_045512.2.fa                                -O wgs -CAP 99999                                 -Q 30                                 -MQ 30                                 

java -jar /app/picard.jar CollectMultipleMetrics -I Cneg.sorted.bam                                      -R NC_045512.2.fa                                     -O metrics
