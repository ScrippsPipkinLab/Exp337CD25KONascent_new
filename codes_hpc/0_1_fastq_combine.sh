#!/bin/bash

in_dir=/gpfs/group/pipkin/Exp337
cd $in_dir

for i in $(seq 1 24)
do
  cat 337-${i}_*_R1_001.fastq >> 337-${i}_R1.fastq
  cat 337-${i}_*_R2_001.fastq >> 337-${i}_R2.fastq
done
