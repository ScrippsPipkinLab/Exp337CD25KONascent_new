#!/bin/bash
#PBS -l mem=8gb
#PBS -t 1-24

in_dir=/gpfs/group/pipkin/Exp337/cb_fastq
cd $in_dir

module load trimgalore

file1=337-${PBS_ARRAYID}_R1.fastq
file2=337-${PBS_ARRAYID}_R2.fastq

TRIM_CMD="trim_galore --paired --length 24 --stringency 3 $file1 $file2"
$TRIM_CMD

