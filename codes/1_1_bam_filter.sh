#!/bin/bash
#PBS -l mem=8gb
#PBS -t 1-24

INDIR=/gpfs/group/pipkin/Exp337/srt_bam
cd $INDIR

##### Load modules
module load samtools

##### Names
bam_name_srt=337-${PBS_ARRAYID}_srt.bam
bam_name_srt_dupr=337-${PBS_ARRAYID}_srt_dupr.bam

sam_name_srt_flt=337-${PBS_ARRAYID}_srt_flt.sam
sam_name_srt_dupr_flt=337-${PBS_ARRAYID}_srt_dupr_flt.sam

bam_name_srt_flt=337-${PBS_ARRAYID}_srt_flt.bam
bam_name_srt_dupr_flt=337-${PBS_ARRAYID}_srt_dupr_flt.bam



samtools view $bam_name_srt | awk '{if ($3 ~/'X'/ || $3 ~/'Y'/ || $3 ~/^[0-9]+$/ ) print $0}' > $sam_name_srt_flt
samtools view $bam_name_srt_dupr | awk '{if ($3 ~/'X'/ || $3 ~/'Y'/ || $3 ~/^[0-9]+$/ ) print $0}' > $sam_name_srt_dupr_flt

samtools view -bS -h $sam_name_srt_flt > $bam_name_srt_flt
samtools view -bS -h $sam_name_srt_dupr_flt > $bam_name_srt_dupr_flt

samtools index $bam_name_srt_flt
samtools index $bam_name_srt_dupr_flt
