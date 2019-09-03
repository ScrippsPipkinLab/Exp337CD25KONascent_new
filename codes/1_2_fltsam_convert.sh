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

sam_name_srt_flt_h=337-${PBS_ARRAYID}_srt_flt_h.sam
sam_name_srt_dupr_flt_h=337-${PBS_ARRAYID}_srt_dupr_flt_h.sam

bam_name_srt_flt=337-${PBS_ARRAYID}_srt_flt.bam
bam_name_srt_dupr_flt=337-${PBS_ARRAYID}_srt_dupr_flt.bam

sp_header=337-${PBS_ARRAYID}.sp_header

samtools view -H $bam_name_srt > $sp_header

cat $sp_header $sam_name_srt_flt > $sam_name_srt_flt_h
cat $sp_header $sam_name_srt_dupr_flt > $sam_name_srt_dupr_flt_h

samtools view -bS -h $sam_name_srt_flt_h > $bam_name_srt_flt
samtools view -bS -h $sam_name_srt_dupr_flt_h > $bam_name_srt_dupr_flt

samtools index $bam_name_srt_flt
samtools index $bam_name_srt_dupr_flt

