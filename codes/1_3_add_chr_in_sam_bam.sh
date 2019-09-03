#!/bin/bash
#PBS -l mem=8gb
#PBS -t 1-24

INDIR=/gpfs/group/pipkin/Exp337/srt_flt_bam
cd $INDIR

##### Load modules
module load samtools

##### Names
bam_name_srt_flt=337-${PBS_ARRAYID}_srt_flt.bam
bam_name_srt_dupr_flt=337-${PBS_ARRAYID}_srt_dupr_flt.bam

sam_name_srt_flt_chr=337-${PBS_ARRAYID}_srt_flt_chr.sam
sam_name_srt_dupr_flt_chr=337-${PBS_ARRAYID}_srt_dupr_flt_chr.sam

bam_name_srt_flt_chr=337-${PBS_ARRAYID}_srt_flt_chr.bam
bam_name_srt_dupr_flt_chr=337-${PBS_ARRAYID}_srt_dupr_flt_chr.bam

samtools view -h $bam_name_srt_flt \
| awk  \
'{if ($3 ~/'X'/ || $3 ~/'Y'/ || $3 ~/^[0-9]+$/ ) \
{print $1 "\t" $2 "\t" "chr"$3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t" $13 "\t" $14 "\t" $15 "\t" $16 "\t" $17 "\t" $18 "\t" $19 "\t" $20 "\t" $21} \
else {print $0}}' \
> $sam_name_srt_flt_chr


samtools view -h $bam_name_srt_dupr_flt \
| awk  \
'{if ($3 ~/'X'/ || $3 ~/'Y'/ || $3 ~/^[0-9]+$/ ) \
{print $1 "\t" $2 "\t" "chr"$3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t" $13 "\t" $14 "\t" $15 "\t" $16 "\t" $17 "\t" $18 "\t" $19 "\t" $20 "\t" $21} \
else {print $0}}' \
> $sam_name_srt_dupr_flt_chr

samtools view -h -bS $sam_name_srt_flt_chr > $bam_name_srt_flt_chr

samtools view -h -bS $sam_name_srt_dupr_flt_chr > $bam_name_srt_dupr_flt_chr
