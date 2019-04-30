#!/bin/bash
#PBS -t 1-24

##### Load modules
module load python

##### Directory
INDIR=/gpfs/group/pipkin/Exp337/srt_flt_bam/non_dupr/stranded_bam
cd $INDIR

##### Names
ref_gff=/gpfs/group/pipkin/Exp337/codes/GRCm38_exon_rmdup.gff

bam_name_srt_flt_pos_srt=337-${PBS_ARRAYID}_srt_flt_pos_srt.bam
bam_name_srt_flt_neg_srt=337-${PBS_ARRAYID}_srt_flt_neg_srt.bam

pos_count=337-${PBS_ARRAYID}_pos.count
neg_count=337-${PBS_ARRAYID}_neg.count

##### Count
htseq-count \
-f bam -t exon -i ID -m intersection-nonempty \
$bam_name_srt_flt_pos_srt $ref_gff > $pos_count

htseq-count \
-f bam -t exon -i ID -m intersection-nonempty \
$bam_name_srt_flt_neg_srt $ref_gff > $neg_count
