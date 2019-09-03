#!/bin/bash
#PBS -t 1-24

##### Load modules
module load python

##### Directory
INDIR=/gpfs/group/pipkin/Exp337/srt_flt_bam/dupr/stranded_bam
cd $INDIR

##### Names
ref_gff=/gpfs/group/pipkin/Exp337/codes/GRCm38_exon_rmdup_srt_cb_srt_dupr.gff

bam_name_srt_dupr_flt_pos_srt=337-${PBS_ARRAYID}_srt_dupr_flt_pos_srt.bam
bam_name_srt_dupr_flt_neg_srt=337-${PBS_ARRAYID}_srt_dupr_flt_neg_srt.bam

dupr_pos_count=337-${PBS_ARRAYID}_dupr_pos.count
dupr_neg_count=337-${PBS_ARRAYID}_dupr_neg.count

##### Count
htseq-count \
-f bam -t exon -i ID -m intersection-nonempty \
$bam_name_srt_dupr_flt_pos_srt $ref_gff > $dupr_pos_count

htseq-count \
-f bam -t exon -i ID -m intersection-nonempty \
$bam_name_srt_dupr_flt_neg_srt $ref_gff > $dupr_neg_count

