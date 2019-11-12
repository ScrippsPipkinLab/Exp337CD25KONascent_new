#!/bin/bash
#PBS -t 1-24

##### Load modules
module load subread

##### Names
ref_gff_f=/gpfs/group/pipkin/Exp337/codes/GRCm38_exon_rmdup_srt_cb_srt_dupr_f.gff
ref_gff_r=/gpfs/group/pipkin/Exp337/codes/GRCm38_exon_rmdup_srt_cb_srt_dupr_r.gff

#####---------- Dupr
INDIR=/gpfs/group/pipkin/Exp337/srt_flt_bam/dupr/stranded_bam
WKDIR=/gpfs/group/pipkin/hdiao/Exp337/featureCounts_dupr
cd $WKDIR

bam_pos_psrt=337-${PBS_ARRAYID}_pos_pSrt.bam
bam_neg_psrt=337-${PBS_ARRAYID}_neg_pSrt.bam

pos_count=337-${PBS_ARRAYID}_dupr_r.txt
neg_count=337-${PBS_ARRAYID}_dupr_f.txt

##### Count
featureCounts -T 8 -t "exon" -f -g ID -p -O -a $ref_gff_r -o $pos_count ${INDIR}/${bam_pos_psrt}
featureCounts -T 8 -t "exon" -f -g ID -p -O -a $ref_gff_f -o $neg_count ${INDIR}/${bam_neg_psrt}

#####---------- Non Dupr
INDIR=/gpfs/group/pipkin/Exp337/srt_flt_bam/non_dupr/stranded_bam
WKDIR=/gpfs/group/pipkin/hdiao/Exp337/featureCounts_nondupr
cd $WKDIR

bam_pos_psrt=337-${PBS_ARRAYID}_pos_pSrt.bam
bam_neg_psrt=337-${PBS_ARRAYID}_neg_pSrt.bam

pos_count=337-${PBS_ARRAYID}_r.txt
neg_count=337-${PBS_ARRAYID}_f.txt

##### Count
featureCounts -T 8 -t "exon" -f -g ID -p -O -a $ref_gff_r -o $pos_count ${INDIR}/${bam_pos_psrt}
featureCounts -T 8 -t "exon" -f -g ID -p -O -a $ref_gff_f -o $neg_count ${INDIR}/${bam_neg_psrt}
