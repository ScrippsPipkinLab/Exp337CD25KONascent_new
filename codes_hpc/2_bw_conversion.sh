#!/bin/bash
#SBATCH --mem=8gb
#SBATCH --array=1-24

module load trimgalore
module load bowtie2
module load samtools
module load bedtools
module load python
module load ucsc_tools
module load bamtools

BAMDIR=/gpfs/group/pipkin/hdiao/Exp337/srt_flt_bam
CODEDIR=/gpfs/group/pipkin/hdiao/Exp337/codes_hpc

name_array=(WT_0h_rep1 WT_6h_rep1 WT_24h_rep1 WT_48h_rep1 KO_0h_rep1 KO_6h_rep1 KO_24h_rep1 KO_48h_rep1 WT_0h_rep2 WT_6h_rep2 WT_24h_rep2 WT_48h_rep2 KO_0h_rep2 KO_6h_rep2 KO_24h_rep2 KO_48h_rep2 WT_0h_rep3 WT_6h_rep3 WT_24h_rep3 WT_48h_rep3 KO_0h_rep3 KO_6h_rep3 KO_24h_rep3 KO_48h_rep3)

bam_name_srt_flt=$BAMDIR/337_${SLURM_ARRAY_TASK_ID}_srt_flt.bam
bam_name_srt_flt_Rm=$BAMDIR/337_${SLURM_ARRAY_TASK_ID}_srt_flt_Rm.bam
bam_name_srt_flt_RmMm=$BAMDIR/337_${SLURM_ARRAY_TASK_ID}_srt_flt_RmMm.bam
bam_name_srt_flt_RmMm_cor=$BAMDIR/337_${SLURM_ARRAY_TASK_ID}_srt_flt_RmMm_cor.bam

bam_name_srt_flt_r1=$BAMDIR/337_${SLURM_ARRAY_TASK_ID}_srt_flt_r1.bam
bam_name_srt_flt_r2=$BAMDIR/337_${SLURM_ARRAY_TASK_ID}_srt_flt_r2.bam
bam_name_srt_flt_r1_F=$BAMDIR/337_${SLURM_ARRAY_TASK_ID}_srt_flt_r1_F.bam
bam_name_srt_flt_r1_R=$BAMDIR/337_${SLURM_ARRAY_TASK_ID}_srt_flt_r1_R.bam
bam_name_srt_flt_r2_F=$BAMDIR/337_${SLURM_ARRAY_TASK_ID}_srt_flt_r2_F.bam
bam_name_srt_flt_r2_R=$BAMDIR/337_${SLURM_ARRAY_TASK_ID}_srt_flt_r2_R.bam

bam_name_srt_flt_pos=$BAMDIR/337_${SLURM_ARRAY_TASK_ID}_srt_flt_pos.bam
bam_name_srt_flt_neg=$BAMDIR/337_${SLURM_ARRAY_TASK_ID}_srt_flt_neg.bam
bam_name_srt_flt_posneg=$BAMDIR/337_${SLURM_ARRAY_TASK_ID}_srt_flt_posneg.bam
bam_name_srt_flt_posneg_srt=$BAMDIR/337_${SLURM_ARRAY_TASK_ID}_srt_flt_posneg_srt.bam
bdg_name_srt_flt_posneg=$BAMDIR/337_${SLURM_ARRAY_TASK_ID}_srt_flt_posneg.bdg
bdg_name_srt_flt_posneg_srt=$BAMDIR/337_${SLURM_ARRAY_TASK_ID}_srt_flt_posneg_srt.bdg
bdg_name_srt_flt_posneg_srt_chr=$BAMDIR/337_${SLURM_ARRAY_TASK_ID}_srt_flt_posneg_srt.chr.bdg
bw_name_old=$BAMDIR/337_${SLURM_ARRAY_TASK_ID}.bw
bw_name=$BAMDIR/${name_array[${SLURM_ARRAY_TASK_ID}]}.bw

#---- Filter r1 / r2
# Remove reads that are unmapped
# Remove reads whose mates are unmapped
# Keep only reads that are properly paired
samtools view -F 0x4 -h -b $bam_name_srt_flt > $bam_name_srt_flt_Rm 
samtools view -F 0x8 -h -b $bam_name_srt_flt_Rm > $bam_name_srt_flt_RmMm 
samtools view -f 0x2 -h -b $bam_name_srt_flt_RmMm > $bam_name_srt_flt_RmMm_cor 

# Filter out r1 / r2 seperatly
samtools view -f 0x40 -h -b $bam_name_srt_flt_RmMm_cor > $bam_name_srt_flt_r1
samtools view -f 0x80 -h -b $bam_name_srt_flt_RmMm_cor > $bam_name_srt_flt_r2

# Filter out r1 that are positive / r2 that are negative
samtools view -f 0x20 -h -b $bam_name_srt_flt_r1 > $bam_name_srt_flt_r1_F
samtools view -f 0x10 -h -b $bam_name_srt_flt_r2 > $bam_name_srt_flt_r2_R
bamtools merge -in $bam_name_srt_flt_r1_F -in $bam_name_srt_flt_r2_R -out $bam_name_srt_flt_pos

# Filter out r1 that are negative / r2 that are positive
samtools view -f 0x10 -h -b $bam_name_srt_flt_r1 > $bam_name_srt_flt_r1_R
samtools view -f 0x20 -h -b $bam_name_srt_flt_r2 > $bam_name_srt_flt_r2_F
bamtools merge -in $bam_name_srt_flt_r1_R -in $bam_name_srt_flt_r2_F -out $bam_name_srt_flt_neg

#----- Merge strands
bamtools merge -in $bam_name_srt_flt_pos -in $bam_name_srt_flt_neg -out $bam_name_srt_flt_posneg
samtools sort $bam_name_srt_flt_posneg -o $bam_name_srt_flt_posneg_srt
samtools index $bam_name_srt_flt_posneg_srt

#----- Converg ba to bdg
bamCoverage --bam $bam_name_srt_flt_posneg_srt -o $bdg_name_srt_flt_posneg --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 2652783500 --extendReads --outFileFormat bedgraph

#----- Sort
LC_COLLATE=C sort -k1,1 -k2,2n $bdg_name_srt_flt_posneg > $bdg_name_srt_flt_posneg_srt

#----- Covert ensembl chr names
python $CODEDIR/Ensembl_bed_to_ucsc_bed.py $bdg_name_srt_flt_posneg_srt

#----- Bdg to bw
chrom_sizes=/gpfs/group/pipkin/hdiao/ref_resources/mm/release102/GRCm38.genome.sizes.withChr
bedGraphToBigWig $bdg_name_srt_flt_posneg_srt_chr $chrom_sizes $bw_name_old
cp $bw_name_old $bw_name
echo $bw_name_old
echo $bw_name






