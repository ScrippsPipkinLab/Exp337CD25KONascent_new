#!/bin/bash
#PBS -l mem=8gb
#PBS -t 1-24

INDIR=/gpfs/group/pipkin/Exp337/srt_flt_bam
cd $INDIR

##### Load modules
module load samtools
module load python
module load ucsc_tools

##### Names
bam_name_srt_flt=337-${PBS_ARRAYID}_srt_flt.bam
bam_name_srt_dupr_flt=337-${PBS_ARRAYID}_srt_dupr_flt.bam

bdg_name_srt_flt=337-${PBS_ARRAYID}_srt_flt.bdg
bdg_name_srt_dupr_flt=337-${PBS_ARRAYID}_srt_dupr_flt.bdg

bdg_name_srt_flt_chr=337-${PBS_ARRAYID}_srt_flt_chr.bdg
bdg_name_srt_dupr_flt_chr=337-${PBS_ARRAYID}_srt_dupr_flt_chr.bdg

bdg_name_srt_flt_chr_srt=337-${PBS_ARRAYID}_srt_flt_chr_srt.bdg
bdg_name_srt_dupr_flt_chr_srt=337-${PBS_ARRAYID}_srt_dupr_flt_chr_srt.bdg

bw_name_srt_flt_chr_srt=337-${PBS_ARRAYID}_srt_flt.bw
bw_name_srt_dupr_flt_chr_srt=337-${PBS_ARRAYID}_srt_dupr_flt.bw

#---- Bam to bdg
bamCoverage --bam $bam_name_srt_flt -o $bdg_name_srt_flt  --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 2652783500 --extendReads --outFileFormat bedgraph 
bamCoverage --bam $bam_name_srt_dupr_flt -o $bdg_name_srt_dupr_flt  --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 2652783500 --extendReads --outFileFormat bedgraph 

#---- Bdg add chr
awk '{print "chr"$1 "\t" $2 "\t" $3 "\t" $4}' $bdg_name_srt_flt > $bdg_name_srt_flt_chr
awk '{print "chr"$1 "\t" $2 "\t" $3 "\t" $4}' $bdg_name_srt_dupr_flt > $bdg_name_srt_dupr_flt_chr

#---- Sort bdg
LC_COLLATE=C sort -k1,1 -k2,2n $bdg_name_srt_flt_chr > $bdg_name_srt_flt_chr_srt
LC_COLLATE=C sort -k1,1 -k2,2n $bdg_name_srt_dupr_flt_chr > $bdg_name_srt_dupr_flt_chr_srt

#---- Bdg to bw
bedGraphToBigWig $bdg_name_srt_flt_chr_srt mm10.chrom.sizes $bw_name_srt_flt_chr_srt
bedGraphToBigWig $bdg_name_srt_dupr_flt_chr_srt mm10.chrom.sizes $bw_name_srt_dupr_flt_chr_srt
