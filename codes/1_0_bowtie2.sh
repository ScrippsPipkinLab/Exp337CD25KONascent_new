#!/bin/bash
#PBS -l nodes=1:ppn=16
#PBS -l mem=8gb
#PBS -l walltime=12:00:00
#PBS -t 1-24

INDIR=/gpfs/group/pipkin/Exp337/trimmed_cb_fq
cd $INDIR

##### Load modules
module load bowtie2
module load samtools
module load bedtools

##### Names
#---- fastq files
trim_fastq_end1=337-${PBS_ARRAYID}_R1_val_1.fq
trim_fastq_end2=337-${PBS_ARRAYID}_R2_val_2.fq

##### Alignment
sam_name=337-${PBS_ARRAYID}.sam

bowtie2_index_m_ensembl="/gpfs/home/hdiao/resources/illumina_ensembl_mus_musculus/Mus_musculus/Ensembl/GRCm38/Sequence/Bowtie2Index/genome"
bowtie2_cmd_1="bowtie2 -p 16 -x $bowtie2_index_m_ensembl -X 2000 --fr"
bowtie2_cmd="$bowtie2_cmd_1 -1 $trim_fastq_end1 -2 $trim_fastq_end2 -S $sam_name"
$bowtie2_cmd


##### Convert/sort/filter
bam_name=337-${PBS_ARRAYID}.bam
bam_name_srt=337-${PBS_ARRAYID}_srt.bam
bam_name_srt_dupr=337-${PBS_ARRAYID}_srt_dupr.bam
samtools view -bS $sam_name > $bam_name
samtools sort $bam_name -o $bam_name_srt
samtools index $bam_name_srt
samtools rmdup -S $bam_name_srt $bam_name_srt_dupr
samtools index $bam_name_srt_dupr
