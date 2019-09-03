########## DEseq2 for Nascent RNA-seq count table ##########
# Author: Huitian (Yolanda) Diao
# April 30th, 2019
# Dependencies:
# -| Compiled count table <- 2_1_collect_count_files.py
# --| HT-seq count file <- HT-seq count <- 2_0_ht-seq_count.dupr.sh
# ---| Strand seperated bam file <- 1_3_splitStrands_bamTobdg_add-chr_bdgTobw_new.sh
# ---| Exon reference GFF3 file
# ----| Biomart: http://useast.ensembl.org/biomart/martview/a4b3d3135f51db16df0294bef537f063
# ----| BioMart-out-csv_To_Gff.sh  &  Gff_rmdup.py


######################################## Imports ########################################
library(DESeq2)

###--- Config
wk.dir <- "/Volumes/EXP337/Exp337CD25KONascent/2_count/3_DE-seq"
op.dir <- "/Volumes/EXP337/Exp337CD25KONascent/2_count/3_DE-seq"
setwd(wk.dir)

######################################## Main ########################################
count.file <- "/Volumes/EXP337/Exp337CD25KONascent/2_count/2_compiled_csv/Exp337_all.csv"
row.names <- c("KO_0h_rep1_dupr", "KO_0h_rep2_dupr", "KO_0h_rep3_dupr", "KO_24h_rep1_dupr", 
               "KO_24h_rep2_dupr", "KO_24h_rep3_dupr", "KO_48h_rep1_dupr", "KO_48h_rep2_dupr", 
               "KO_48h_rep3_dupr", "KO_6h_rep1_dupr", "KO_6h_rep2_dupr", "WT_0h_rep1_dupr", 
               "WT_0h_rep2_dupr", "WT_0h_rep3_dupr", "WT_24h_rep1_dupr", "WT_24h_rep2_dupr", 
               "WT_24h_rep3_dupr", "WT_48h_rep1_dupr", "WT_48h_rep2_dupr", "WT_48h_rep3_dupr", 
               "WT_6h_rep1_dupr", "WT_6h_rep2_dupr", "WT_6h_rep3_dupr")
conds <- c("KO_0h", "KO_0h", "KO_0h", "KO_24h", "KO_24h", "KO_24h", "KO_48h", "KO_48h", "KO_48h", "KO_6h", "KO_6h", 
           "WT_0h", "WT_0h", "WT_0h", "WT_24h", "WT_24h", "WT_24h", "WT_48h", "WT_48h", "WT_48h", "WT_6h", "WT_6h", "WT_6h")






