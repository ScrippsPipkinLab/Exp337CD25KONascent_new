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
count.file <- "/Volumes/EXP337/Exp337CD25KONascent/2_count/2_Compiled_csv/Exp337_dupr_all_count_c5.csv"
row.names <- c("KO_0h_rep1_dupr", "KO_0h_rep2_dupr", "KO_0h_rep3_dupr", "KO_24h_rep1_dupr", 
               "KO_24h_rep2_dupr", "KO_24h_rep3_dupr", "KO_48h_rep1_dupr", "KO_48h_rep2_dupr", 
               "KO_48h_rep3_dupr", "KO_6h_rep1_dupr", "KO_6h_rep2_dupr", "WT_0h_rep1_dupr", 
               "WT_0h_rep2_dupr", "WT_0h_rep3_dupr", "WT_24h_rep1_dupr", "WT_24h_rep2_dupr", 
               "WT_24h_rep3_dupr", "WT_48h_rep1_dupr", "WT_48h_rep2_dupr", "WT_48h_rep3_dupr", 
               "WT_6h_rep1_dupr", "WT_6h_rep2_dupr", "WT_6h_rep3_dupr")
conds <- c("KO_0h", "KO_0h", "KO_0h", "KO_24h", "KO_24h", "KO_24h", "KO_48h", "KO_48h", "KO_48h", "KO_6h", "KO_6h", 
           "WT_0h", "WT_0h", "WT_0h", "WT_24h", "WT_24h", "WT_24h", "WT_48h", "WT_48h", "WT_48h", "WT_6h", "WT_6h", "WT_6h")

#---- Read input file
DeseqData <- read.csv(count.file, header=TRUE, row.names=1)
#---- Meta data
DeseqDesign <- data.frame(row.names=row.names,condition=as.factor(conds))
#---- Generate DEseq matrix
DESmat <- DESeqDataSetFromMatrix(countData = DeseqData, colData = DeseqDesign,design = ~ condition)
#---- Run DEseq
DESmat <- DESeq(DESmat)
#---- Output
setwd(op.dir)

WT_6_0 <- results(DESmat, contrast=c("condition", "WT_6h", "WT_0h"))
WT_24_0 <- results(DESmat, contrast=c("condition", "WT_24h", "WT_0h"))
WT_48_0 <- results(DESmat, contrast=c("condition", "WT_48h", "WT_0h"))
write.csv(WT_6_0, "WT-6h_vs_WT-0h.csv")
write.csv(WT_24_0, "WT-24h_vs_WT-0h.csv")
write.csv(WT_48_0, "WT-48h_vs_WT-0h.csv")

KO_6_0 <- results(DESmat, contrast=c("condition", "KO_6h", "KO_0h"))
KO_24_0 <- results(DESmat, contrast=c("condition", "KO_24h", "KO_0h"))
KO_48_0 <- results(DESmat, contrast=c("condition", "KO_48h", "KO_0h"))
write.csv(KO_6_0, "KO-6h_vs_KO-0h.csv")
write.csv(KO_24_0, "KO-24h_vs_KO-0h.csv")
write.csv(KO_48_0, "KO-48h_vs_KO-0h.csv")

WT0_KO0 <- results(DESmat, contrast=c("condition", "WT_0h", "KO_0h"))
WT6_KO6 <- results(DESmat, contrast=c("condition", "WT_6h", "KO_6h"))
WT24_KO24 <- results(DESmat, contrast=c("condition", "WT_24h", "KO_24h"))
WT48_KO48 <- results(DESmat, contrast=c("condition", "WT_48h", "KO_48h"))
write.csv(WT0_KO0, "WT-0h_vs_KO-0h.csv")
write.csv(WT6_KO6, "WT-6h_vs_KO-6h.csv")
write.csv(WT24_KO24, "WT-24h_vs_KO-24h.csv")
write.csv(WT48_KO48, "WT-48h_vs_KO-48h.csv")