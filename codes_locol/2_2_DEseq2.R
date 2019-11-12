########## DEseq2 for Nascent RNA-seq count table ##########
# Author: Huitian (Yolanda) Diao
# April 30th, 2019
# Dependencies:
# -| Compiled count table <- 2_0_collect_featureCounts.py
# --| Count files <- featureCount <- 2_0_featureCount.sh
# ---| Strand seperated bam file <- 1_3_splitStrands_bamTobdg_add-chr_bdgTobw_new.sh
# ---| Exon reference GFF3 file
# ----| Biomart: http://useast.ensembl.org/biomart/martview/a4b3d3135f51db16df0294bef537f063
# ----| BioMart-out-csv_To_Gff.sh  &  Gff_rmdup.py

######################################## Imports ########################################
library(DESeq2)
library(dplyr)
library(tidyverse)


######################################## Functions ########################################
# Replace vector elements into new names
cvt_spNames <- function(vec_x, spNumber, spName, spCond){
  vec_x_out <- c()
  vec_cond_out <- c()
  spNumber <- paste(spNumber, "_", sep="")
  for (i in vec_x) {
    new_i <- gsub("_dupr", "", i)
    new_i <- gsub("_f", "", new_i)
    new_i <- gsub("_r", "", new_i)
    new_i <- paste(new_i, "_", sep="")
    if (new_i %in% spNumber){
      i_idx <- match(new_i, spNumber)
      vec_x_out <- c(vec_x_out, spName[i_idx])
      vec_cond_out <- c(vec_cond_out, spCond[i_idx])
    } else {
      print(paste(as.character(i), "not found!!"))
    }
  }
  return(list(vec_x_out, vec_cond_out))
}

# Convert names to conditions
cvt_nameToCond <- function(vec_x){
  out_vec <- c()
  for (i in vec_x) {
    out_vec <- c(out_vec, strsplit(i, "_rep")[[1]][1])
  }
  return(out_vec)
}


######################################## Config ########################################
##########---------- Dupr ----------##########
#wk.dir <- "/Volumes/Yolanda/Exp337CD25KONascent/3_DE-seq/dupr"
#count.file.f <- "/Volumes/Yolanda/Exp337CD25KONascent/2_count/dupr/2_compiled_csv/Exp337_dupr_F_count_c10.csv"
#count.file.r <- "/Volumes/Yolanda/Exp337CD25KONascent/2_count/dupr/2_compiled_csv/Exp337_dupr_R_count_c10.csv"
#sp.info.file <- "/Volumes/Yolanda/Exp337CD25KONascent/Info/sample_sheet.csv"
#setwd(wk.dir)

##########---------- nonDupr ----------##########
wk.dir <- "/Volumes/Yolanda/Exp337CD25KONascent/3_DE-seq/nondupr"
count.file.f <- "/Volumes/Yolanda/Exp337CD25KONascent/2_count/nondupr/2_compiled_csv/Exp337_F_count_c10.csv"
count.file.r <- "/Volumes/Yolanda/Exp337CD25KONascent/2_count/nondupr/2_compiled_csv/Exp337_R_count_c10.csv"
sp.info.file <- "/Volumes/Yolanda/Exp337CD25KONascent/Info/sample_sheet.csv"
setwd(wk.dir)

######################################## Main ########################################
sp.tb <- read_csv(sp.info.file)
f.tb <- read_csv(count.file.f)
r.tb <- read_csv(count.file.r)

###----- Convert sample number to sample names
f.oldnames <- colnames(f.tb)[2:length(colnames(f.tb))]
f.names.conds <- cvt_spNames(f.oldnames, sp.tb$sp_order, sp.tb$sp_name, sp.tb$sp_cond)
colnames(f.tb) <- c("name", f.names.conds[[1]])
f.tb$KO_6h_rep3 <- NULL

r.oldnames <- colnames(r.tb)[2:length(colnames(r.tb))]
r.names.conds <- cvt_spNames(r.oldnames, sp.tb$sp_order, sp.tb$sp_name, sp.tb$sp_cond)
colnames(r.tb) <- c("name", r.names.conds[[1]])
r.tb$KO_6h_rep3 <- NULL

###----- Check if f and r names match
match <- (f.names.conds[[1]] == r.names.conds[[1]])
if (FALSE %in% match) {
  fr.tb <- NULL
  print("Forward and reverse strand count summary files do not match (sample orders different)")
} else {
  fr.tb <- bind_rows(f.tb, r.tb) %>% column_to_rownames(var = "name")
}

###----- Construct DEseq objects
DeseqDesign <- data.frame(names = colnames(fr.tb), condition = as.factor(cvt_nameToCond(colnames(fr.tb))))
DESmat <- DESeqDataSetFromMatrix(countData = fr.tb, colData = DeseqDesign, design = ~ condition)
DESmat <- DESeq(DESmat)


###----- Write outputs
cts.use <- list(c("KO_6h","KO_0h"), c("KO_24h", "KO_0h"), c("KO_48h", "KO_0h"),
                c("WT_6h","WT_0h"), c("WT_24h", "WT_0h"), c("WT_48h", "WT_0h"),
                c("WT_0h", "KO_0h"), c("WT_6h", "KO_6h"), 
                c("WT_24h", "KO_24h"), c("WT_48h", "KO_48h"))
for (i in c(1 : length(cts.use))) {
  cts <- cts.use[[i]]
  cts.name <- paste(cts[1], "_vs_", cts[2], ".csv", sep = "")
  cts.res <- results(DESmat, contrast = c("condition", cts))
  write.csv(cts.res, cts.name)
}






