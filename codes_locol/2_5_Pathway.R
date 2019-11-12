########## Pathway analysis ##########
# Author: Huitian (Yolanda) Diao
# April 30th, 2019
# Dependencies:
# -| DEseq2 output
# --| Compiled count table <- 2_1_collect_count_files.py
# ---| HT-seq count file <- HT-seq count <- 2_0_ht-seq_count.dupr.sh
# ----| Strand seperated bam file <- 1_3_splitStrands_bamTobdg_add-chr_bdgTobw_new.sh
# ----| Exon reference GFF3 file
# -----| Biomart: http://useast.ensembl.org/biomart/martview/a4b3d3135f51db16df0294bef537f063
# -----| BioMart-out-csv_To_Gff.sh  &  Gff_rmdup.py

######################################## Imports ########################################
library(ggplot2)
library(ggrepel)
library("clusterProfiler")
#BiocManager::install("org.Mm.eg.db")
library("org.Mm.eg.db")
library(dplyr)
library(tidyverse)

GO_run <- function(genes.i, i){
  #file.i <-  "/Volumes/Yolanda/Exp334CD25KOSc/4_D3E/0_clusterComparison/2_updn/C0_vs_C1_fix_pval0.05_dn.csv" # For testing
  #log2fc_cutoff <- 1.5 # For testing
  print(paste(i, "    Gene number: ", as.character(length(genes.i), sep="")))
  genes.i.id <- AnnotationDbi::select(org.Mm.eg.db, genes.i, c("ENTREZID"), "ALIAS")
  
  egoBP <- enrichGO(gene = genes.i.id$ENTREZID, keyType = 'ENTREZID', OrgDb = org.Mm.eg.db, ont = "BP", pvalueCutoff = 0.05, readable = TRUE) #pAdjustMethod = "none"
  egoCC <- enrichGO(gene = genes.i.id$ENTREZID, keyType = 'ENTREZID', OrgDb = org.Mm.eg.db, ont = "CC", pvalueCutoff = 0.05, readable = TRUE)
  egoMF <- enrichGO(gene = genes.i.id$ENTREZID, keyType = 'ENTREZID', OrgDb = org.Mm.eg.db, ont = "MF", pvalueCutoff = 0.05, readable = TRUE)
  
  # Dotplot visualization
  if (!is.null(egoBP)){
    pdf.name <- paste(i,"_BP_dotplot.pdf",sep="")
    csv.name <- paste(i,"_BP_dotplot.csv",sep="")
    write.csv(egoBP@result, file=csv.name, row.names=FALSE)
    egoBP.dotplot <- dotplot(egoBP, x="count", showCategory=25)
    ggsave(pdf.name, egoBP.dotplot, device = "pdf", width = 30, height = 20, units = "cm")  
    
  }
  if(!is.null(egoCC)){
    csv.name <- paste(i,"_CC_dotplot.csv",sep="")
    pdf.name <- paste(i,"_CC_dotplot.pdf",sep="")
    write.csv(egoCC@result, file=csv.name, row.names=FALSE)
    egoCC.dotplot <- dotplot(egoCC, x="count", showCategory=25)
    ggsave(pdf.name, egoCC.dotplot, device = "pdf", width = 30, height = 20, units = "cm")  
  }
  if(!is.null(egoMF)){
    csv.name <- paste(i,"_MF_dotplot.csv",sep="")
    pdf.name <- paste(i,"_MF_dotplot.pdf",sep="")
    write.csv(egoMF@result, file=csv.name, row.names=FALSE)
    egoMF.dotplot <- dotplot(egoMF, x="count", showCategory=25)
    ggsave(paste(i,"_MF_dotplot.pdf",sep=""), egoMF.dotplot, device = "pdf", width = 30, height = 20, units = "cm")  
  }
}

##########---------- Configue
### Dupr
in.dir <- "/Volumes/Yolanda/Exp337CD25KONascent/3_DE-seq/0.1_original_GN/dupr"
wk.dir <- "/Volumes/Yolanda/Exp337CD25KONascent/3_DE-seq/2_Pathway/dupr"

##########---------- Main
in.file.vec <- list.files(path=in.dir, pattern="*.csv")
setwd(wk.dir)

for (i in in.file.vec) {
  if (grepl("KO", i) &&   grepl("WT", i)) {
    print("WT vs KO")
    padj.cutoff <- 0.05
    pval.cutoff <- 1
    log2fc.cutoff <- 0
  } else {
    padj.cutoff <- 0.01
    pval.cutoff <- 1
    log2fc.cutoff <- 3
  }
  
  i.simp.name <- gsub("_addGN.csv", "", i)
  print(paste("Start analysis of:", i, sep=" "))
  
  in.df <- read_csv(paste(in.dir, i, sep="/")) %>% filter(padj <= padj.cutoff) %>% filter(pvalue <= pval.cutoff)
  
  ###----- Up in first group
  in.df.up <- in.df %>% filter(log2FoldChange >= log2fc.cutoff)
  genes.i.up <- in.df.up$gene_name
  out.name.base <- paste(i.simp.name, "_up", sep="")
  GO_run(genes.i.up, out.name.base)
  
  ###----- Dn in first group
  in.df.dn <- in.df %>% filter(log2FoldChange <= -log2fc.cutoff)
  genes.i.dn <- in.df.dn$gene_name
  out.name.base <- paste(i.simp.name, "_dn", sep="")
  GO_run(genes.i.dn, out.name.base)
}


