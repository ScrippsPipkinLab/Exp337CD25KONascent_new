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

##########---------- Main
in.dir <- "/Volumes/EXP337/Exp337CD25KONascent/3_DE-seq/DEseq2_out/1.2.1_pval_GN"
wk.dir <- "/Volumes/EXP337/Exp337CD25KONascent/3_DE-seq/DEseq2_out/1.2.1.1_pval_GN_pathway"

in.file.vec <- list.files(path=in.dir, pattern="*.csv")
in.file.vec <- tail(in.file.vec, n=6)
in.file.vec

### Go-term dot plot
if (FALSE) {
  for (i in in.file.vec) {
    #i <- in.file.vec[1]
    print(paste("Start analysis of:", i, sep=" "))
    setwd(wk.dir)
    in.df <- read.csv(paste(in.dir, i, sep="/"))
    genes.i <- as.character(unlist(in.df$gene_name))
    genes.i <- unique(genes.i)
    genes.i.id <- select(org.Mm.eg.db, genes.i, c("ENTREZID"), "ALIAS")
    #genes.i.id$ENTREZID

    egoBP <- enrichGO(gene = genes.i.id$ENTREZID, keyType = 'ENTREZID', OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "none", pvalueCutoff = 0.05, readable = TRUE)
    egoCC <- enrichGO(gene = genes.i.id$ENTREZID, keyType = 'ENTREZID', OrgDb = org.Mm.eg.db, ont = "CC", pAdjustMethod = "none", pvalueCutoff = 0.05, readable = TRUE)
    egoMF <- enrichGO(gene = genes.i.id$ENTREZID, keyType = 'ENTREZID', OrgDb = org.Mm.eg.db, ont = "MF", pAdjustMethod = "none", pvalueCutoff = 0.05, readable = TRUE)
    
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
}


