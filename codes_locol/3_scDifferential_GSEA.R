######################################## GSEA analysis of single cell paga cluster differential genes ########################################
# Author: Huitian (Yolanda) Diao
# Nov 12th, 2019

######################################## Libraries ########################################
library("org.Mm.eg.db")
library(DOSE)
require(clusterProfiler)
library(dplyr)
library(enrichplot)
library(tidyverse)

######################################## Self-defined functions ########################################
GSEA_analysis_deseq2 <- function(de_df, out_name, gs_file, log2fc_cutoff, pval_cutoff) {
  #d3e_file <- d3e.files[1]
  #gs_file <- "/Volumes/Yolanda/Exp334CD25KOSc/source/GSEA/all_GSEA.csv"
  #log2fc_cutoff <- 0.05
  #pval_cutoff <- 1
  
  #####---------- Sample gsea analysis with custome gene list
  d3e.df <- de_df
  d3e.df <- d3e.df %>% filter(pvalue <= pval_cutoff )
  gene.list <- d3e.df$log2FoldChange
  names(gene.list) <- as.character(d3e.df$X1)
  gene.list <- sort(gene.list, decreasing = TRUE)
  deg.list <- names(gene.list)[abs(gene.list) > log2fc_cutoff]
  
  file.name.simp <- out_name
  
  print(file.name.simp)
  print(length(deg.list))
  
  #####---------- Read GSEA reference dataset
  gs.tb <- read_csv(gs_file)
  unique( gs.tb$gs_name)
  
  gs.name.simp <- unlist(strsplit(gs_file, "/"))
  gs.name.simp <- tail(gs.name.simp, 1)
  gs.name.simp <- gsub(".csv", "", gs.name.simp)
  
  #####---------- RUN GSEA
  em <- enricher(deg.list, TERM2GENE=gs.tb)
  em2 <- GSEA(gene.list, TERM2GENE=gs.tb,  nPerm = 10000, minGSSize = 1, maxGSSize = 5000,  pvalueCutoff = 1, by="DOSE")
  
  #####---------- Export results
  tb.name <- paste(file.name.simp, "---", gs.name.simp, ".csv", sep="")
  results.tb <- em2@result
  write_csv(results.tb, tb.name)
  gs <- em2@result$ID
  for (i in c(1:length(gs))){
    gsplot.i.name <- paste(file.name.simp, "---", gs.name.simp,"___",  gs[i], ".pdf", sep="")
    gsplot.i <- gseaplot2(em2, geneSetID = i, title = gs[i])
    ggsave(gsplot.i.name, gsplot.i, device="pdf", width=15, height=10, units="cm")
  }
  
}


######################################## Main ########################################

  
###----- Use nondupr
gs.file <- "/Volumes/Yolanda/Exp337CD25KONascent/2_DE-seq/3_GSEA_source/gsea_WT-vs-CD25KO_nondupr.csv"

###----- With DEseq2
if (FALSE) {
  wk.dir <- "/Volumes/Yolanda/Exp337CD25KONascent/3_scGSEA/nondupr_scDEseq2"
  in.dir <- "/Volumes/Yolanda/Exp334CD25KOSc/4_DEseq2/0_ClusterComparison_WTonly/0_DEseqOutput"
  setwd(wk.dir)
  
  cp <- "p0_vs_p1"
  cp.file <- "C0_vs_C3.csv"
  cp.tb <- read_csv(paste(in.dir, cp.file, sep="/"))
  GSEA_analysis_deseq2(cp.tb, cp, gs.file, 0, 0.05)
  
  cp <- "p2_vs_p1"
  cp.file <- "C9_vs_C3.csv"
  cp.tb <- read_csv(paste(in.dir, cp.file, sep="/"))
  GSEA_analysis_deseq2(cp.tb, cp, gs.file, 0, 0.05)
  
  cp <- "p2_vs_p3"
  cp.file <- "C9_vs_C8.csv"
  cp.tb <- read_csv(paste(in.dir, cp.file, sep="/"))
  GSEA_analysis_deseq2(cp.tb, cp, gs.file, 0, 0.05)
  
  cp <- "p2_vs_p4"
  cp.file <- "C9_vs_C4.csv"
  cp.tb <- read_csv(paste(in.dir, cp.file, sep="/"))
  GSEA_analysis_deseq2(cp.tb, cp, gs.file, 0, 0.05)
}

###----- With D3E