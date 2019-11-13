######################################## GSEA - make table ########################################
# Author: Huitian (Yolanda) Diao
# July 22nd, 2019
# Create GSEA tables

######################################## Libraries ########################################
library(dplyr)
library(tidyverse)

######################################## Self-defined function ########################################

wk.dir <- "/Volumes/Yolanda/Exp337CD25KONascent/3_DE-seq/3_GSEA"
setwd(wk.dir)

###----- Use non-dupr...
in.base <- "/Volumes/Yolanda/Exp337CD25KONascent/3_DE-seq/0.1_original_GN/nondupr/nondupr_"
wt.ko.out.name <- "gsea_WT-vs-CD25KO_nondupr.csv"

### WT_vs_KO
log2fc.cutoff <- 1
padj.cutoff <- 0.05
pval.cutoff <- 1

if (TRUE) {
  WT_vs_KO_0h <- paste(in.base, "WT_0h_vs_KO_0h_addGN.csv", sep="")
  WT_vs_KO_6h <- paste(in.base, "WT_6h_vs_KO_6h_addGN.csv", sep="") 
  WT_vs_KO_24h <- paste(in.base, "WT_24h_vs_KO_24h_addGN.csv", sep="") 
  WT_vs_KO_48h <- paste(in.base, "WT_48h_vs_KO_48h_addGN.csv", sep="") 
  
  tb.0h <- read_csv(WT_vs_KO_0h) %>% filter(abs(log2FoldChange) >= log2fc.cutoff) %>% filter(pvalue <= pval.cutoff) %>% filter(padj <= padj.cutoff)
  tb.6h <- read_csv(WT_vs_KO_6h) %>% filter(abs(log2FoldChange) >= log2fc.cutoff) %>% filter(pvalue <= pval.cutoff) %>% filter(padj <= padj.cutoff)
  tb.24h <- read_csv(WT_vs_KO_24h) %>% filter(abs(log2FoldChange) >= log2fc.cutoff) %>% filter(pvalue <= pval.cutoff) %>% filter(padj <= padj.cutoff)
  tb.48h <- read_csv(WT_vs_KO_48h) %>% filter(abs(log2FoldChange) >= log2fc.cutoff) %>% filter(pvalue <= pval.cutoff) %>% filter(padj <= padj.cutoff)
  
  wt.up.0h <- tb.0h %>% filter(log2FoldChange > 0) %>% .$gene_name %>% unique()
  wt.dn.0h <- tb.0h %>% filter(log2FoldChange < 0) %>% .$gene_name %>% unique()
  wt.up.6h <- tb.6h %>% filter(log2FoldChange > 0) %>% .$gene_name %>% unique()
  wt.dn.6h <- tb.6h %>% filter(log2FoldChange < 0) %>% .$gene_name %>% unique()
  wt.up.24h <- tb.24h %>% filter(log2FoldChange > 0) %>% .$gene_name %>% unique()
  wt.dn.24h <- tb.24h %>% filter(log2FoldChange < 0) %>% .$gene_name %>% unique()
  wt.up.48h <- tb.48h %>% filter(log2FoldChange > 0) %>% .$gene_name %>% unique()
  wt.dn.48h <- tb.48h %>% filter(log2FoldChange < 0) %>% .$gene_name %>% unique()
  
  gs.names <- c(rep("WTvsKO_0h_up", length(wt.up.0h)),
                rep("WTvsKO_0h_dn", length(wt.dn.0h)),
                rep("WTvsKO_6h_up", length(wt.up.6h)),
                rep("WTvsKO_6h_dn", length(wt.dn.6h)),
                rep("WTvsKO_24h_up", length(wt.up.24h)),
                rep("WTvsKO_24h_dn", length(wt.dn.24h)),
                rep("WTvsKO_48h_up", length(wt.up.48h)),
                rep("WTvsKO_48h_dn", length(wt.dn.48h)))
  gene.symbols <- c(wt.up.0h, wt.dn.0h, wt.up.6h, wt.dn.6h, wt.up.24h, wt.dn.24h, wt.up.48h, wt.dn.48h)
  
  wt.vs.ko.gsea.tb <- tibble(gs_name = gs.names, gene_symbol = gene.symbols)
  
  write_csv(wt.vs.ko.gsea.tb, wt.ko.out.name)
}

