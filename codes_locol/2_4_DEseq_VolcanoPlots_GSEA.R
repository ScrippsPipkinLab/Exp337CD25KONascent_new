########## DEseq Volcano plots ##########

######################################## Libraries ########################################
library(dplyr)
library(tidyverse)
library(ggrepel)
library(grid)
library(DOSE)
library(enrichplot)
require(clusterProfiler)


######################################## Imports ########################################
###----- In this function:
### used p value for plotting & cutoffs
vol_plot <- function(in_file, mlog10p_cutoff, log2fc_cutoff, 
                     use_genes,log2fc_lim, mlog10p_lim, 
                     compare_pair, wid, hei, 
                     virtual_p_add, use_density) {
  # Test parameters
  if (FALSE) {
    mlog10p_cutoff <- 1
    log2fc_cutoff <- 1
    in_file <- "/Volumes/Yolanda/Exp337CD25KONascent/2_DE-seq/0.1_original_GN/nondupr/nondupr_WT_0h_vs_KO_0h_addGN.csv"
    log2fc_lim <- FALSE
    mlog10p_lim <- FALSE
    use_genes <- c("Tbx21", "Prdm1", "Id2", "Id3", "Slamf6", "Tcf7", "Il2ra", "Sell", "Cd44", "Tox", "Klrg1")
    compare_pair <- c("KO_24h", "KO_0h")
    wid <- 4.5
    hei <- 4.5
    virtual_p_add <- TRUE 
    use_density <- FALSE
  }
  
  use_col_vec <- c( "firebrick", "cyan4")
  ###----- Read input file and tranform data
  if (TRUE) {
    in_file_nopath <- basename(in_file)
    out_file <- gsub(".csv", "_vplot", in_file_nopath)
    gsea_file <- gsub(".csv", "_GSEA", in_file_nopath)
    out_file_nolegend <- gsub(".csv", "_vplot_nolegend", in_file_nopath)
    out_file_noAnno <- gsub(".csv", "_vplot_noAnno", in_file_nopath)
    in_tab <- read_csv(in_file)
    in_tab <- in_tab %>% mutate(padj = replace_na(padj, 1)) %>% mutate(pvalue = replace_na(pvalue, 1))
    colnames(in_tab) <- c("transcript_name", colnames(in_tab)[2:length(colnames(in_tab))])
    
    if (virtual_p_add) {
      in_tab <- in_tab %>% mutate(mlog10p = - log10(pvalue + 0.00000000001)) %>% mutate(sqrt_mlog10p = sqrt(mlog10p))
    } else {
      in_tab <- in_tab %>% mutate(mlog10p = - log10(pvalue)) %>% mutate(sqrt_mlog10p = sqrt(mlog10p))
    }
    in_tab <- in_tab %>% mutate(sqrt_log2fc = (sqrt(abs(log2FoldChange)) * log2FoldChange / abs(log2FoldChange)))
  }
  
  ###----- Summarize the upper and lower limits of input
  if (TRUE) {
    max_mlog10p <- max(in_tab$mlog10p, na.rm = TRUE)
    min_mlog10p <- min(in_tab$mlog10p, na.rm = TRUE)
    max_log2FC <- max(in_tab$log2FoldChange, na.rm = TRUE)
    min_log2FC <- min(in_tab$log2FoldChange, na.rm = TRUE)
    max_sqrt_mlog10p <- max(in_tab$sqrt_mlog10p, na.rm = TRUE)
    min_sqrt_mlog10p <- min(in_tab$sqrt_mlog10p, na.rm = TRUE)
    lims_vec <- c(max_mlog10p, min_mlog10p,
                  max_log2FC, min_log2FC,
                  max_sqrt_mlog10p, min_sqrt_mlog10p)
    lims_names <- c("max_mlog10p", "min_mlog10p", 
                    "max_log2FC", "min_log2FC", 
                    "max_sqrt_mlog10p", "min_sqrt_mlog10p")
    lims_tb <- tibble(limits=lims_names, vals=lims_vec)
    colnames(lims_tb) <- c("limits", gsub(".csv", "", in_file_nopath))
  }
  
  ###----- Determine limits of graph
  if (log2fc_lim==FALSE) {
    log2fc_lim <- ceiling(max(max_log2FC, abs(min_log2FC)))
  }
  if (mlog10p_lim==FALSE){
    mlog10p_lim <- ceiling(max_sqrt_mlog10p)
  }
  
  ###----- Separate into non-significant & up / down reaulgated groups
  if (TRUE) {
    cp1_name <- paste(compare_pair[1], "up", sep="_")
    cp2_name <- paste(compare_pair[2], "up", sep="_")
    in_tb_cp1 <- in_tab %>% filter(mlog10p >= mlog10p_cutoff) %>% filter(log2FoldChange >= log2fc_cutoff)
    in_tb_cp2 <- in_tab %>% filter(mlog10p >= mlog10p_cutoff) %>% filter(log2FoldChange <= -log2fc_cutoff)
    in_tb_cp1$group <- rep(cp1_name, nrow(in_tb_cp1))
    in_tb_cp2$group <- rep(cp2_name, nrow(in_tb_cp2))
    in_tb_cp1 <- in_tb_cp1 %>% dplyr::select(c("transcript_name", "group"))
    in_tb_cp2 <- in_tb_cp2 %>% dplyr::select(c("transcript_name", "group"))
    in_tb_cp <- bind_rows(in_tb_cp1, in_tb_cp2)
    in_tab <- in_tab %>% left_join(in_tb_cp, by="transcript_name")
    in_tab <- in_tab %>% mutate(group = replace_na(group, "n.s."))
    in_tb_ns <- in_tab %>% filter(group == "n.s.")
    in_tb_cp1 <- in_tab %>% filter(group == cp1_name)
    in_tb_cp2 <- in_tab %>% filter(group == cp2_name)
    cp1_nu <- nrow(in_tb_cp1)
    cp2_nu <- nrow(in_tb_cp2)
    in_tb_cp <- bind_rows(in_tb_cp1, in_tb_cp2)
    text_annotations_1 <- data.frame(xpos = c(Inf), ypos = c(Inf), 
                                     annotate_Text = c(as.character(cp1_nu)),
                                     hjustvar = c(1), vjustvar = c(1)) #<- adjust
    text_annotations_2 <- data.frame(xpos = c(-Inf), ypos =  c(Inf),
                                     annotate_Text = c(as.character(cp2_nu)),
                                     hjustvar = c(0), vjustvar = c(1)) #<- adjust
    grob1 <- grobTree(textGrob(as.character(cp1_nu), x=0.8,  y=0.95, hjust=0,
                              gp=gpar(col=use_col_vec[1], fontsize=10)))
    grob2 <- grobTree(textGrob(as.character(cp2_nu), x=0.05,  y=0.95, hjust=0,
                               gp=gpar(col=use_col_vec[2], fontsize=10)))
    print(paste(compare_pair[1], " selected gene number: ", as.character(cp1_nu)), sep="")
    print(paste(compare_pair[2], " selected gene number: ", as.character(cp2_nu)), sep="")
    #anno_tab <- in_tab %>% filter(mlog10p > mlog10p_cutoff) %>% filter(abs(log2FoldChange) > log2fc_cutoff)
    #anno_tab2 <- in_tab %>% filter(mlog10p > 4)
    #anno_tab <- anno_tab %>% bind_rows(anno_tab2) %>% distinct()
    anno_tab <- in_tab %>% dplyr::select(one_of(c("gene_name", "log2FoldChange", "pvalue", "mlog10p", "sqrt_mlog10p", "sqrt_log2fc", "group"))) %>% na.omit() %>%
      filter(group != "n.s.") %>% filter(gene_name %in% use_genes) 
    
    if (nrow(anno_tab) > 0) {
      anno_tab_KO <- anno_tab %>% filter(log2FoldChange > 0)
      anno_tab_KO$use_color <- rep(use_col_vec[2], nrow(anno_tab_KO))
      anno_tab_WT <- anno_tab %>% filter(log2FoldChange < 0)
      anno_tab_WT$use_color <- rep(use_col_vec[1], nrow(anno_tab_WT))
      anno_tab <- anno_tab_WT %>% bind_rows(anno_tab_KO)
    }
    
  }

  ###----- Density plot
  if (TRUE) {
    if (use_density) {
      volcano_plot <- ggplot(in_tab, aes(x=sqrt_log2fc, y=mlog10p)) +
        geom_vline(xintercept=0, size=0.3, alpha=1) +
        geom_vline(xintercept=1, size=0.1, alpha=0.7) +
        geom_vline(xintercept=-1, size=0.1, alpha=0.7) +
        geom_hline(yintercept=1, size=0.1, alpha=0.7) +
        geom_vline(xintercept=-sqrt(log2fc_cutoff), size=0.1, alpha=0.5) +
        geom_vline(xintercept=sqrt(log2fc_cutoff), size=0.1, alpha=0.5) +
        geom_hline(yintercept=mlog10p_cutoff, size=0.1, alpha=0.7) +
        stat_density_2d(aes(alpha= sqrt(..density..)), geom = "raster", contour = FALSE)
    } else {
      volcano_plot <- ggplot(in_tab, aes(x=sqrt_log2fc, y=mlog10p)) +
        geom_vline(xintercept=0, size=0.3, alpha=1) +
        geom_vline(xintercept=1, size=0.1, alpha=0.7) +
        geom_vline(xintercept=-1, size=0.1, alpha=0.7) +
        geom_hline(yintercept=1, size=0.1, alpha=0.7) +
        geom_vline(xintercept=-sqrt(log2fc_cutoff), size=0.1, alpha=0.5) +
        geom_vline(xintercept=sqrt(log2fc_cutoff), size=0.1, alpha=0.5) +
        geom_hline(yintercept=mlog10p_cutoff, size=0.1, alpha=0.7)
    }
  }
  
  ###----- Volcano plot
  if (TRUE) {
    volcano_plot_basic <- volcano_plot +
      geom_point(data=in_tb_ns, aes(x=sqrt_log2fc, y=mlog10p), stroke=0, alpha=0.5, size=0.3, color="gray75") +
      geom_point(data=in_tb_cp1, aes(x=sqrt_log2fc, y=mlog10p), stroke=0, alpha=0.5, size=0.3, color=use_col_vec[1]) +
      geom_point(data=in_tb_cp2, aes(x=sqrt_log2fc, y=mlog10p), stroke=0, alpha=0.5, size=0.3, color=use_col_vec[2]) +
      scale_x_continuous(limits = c(-log2fc_lim, log2fc_lim), expand = c(0, 0)) +
      scale_y_sqrt(limits = c(0, mlog10p_lim), expand = c(0, 0)) +
      annotation_custom(grob1) + annotation_custom(grob2) +
      #geom_text(data = text_annotations_1, aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotate_Text, color=use_col_vec[1])) +
      #geom_text(data = text_annotations_2, aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotate_Text, color=use_col_vec[2])) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank())
    #volcano_plot
    if (nrow(anno_tab) > 0) {
      volcano_plot <- volcano_plot_basic +     
        geom_point(data = anno_tab, aes(x=sqrt_log2fc, y=mlog10p, color=anno_tab$use_color, stroke=0),  size=0.5) +
        scale_color_manual(values = use_col_vec) +
        geom_text_repel(data = anno_tab, label=anno_tab$gene_name, max.iter = 40000, force = 1, aes(color=anno_tab$use_color), size=3, segment.alpha = 0.3, segment.size = 0.5)
    } else {
      print("No annotation")
      volcano_plot <- volcano_plot_basic
    }
    #volcano_plot
    
    volcano_plot_nolegend <- volcano_plot +
      theme(legend.position="none",
            axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.y=element_blank())
    
    volcano_plot_noAnno <- volcano_plot_basic +
      theme(legend.position="none",
            axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.y=element_blank())
  }
  
  ###----- Save
  ggsave(paste(out_file, ".pdf", sep=""), volcano_plot, device="pdf", width=(wid+6), height=hei, units=('cm'), dpi=300)
  ggsave(paste(out_file_nolegend, ".pdf", sep=""), volcano_plot_nolegend, device="pdf", width=wid, height=hei-1, units=('cm'), dpi=300)
  ggsave(paste(out_file_noAnno, ".pdf", sep=""), volcano_plot_noAnno, device="pdf", width=wid, height=hei-1, units=('cm'), dpi=300)
  ggsave(paste(out_file, ".png", sep=""), volcano_plot, device="png", width=(wid+6), height=hei, units=('cm'), dpi=300)
  ggsave(paste(out_file_nolegend, ".png", sep=""), volcano_plot_nolegend, device="png", width=wid, height=hei-1, units=('cm'), dpi=300)
  ggsave(paste(out_file_noAnno, ".png", sep=""), volcano_plot_noAnno, device="png", width=wid, height=hei-1, units=('cm'), dpi=300)
  
  ###----- GSEA analysis
  if (FALSE) {
    in_tb_sig <- in_tab %>% filter(padj <= 0.05)
    gsea.gene.list <- in_tb_sig$log2FoldChange
    names(gsea.gene.list) <- as.character(in_tb_sig$gene_name)
    gsea.gene.list <- sort(gsea.gene.list, decreasing = TRUE)
    
    print(paste("Gene number for GSEA analysis: ", as.character(length(gsea.gene.list)), sep=""))
    GSEA_analysis(gsea.gene.list, gsea_file)
  }
  
  return(lims_tb)
}

###----- In GSEA analysis:
### used padj for setting gene cutoff
GSEA_analysis <- function(gene_list, out_name, gs_file) {
  #gene_list <- gsea.gene.list
  
  #####---------- Sample gsea analysis with custome gene list
  gene.list <- gene_list
  file.name.simp <- out_name
  
  #####---------- Read GSEA reference dataset
  gs.tb <- read_csv(gs_file)
  unique( gs.tb$gs_name)
  
  gs.name.simp <- unlist(strsplit(gs_file, "/"))
  gs.name.simp <- tail(gs.name.simp, 1)
  gs.name.simp <- gsub(".csv", "", gs.name.simp)
  
  #####---------- RUN GSEA
  em2 <- GSEA(gene.list, TERM2GENE=gs.tb,  nPerm = 10000, minGSSize = 1, maxGSSize = 5000,  pvalueCutoff = 1, by="DOSE")
  
  #####---------- Export results
  tb.name <- paste(file.name.simp, "---", gs.name.simp, ".csv", sep="")
  results.tb <- em2@result
  write_csv(results.tb, tb.name)
  gs <- em2@result$ID
  for (i in c(1:length(gs))){
    gsplot.i.name.pdf <- paste(file.name.simp, "---", gs.name.simp,"___",  gs[i], ".pdf", sep="")
    gsplot.i.name.png <- paste(file.name.simp, "---", gs.name.simp,"___",  gs[i], ".png", sep="")
    gsplot.i <- gseaplot2(em2, geneSetID = i, title = gs[i])
    ggsave(gsplot.i.name.pdf, gsplot.i, device="pdf", width=15, height=10, units="cm")
    ggsave(gsplot.i.name.png, gsplot.i, device="png", width=15, height=10, units="cm")
  }
  
}

###----- Designed for naming in this experiment only...
simp_name_cvt <- function(in_name) {
  #in_name <- files[1]
  out_name <- strsplit(basename(in_name), split="\\.")[[1]][1]
  out_name <- strsplit(out_name, split="_GSEA")[[1]][1]
  out_name <- gsub("nondupr_", "", out_name)
  out_name <- gsub("dupr_", "", out_name)
  out_name <- gsub("_addGN", "", out_name)
  return(out_name)
}

le_to_le_sig <- function(vec_x) {
  out_vec <- c()
  for (i in vec_x){
    #i <- "tags=58%, list=12%, signal=52%"
    i_vec <- unlist(strsplit(i, ", "))
    i_sig <- tail(i_vec, n=1)
    i_sig <- gsub("signal=", "", i_sig)
    i_sig <- gsub("%", "", i_sig)
    i_sig <- as.numeric(i_sig)
    out_vec <- c(out_vec, i_sig)
  }
  return(out_vec)
}

GSEA_sum_new <- function(file_list, outname, wid, hei, input_order_vec, path_order_vec) {
  #file_list <- gsea.files.use
  #outname <- "louvain_expr_GSEA"
  #wid <- 30
  #hei <- 20
  #input_order_vec <- order_comparisons
  #comp_new_name <- order_cp_newnames
  #path_order_vec <- c("NAV_up", "NAV_dn", "Act48h_up", "Act48h_dn", 
  #                    "EE_up", "EE_dn", "MP_up", "MP_dn", "DP_up", "DP_dn", "TE_up", "TE_dn",
  #                    "MEM_up", "MEM_dn", "EXH_up", "EXH_dn")
  
  #file_list <- gsea.files.use
  cp.vec <- c()
  pw.vec <- c()
  nes.vec <- c()
  padj.vec <- c()
  le.sig.vec <- c()
  for (filex in file_list) {
    #filex <- gsea.files.use[1]
    cpx <- simp_name_cvt(filex)
    filex.tb <- read_csv(filex)
    cp.vec <- c(cp.vec, rep(cpx, nrow(filex.tb)))
    pw.vec <- c(pw.vec, filex.tb$ID)
    nes.vec <- c(nes.vec, filex.tb$NES)
    padj.vec <- c(padj.vec, filex.tb$p.adjust)
    le.x <- filex.tb$leading_edge
    le.sig.x <- le_to_le_sig(le.x)
    le.sig.vec <- c(le.sig.vec, le.sig.x)
  }
  
  plot.tb <- tibble(comparison=cp.vec,
                    pathway=pw.vec,
                    NES=nes.vec,
                    padj=padj.vec,
                    leadingEdge_signal=le.sig.vec)
  plot.tb$mlog10padj <- -log10(plot.tb$padj)
  plot.tb$comparison <- factor(plot.tb$comparison, levels=input_order_vec)
  new.factor <- factor(plot.tb$pathway, levels=path_order_vec)
  plot.tb$pathway <- new.factor
  
  
  
  bbplot <- ggplot(plot.tb, aes(pathway, comparison)) +
    geom_point(aes(size=mlog10padj, color=NES)) +
    scale_color_gradient2(low="blue", mid="white", high="red", midpoint=0) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  outname_csv <- paste(outname, ".csv", sep="")
  outname_pdf <- paste(outname, ".pdf", sep="")
  write_csv(plot.tb, outname_csv)
  ggsave(outname_pdf, bbplot, device="pdf", width=wid, height=hei, units="cm")
}

######################################## Main ########################################

###----- Volcano plots
if (TRUE) {
  wk.dir <- "/Volumes/Yolanda/Exp337CD25KONascent/2_DE-seq/4_DE_volcanoPlots_GSEA/nondupr/0_vplots"
  setwd(wk.dir)
  
  in.dir <- "/Volumes/Yolanda/Exp337CD25KONascent/2_DE-seq/0.1_original_GN/nondupr"
  files <- list.files(path=in.dir, pattern="*.csv", full.names=TRUE, recursive=FALSE)
  genes.use <- c("Tbx21", "Prdm1", "Id2", "Id3", "Slamf6", "Tcf7", "Il2ra", "Sell", "Cd44", "Tox", "Klrg1")
  
  lims.list <- list()
  for (i in c(1:length(files))){
    #i <- 1
    file.i <- files[i]
    i.name <- basename(file.i)
    i.name <- gsub("nondupr_", "", i.name)
    i.name <- gsub("dupr_", "", i.name)
    i.name <- gsub("_addGN.csv", "", i.name)
    comp.pair <- strsplit(i.name, split="_vs_")
    comp.pair <- unlist(comp.pair)
    
    
    # in_file, mlog10p_cutoff, log2fc_cutoff
    # use_genes,log2fc_lim, mlog10p_lim
    # compare_pair, wid, hei, 
    # virtual_p_add, use_density
    if (grepl("KO", file.i) &&   grepl("WT", file.i)) {
      lims.i <- vol_plot(file.i, 1.30103, 1, 
                         genes.use, 5, 15, 
                         comp.pair, 4.5, 4.5,
                         TRUE, FALSE)
    }  else {
      lims.i <- vol_plot(file.i, 2, 2, 
                         genes.use, 5, 15, 
                         comp.pair, 4.5, 4.5,
                         TRUE, FALSE)
      }
    lims.list[[i]] <- lims.i
  }
  
  lims.tb <- lims.list[[1]]
  for (i in c(2:length(files))) {
    lims.tb <- lims.tb %>% left_join(lims.list[[i]])
  }
  
  write_csv(lims.tb, "lims.csv")
}

###----- GSEA only
if (FALSE) {
  wk.dir <- "/Volumes/Yolanda/Exp337CD25KONascent/2_DE-seq/4_DE_volcanoPlots_GSEA/nondupr"
  setwd(wk.dir)
  
  in.dir <- "/Volumes/Yolanda/Exp337CD25KONascent/2_DE-seq/0.1_original_GN/nondupr"
  files <- list.files(path=in.dir, pattern="*.csv", full.names=TRUE, recursive=FALSE)
  
  gs_file1 <- "/Volumes/Yolanda/Exp334CD25KOSc/source/GSEA/Best_JJM_GSEA.csv"
  gs_file2 <- "/Volumes/Yolanda/Exp334CD25KOSc/source/GSEA/GSE88987-Exp276_slt_GSEA.csv"
  gs_file3 <- "/Volumes/Yolanda/Exp334CD25KOSc/source/GSEA/Yao_Tox_sc_GSEA.csv"
  
  gs_file <- gs_file3
  
  for (i in files) {
    tb.i <- read_csv(i)
    i.base <- basename(i)
    gsea_file <- gsub(".csv", "_GSEA", i.base)
    
    in_tb_sig <- tb.i %>% filter(padj <= 0.05)
    gsea.gene.list <- in_tb_sig$log2FoldChange
    names(gsea.gene.list) <- as.character(in_tb_sig$gene_name)
    gsea.gene.list <- sort(gsea.gene.list, decreasing = TRUE)
    print(paste("Gene number for GSEA analysis: ", as.character(length(gsea.gene.list)), sep=""))
    GSEA_analysis(gsea.gene.list, gsea_file, gs_file)
  }
}

###----- GSEA summary: WT versus KO -- GSE88987-Exp276
if (FALSE) {
  in.dir <- "/Volumes/Yolanda/Exp337CD25KONascent/2_DE-seq/4_DE_volcanoPlots_GSEA/nondupr/1_GSEA_GSE88987-Exp276/WT_KO"
  wk.dir <- "/Volumes/Yolanda/Exp337CD25KONascent/2_DE-seq/4_DE_volcanoPlots_GSEA/nondupr"
  setwd(wk.dir)
  out.name <- "WTvsKO_GSE88987-Exp276_GSEA_sum"
  
  files <- list.files(path=in.dir, pattern="*.csv", full.names=TRUE, recursive=FALSE)
  # Check sample simplified names & all possible pathways
  # For setting order lists
  if (TRUE) {
    file.simp.vec <- c()
    path.all.vec <- c()
    for (i in files) {
      #i <- files[1]
      tb.i <- read_csv(i)
      path.all.vec <- tb.i$ID
      file.simp.vec <- c(file.simp.vec, simp_name_cvt(i))
    }
    path.all.vec <- unique(path.all.vec)
  }
  print(file.simp.vec)
  print(path.all.vec)
  
  file.ordered.vec <- rev(c("WT_0h_vs_KO_0h", "WT_6h_vs_KO_6h", "WT_24h_vs_KO_24h", "WT_48h_vs_KO_48h"))
  path.ordered.vec <- c("NAV_up", "NAV_dn","Act48h_up", "Act48h_dn", 
                        "TE_up", "TE_dn", "EE_up", "EE_dn", "DP_up", "DP_dn", "MP_up", "MP_dn", 
                        "MEM_up", "MEM_dn", "EXH_up", "EXH_dn")
  
  GSEA_sum_new(files, out.name, 20, 8, file.ordered.vec, path.ordered.vec)
  
}

###----- GSEA summary: WT versus KO -- Best
if (FALSE) {
  in.dir <- "/Volumes/Yolanda/Exp337CD25KONascent/2_DE-seq/4_DE_volcanoPlots_GSEA/nondupr/1_GSEA_Best_JJM/WT_KO"
  wk.dir <- "/Volumes/Yolanda/Exp337CD25KONascent/2_DE-seq/4_DE_volcanoPlots_GSEA/nondupr"
  setwd(wk.dir)
  out.name <- "WTvsKO_1_GSEA_Best_JJM_sum"
  
  files <- list.files(path=in.dir, pattern="*.csv", full.names=TRUE, recursive=FALSE)
  # Check sample simplified names & all possible pathways
  # For setting order lists
  if (TRUE) {
    file.simp.vec <- c()
    path.all.vec <- c()
    for (i in files) {
      #i <- files[1]
      tb.i <- read_csv(i)
      path.all.vec <- tb.i$ID
      file.simp.vec <- c(file.simp.vec, simp_name_cvt(i))
    }
    path.all.vec <- unique(path.all.vec)
  }
  print(file.simp.vec)
  print(path.all.vec)
  
  file.ordered.vec <- rev(c("WT_0h_vs_KO_0h", "WT_6h_vs_KO_6h", "WT_24h_vs_KO_24h", "WT_48h_vs_KO_48h"))
  path.ordered.vec <- c("best_cluster_1", "best_cluster_2", "best_cluster_3", "best_cluster_4", "best_cluster_5",
                        "best_cluster_6", "best_cluster_7", "best_cluster_8", "best_cluster_9", "best_cluster_10",
                        "JJM_CirculatingSignature", "JJM_TrmSignature")
  
  GSEA_sum_new(files, out.name, 20, 10, file.ordered.vec, path.ordered.vec)
  
}

###----- GSEA summary: WT versus KO -- Yao
if (TRUE) {
  in.dir <- "/Volumes/Yolanda/Exp337CD25KONascent/2_DE-seq/4_DE_volcanoPlots_GSEA/nondupr/1_GSEA_Yao/WT_KO"
  wk.dir <- "/Volumes/Yolanda/Exp337CD25KONascent/2_DE-seq/4_DE_volcanoPlots_GSEA/nondupr"
  setwd(wk.dir)
  out.name <- "WTvsKO_1_GSEA_Yao"
  
  files <- list.files(path=in.dir, pattern="*.csv", full.names=TRUE, recursive=FALSE)
  # Check sample simplified names & all possible pathways
  # For setting order lists
  if (TRUE) {
    file.simp.vec <- c()
    path.all.vec <- c()
    for (i in files) {
      #i <- files[1]
      tb.i <- read_csv(i)
      path.all.vec <- tb.i$ID
      file.simp.vec <- c(file.simp.vec, simp_name_cvt(i))
    }
    path.all.vec <- unique(path.all.vec)
  }
  print(file.simp.vec)
  print(path.all.vec)
  
  file.ordered.vec <- rev(c("WT_0h_vs_KO_0h", "WT_6h_vs_KO_6h", "WT_24h_vs_KO_24h", "WT_48h_vs_KO_48h"))
  path.ordered.vec <- c( "Yao_MP", "Yao_PRO", "Yao_EXH")
  
  GSEA_sum_new(files, out.name, 10, 8, file.ordered.vec, path.ordered.vec)
  
}