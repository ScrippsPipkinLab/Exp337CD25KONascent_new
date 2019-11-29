#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 15:18:02 2019

@author: yolandatiao
"""

########## Create venn diagrams for DEseq output ##########
# Author: Huitian (Yolanda) Diao
# Nov 12, 2019
# Add reference size for plotting comparable size circles

########## Import ##########
import os # For changing directory
import pandas as pd # For using pandas in seaborn plot
from matplotlib_venn import venn3, venn3_circles, venn2
import matplotlib.pyplot as plt

########## Self defined functions ##########
def venn2_de_ref1000(in_dict, out_name, pval_cutoff, log2fc_cutoff):
    #padj_cutoff = 0.05
    #log2fc_cutoff = 1
    #in_dict = {"wt6h_0h":wt6h_0h, "ko6h_0h":ko6h_0h}
    #out_name = "WT6h-0h__vs__KO6h-0h.pdf"
    
    f1_name, f2_name = list(in_dict.keys())
    f1, f2 = in_dict[f1_name], in_dict[f2_name]
    
    tb1 = pd.read_csv(f1)
    tb1_sig = tb1[tb1['pvalue'] <= pval_cutoff]
    tb1_sig_up = tb1_sig[tb1_sig['log2FoldChange'] >= log2fc_cutoff]
    tb1_sig_dn = tb1_sig[tb1_sig['log2FoldChange'] <= - log2fc_cutoff]
    
    tb2 = pd.read_csv(f2)
    tb2_sig = tb2[tb2['pvalue'] <= pval_cutoff]
    tb2_sig_up = tb2_sig[tb2_sig['log2FoldChange'] >= log2fc_cutoff]
    tb2_sig_dn = tb2_sig[tb2_sig['log2FoldChange'] <= - log2fc_cutoff]
    
    up1 = set(list(tb1_sig_up.iloc[:, 0]))
    dn1 = set(list(tb1_sig_dn.iloc[:, 0]))
    up2 = set(list(tb2_sig_up.iloc[:, 0]))
    dn2 = set(list(tb2_sig_dn.iloc[:, 0]))
    
    set_ref = set(["n%s"%str(i) for i in range(0,1000)])
    
    plt.figure()
    up_name = out_name + "_up.pdf"
    v3c = venn3([up1, up2, set_ref], (f1_name, f2_name, "ref1000"))
    plt.savefig(up_name, transparent=True)
    
    plt.figure()
    dn_name = out_name + "_dn.pdf"
    v3c = venn3([dn1, dn2, set_ref], (f1_name, f2_name, "ref1000"))
    plt.savefig(dn_name, transparent=True)
    


########## Main ##########
wk_dir = "/Volumes/Yolanda/Exp337CD25KONascent/2_DE-seq/4_DE_Venn"
os.chdir(wk_dir)
in_base = "/Volumes/Yolanda/Exp337CD25KONascent/2_DE-seq/0.1_original_GN/nondupr/nondupr_"




ko6h_0h = in_base + "KO_6h_vs_KO_0h_addGN.csv"
ko24h_0h = in_base + "KO_24h_vs_KO_0h_addGN.csv"
ko48h_0h = in_base + "KO_48h_vs_KO_0h_addGN.csv"
wt6h_0h = in_base + "WT_6h_vs_WT_0h_addGN.csv"
wt24h_0h = in_base + "WT_24h_vs_WT_0h_addGN.csv"
wt48h_0h = in_base + "WT_48h_vs_WT_0h_addGN.csv"
wt_ko_6h = in_base + "WT_6h_vs_KO_6h_addGN.csv"
wt_ko_24h = in_base + "WT_24h_vs_KO_24h_addGN.csv"
wt_ko_48h = in_base + "WT_48h_vs_KO_48h_addGN.csv"


use_dict = {"wt6h_0h":wt6h_0h, "ko6h_0h":ko6h_0h}
plot_name = "WT6h-0h__vs__KO6h-0h"
venn2_de_ref1000(use_dict, plot_name, 0.05, 1)

use_dict = {"wt24h_0h":wt24h_0h, "ko24h_0h":ko24h_0h}
plot_name = "WT24h-0h__vs__KO24h-0h"
venn2_de_ref1000(use_dict, plot_name, 0.05, 1)

use_dict = {"wt48h_0h":wt48h_0h, "ko48h_0h":ko48h_0h}
plot_name = "WT48h-0h__vs__KO48h-0h"
venn2_de_ref1000(use_dict, plot_name, 0.05, 1)

###----- WT vs KO
use_dict = {"wt6h_0h":wt6h_0h, "wt6h_ko6h":wt_ko_6h}
plot_name = "WT6h-0h__vs__WT6h-KO6h"
venn2_de_ref1000(use_dict, plot_name, 0.05, 1)

use_dict = {"wt24h_0h":wt24h_0h, "wt24h_ko24h":wt_ko_24h}
plot_name = "WT24h-0h__vs__WT24h-KO24h"
venn2_de_ref1000(use_dict, plot_name, 0.05, 1)

use_dict = {"wt48h_0h":wt48h_0h, "wt48h_ko48h":wt_ko_48h}
plot_name = "WT48h-0h__vs__WT48h-KO48h"
venn2_de_ref1000(use_dict, plot_name, 0.05, 1)
