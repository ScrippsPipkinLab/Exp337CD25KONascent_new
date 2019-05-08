#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 16:30:39 2019

@author: yolandatiao
"""

########## ExonID to Genename ##########
# Author: Huitian (Yolanda) Diao
# Apr 30, 2019

########## Import ##########
import os
import csv
import glob
from astropy.io import ascii

########## Self defined functions ##########
def TIDtoGN(fileX, idList, gnList):    
    outName = fileX.replace(".csv", "_GN.csv")
    with open(fileX, "r") as fin:
        with open(outName, "w") as fout:
            rfin = csv.reader(fin, delimiter=",")
            wfout = csv.writer(fout, delimiter=",")
            header = next(rfin)
            header[0] = "exonID"
            header.append("gene_name")
            wfout.writerow(header)
            for row in rfin:
                gn_idx = idList.index(row[0])
                row.append(gnList[gn_idx])
                wfout.writerow(row)

########## Main ##########
wk_dir = "/Volumes/EXP337/Exp337CD25KONascent/3_DE-seq/0_DEseq2_out/pval_0.05"
os.chdir(wk_dir)

GRCm38_ref = "/Volumes/EXP337/Exp337CD25KONascent/References/GRCm38_exonid_genename.csv"
reftab = ascii.read(GRCm38_ref)
id_list = list(reftab.columns[0])
gn_list = list(reftab.columns[1])

'''
for file in glob.glob("*.csv"):
    TIDtoGN(file, id_list, gn_list)
'''

os.chdir("/Volumes/EXP337/Exp337CD25KONascent/2_count/2_Compiled_csv")
in_file = "Exp337_dupr_all_count_c5.csv"
TIDtoGN(in_file, id_list, gn_list)







