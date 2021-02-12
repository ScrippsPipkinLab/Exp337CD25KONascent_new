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
import pandas as pd

########## Self defined functions ##########
def TIDtoGN(fileX, idList, gnList):    
    outName = fileX.replace(".csv", "_GN.csv")
    with open(fileX, "r") as fin:
        with open(outName, "w") as fout:
            rfin = csv.reader(fin, delimiter=",")
            wfout = csv.writer(fout, delimiter=",")
            header = next(rfin)
            header[0] = "gene_name"
            wfout.writerow(header)
            for row in rfin:
                gn_idx = idList.index(row[0])
                row[0] = gnList[gn_idx]
                wfout.writerow(row)

########## Main ##########
GRCm38_ref = "/media/pipkin/Yolanda/Exp337CD25KONascent/References/GRCm38_exonid_genename.csv"
reftab = pd.read_csv(GRCm38_ref)
id_list = reftab["Exon stable ID"].tolist()
gn_list = reftab["Gene name"].tolist()

###----- Duplicate removed
if False:
    wk_dir = "/media/pipkin/Yolanda/Exp337CD25KONascent/3_DE-seq/dupr"
    os.chdir(wk_dir)
    filelist = glob.glob("*h.csv")
    TIDtoGN(filelist[0], id_list, gn_list)

    tb_0 = pd.read_csv(filelist[0])
    names_0 = tb_0.iloc[:, 0]
    tb_0_new = pd.read_csv(filelist[0].replace(".csv", "_GN.csv"))
    genes_0 = list(tb_0_new.iloc[:, 0])
    for i in range(0, len(filelist)):
        out_name = filelist[i].replace(".csv", "_addGN.csv")
        tb_i = pd.read_csv(filelist[i])
        names_i = tb_i.iloc[:, 0]
        names_match = list(names_0 == names_i)
        if False not in names_match:
            tb_i["gene_name"] = genes_0
            tb_i.to_csv(out_name, index=False)
        else:
            print("File %s has different name order..." %i)


###----- Duplicate not removed
wk_dir = '/media/pipkin/Yolanda/Exp337CD25KONascent/2_DE-seq/0_original/nondupr/0_original'
os.chdir(wk_dir)
filelist = glob.glob("*h.csv")
TIDtoGN(filelist[0], id_list, gn_list)

tb_0 = pd.read_csv(filelist[0])
names_0 = tb_0.iloc[:, 0]
tb_0_new = pd.read_csv(filelist[0].replace(".csv", "_GN.csv"))
genes_0 = list(tb_0_new.iloc[:, 0])
for i in range(0, len(filelist)):
    out_name = filelist[i].replace(".csv", "_addGN.csv")
    tb_i = pd.read_csv(filelist[i])
    names_i = tb_i.iloc[:, 0]
    names_match = list(names_0 == names_i)
    if False not in names_match:
        tb_i["gene_name"] = genes_0
        tb_i.to_csv(out_name, index=False)
    else:
        print("File %s has different name order..." %i)


