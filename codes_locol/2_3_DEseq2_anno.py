#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 11:25:57 2019

@author: yolandatiao
"""

########## Collect count files ##########
# Author: Huitian (Yolanda) Diao
# Apr 30, 2019
# Dependencies:
# - HTseq-count output count files
# - DEseq2

########## Import ##########
import os
import csv
import glob
from astropy.io import ascii

########## Self defined functions ##########
def CompileTables(wkDir, outName, colN):
    filelist=[]
    for root, dirs, files in os.walk(wkDir):
        for file in files:
            if file.endswith("csv"):
                filelist.append(os.path.join(root, file))
    fnames = [x.split("/")[-1].replace(".csv", "") for x in filelist]
    print(fnames)
    
    with open(outName, "w") as fout:
        wfout = csv.writer(fout, delimiter=",")
        header = ["name"] + fnames
        wfout.writerow(header)
        #print(header)
        finlist = []
        rfinlist = []
        for x in range(0, len(filelist)):
            finlist.append(open(filelist[x], "r"))
        for x in range(0, len(filelist)):
            rfinlist.append(csv.reader(finlist[x], delimiter=","))
        for x in range(0, len(filelist)):
            #print(x)
            next(rfinlist[x])
        for row in rfinlist[0]:
            newrow = [row[0], row[colN]]
            for x in range(1, len(filelist)):
                newrow.append(next(rfinlist[x])[colN])
            wfout.writerow(newrow)

def fltTables(fileList, sltFile, cutoff, direction, apdx):
    keepRows = []
    totalRowN = 0
    with open(sltFile, "r") as fin:
        rfin = csv.reader(fin, delimiter=",")
        next(rfin)
        rowN = 0
        for row in rfin:
            rowN += 1
            rowNumbers = row[1:]
            rowNumbers = [float(x) for x in rowNumbers if x!="NA"]
            if rowNumbers != []:
                if direction == "more":
                    if max(rowNumbers) >= cutoff:
                        keepRows.append(rowN)
                elif direction == "less":
                    if min(rowNumbers) <= cutoff:
                        keepRows.append(rowN)
                elif direction == "abs_more":
                    rowNumbers = [abs(x) for x in rowNumbers]
                    if max(rowNumbers) >= cutoff:
                        keepRows.append(rowN)
        totalRowN = rowN
    for f in fileList:
        outName = f.replace(".csv", "_%s-c%s.csv"%(apdx, cutoff))
        with open(f, "r") as fin:
            with open(outName, "w") as fout:
                rfin = csv.reader(fin, delimiter=",")
                wfout = csv.writer(fout, delimiter=",")
                wfout.writerow(next(rfin))
                rowN = 0
                for row in rfin:
                    rowN += 1
                    if rowN in keepRows:
                        wfout.writerow(row)
    print("Total row number: %s" %totalRowN)
    print("Kept row number: %s" %len(keepRows))
    print("Kept percentage: %s percent" %str(float(len(keepRows))/totalRowN*100))

def fltDEseqOut(fileX, log2fcC, pvalC, padjC, updn):
    outName = fileX.replace(".csv", "_log2fc-c%s-%s_pval-c%s_padj-c%s.csv"%(log2fcC, updn, pvalC, padjC))
    with open(fileX, "r") as fin:
        with open(outName, "w") as fout:
            rfin = csv.reader(fin, delimiter=",")
            wfout = csv.writer(fout, delimiter=",")
            wfout.writerow(next(rfin))
            for row in rfin:
                nu_row = row
                if row[2] == "NA":
                    nu_row[2] = 0
                if row[5] == "NA":
                    nu_row[5] = 1
                if row[6] == "NA":
                    nu_row[6] = 1
                if log2fcC == "NA" or (log2fcC != "NA" and updn == "up" and float(nu_row[2]) > log2fcC) or (log2fcC != "NA" and updn == "dn" and float(nu_row[2]) < log2fcC):
                    if pvalC == "NA" or (pvalC != "NA" and float(nu_row[5]) <= pvalC):
                        if padjC == "NA" or (padjC != "NA" and float(nu_row[6]) <= padjC):
                            wfout.writerow(row)
                        

########## Main ##########
###----- Compile DE-seq2 results
wk_dir = "/Volumes/EXP337/Exp337CD25KONascent/3_DE-seq"
os.chdir(wk_dir)

DEseq2out_dir = "/Volumes/EXP337/Exp337CD25KONascent/3_DE-seq/DEseq2_out"
CompileTables(DEseq2out_dir, "DEseq_pval.csv", 5)
CompileTables(DEseq2out_dir, "DEseq_padj.csv", 6)
CompileTables(DEseq2out_dir, "DEseq_log2fc.csv", 2)
CompileTables(DEseq2out_dir, "DEseq_baseMean.csv", 1)

###----- Filter Compiled outputs
wk_dir = "/Volumes/EXP337/Exp337CD25KONascent/3_DE-seq/1_Compiled_out"
os.chdir(wk_dir)

file_list = ["DEseq_baseMean.csv", "DEseq_log2fc.csv", "DEseq_padj.csv", "DEseq_pval.csv"]
fltTables(file_list, "DEseq_padj.csv", 0.05, "less", "padj")

file_list = ["DEseq_baseMean.csv", "DEseq_log2fc.csv", "DEseq_padj.csv", "DEseq_pval.csv"]
fltTables(file_list, "DEseq_pval.csv", 0.05, "less", "pval")

file_list = ["DEseq_baseMean.csv", "DEseq_log2fc.csv", "DEseq_padj.csv", "DEseq_pval.csv"]
fltTables(file_list, "DEseq_log2fc.csv", 1, "abs_more", "log2fc")

file_list = ["DEseq_baseMean.csv", "DEseq_log2fc.csv", "DEseq_padj.csv", "DEseq_pval.csv"]
fltTables(file_list, "DEseq_log2fc.csv", 2, "abs_more", "log2fc")

file_list = ["DEseq_baseMean_log2fc-c1.csv", "DEseq_log2fc_log2fc-c1.csv", 
             "DEseq_padj_log2fc-c1.csv", "DEseq_pval_log2fc-c1.csv"]
fltTables(file_list,  "DEseq_pval_log2fc-c1.csv", 0.05, "less", "log2fc")

###----- Filter each file
wk_dir = "/Volumes/EXP337/Exp337CD25KONascent/3_DE-seq/DEseq2_out"
os.chdir(wk_dir)

for file in glob.glob("*h.csv"):
    fltDEseqOut(file, 0, 0.05, "NA", "up")
    fltDEseqOut(file, 0, 0.05, "NA", "dn")
    fltDEseqOut(file, 0, "NA", 0.05, "up")
    fltDEseqOut(file, 0, "NA", 0.05, "dn")





