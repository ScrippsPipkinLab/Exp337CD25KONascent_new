#!/usr/bin/env python3
# -*- coding: utf-8 -*-

########## Collect count files ##########
# Author: Huitian (Yolanda) Diao
# Apr 30, 2019
# Dependencies:
# - HTseq-count output count files

########## Import ##########
import os
import csv
import glob
from astropy.io import ascii

########## Self defined functions ##########
def count_to_csv(countFile):
    #countFile = "/Volumes/Yolanda/Exp337CD25KONascent/2_count/dupr/0_featureCounts/337-1_dupr_f.txt"
    csvFile = countFile.replace(".txt", ".csv")
    with open(countFile, "r") as fin:
        with open(csvFile, "w") as fout:
            rfin = csv.reader(fin, delimiter="\t")
            wfout = csv.writer(fout, delimiter=",")
            next(rfin)
            next(rfin)
            newheader = ["ID", "count"]
            wfout.writerow(newheader)
            for row in rfin:
                newrow = [row[0], row[6]]
                wfout.writerow(newrow)

def CompileCount(wkDir, outName):
    #wkDir = "/Volumes/Yolanda/Exp337CD25KONascent/2_count/dupr/1_csv/F"
    #outName = "test.out.csv"
    filelist=[]
    for root, dirs, files in os.walk(wkDir):
        for file in files:
            if not file.startswith("._"):
                if file.endswith("f.csv") or file.endswith("r.csv"):
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
            newrow = [row[0], int(row[1])]
            for x in range(1, len(filelist)):
                newrow.append(int(next(rfinlist[x])[1]))
            wfout.writerow(newrow)

def filterWrongDir(inFile, inStrand, refFile):
    #refFile = "/Volumes/EXP337/Exp337CD25KONascent/References/GRCm38_exon_rmdup_srt_cb_srt_dupr.csv"
    #inStrand = -1
    outFile = inFile.replace(".csv", "_flt.csv")
    
    refTab = ascii.read(refFile)        
    refIDs = list(refTab.columns[8])
    refIDs = [x.replace("ID=", "") for x in refIDs]
    refstrand = list(refTab.columns[6])
    refstrand_cor = [idx for idx, x in enumerate(refstrand) if x==inStrand]
    refIDs_cor = [refIDs[x] for x in refstrand_cor]
    refIDs_cor[0]
    
    with open(inFile, "r") as fin:
        with open(outFile, "w") as fout:
            rfin = csv.reader(fin, delimiter=",")
            wfout = csv.writer(fout, delimiter=",")
            wfout.writerow(next(rfin))
            rowN = 0
            filterN = 0
            keepN = 0
            for row in rfin:
                rowN += 1
                if row[0] in refIDs_cor:
                    wfout.writerow(row)
                    keepN += 1
                else:
                    filterN += 1
    flt_pctg = float(filterN)/rowN*100
    print("Total row number: %s" %rowN)
    print("Filterd out %s percent" %str(flt_pctg))
    print("Total row number left: %s" %keepN)
            
def filterCount(inFile, cutoff):
    outFile = inFile.replace(".csv", "_c%s.csv"%cutoff)
    with open(inFile, "r") as fin:
        with open(outFile, "w") as fout:
            rfin = csv.reader(fin, delimiter=",")
            wfout = csv.writer(fout, delimiter=",")
            wfout.writerow(next(rfin))
            rowN = 0
            filterN = 0
            keepN = 0
            for row in rfin:
                rowN += 1
                rownum = [int(i) for i in row[1:]]
                if max(rownum) > cutoff:
                    wfout.writerow(row)
                    keepN += 1
                else:
                    filterN += 1
    flt_pctg = float(filterN)/rowN*100
    print("Total row number: %s" %rowN)
    print("Filterd out %s percent" %str(flt_pctg))
    print("Total row number left: %s" %keepN)
    

########## Main ##########

##### Convert counts
wk_dir = "/Volumes/Yolanda/Exp337CD25KONascent/2_count/dupr/0_featureCounts"
os.chdir(wk_dir)

for file in glob.glob("*.txt"):
    count_to_csv(file)

wk_dir = "/Volumes/Yolanda/Exp337CD25KONascent/2_count/nondupr/0_featureCounts"
os.chdir(wk_dir)

for file in glob.glob("*.txt"):
    count_to_csv(file)



##### Compile counts
wk_dir = "/Volumes/Yolanda/Exp337CD25KONascent/2_count/dupr/1_csv/F"
os.chdir(wk_dir)
out_name = "Exp337_dupr_F_count.csv"
CompileCount(wk_dir, out_name)

wk_dir = "/Volumes/Yolanda/Exp337CD25KONascent/2_count/dupr/1_csv/R"
os.chdir(wk_dir)
out_name = "Exp337_dupr_R_count.csv"
CompileCount(wk_dir, out_name)

wk_dir = "/Volumes/Yolanda/Exp337CD25KONascent/2_count/nondupr/1_csv/F"
os.chdir(wk_dir)
out_name = "Exp337_F_count.csv"
CompileCount(wk_dir, out_name)

wk_dir = "/Volumes/Yolanda/Exp337CD25KONascent/2_count/nondupr/1_csv/R"
os.chdir(wk_dir)
out_name = "Exp337_R_count.csv"
CompileCount(wk_dir, out_name)

##### Filter out wrong direction matches
# Turns out c5 files did not contain any wrong direction mapping, skip
'''
wk_dir = "/Volumes/EXP337/Exp337CD25KONascent/2_count/2_Compiled_csv"
os.chdir(wk_dir)
refFile = "/Volumes/EXP337/Exp337CD25KONascent/References/GRCm38_exon_rmdup_srt_cb_srt_dupr.csv"
filterWrongDir("Exp337_dupr_F_count.csv", -1, refFile)
filterWrongDir("Exp337_dupr_R_count.csv", 1, refFile)
filterWrongDir("Exp337_dupr_F_count_c5.csv", "-", refFile)
filterWrongDir("Exp337_dupr_R_count_c5.csv", "+", refFile)
'''

##### Filter count table
wk_dir = "/Volumes/Yolanda/Exp337CD25KONascent/2_count/dupr/2_compiled_csv"
os.chdir(wk_dir)

filterCount("Exp337_dupr_F_count.csv", 10)
filterCount("Exp337_dupr_R_count.csv", 10)


wk_dir = "/Volumes/Yolanda/Exp337CD25KONascent/2_count/nondupr/2_compiled_csv"
os.chdir(wk_dir)

filterCount("Exp337_F_count.csv", 10)
filterCount("Exp337_R_count.csv", 10)


























