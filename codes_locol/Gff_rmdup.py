# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

########## Import ##########
import csv
import os

wk_dir = "/Volumes/EXP337/Exp337CD25KONascent/References"
os.chdir(wk_dir)

in_file = "GRCm38_exon.gff"
out_file = "GRCm38_exon_rmdup.gff"

with open(in_file, "r") as fin:
    with open(out_file, "w") as fout:
        rfin = csv.reader(fin, delimiter = "\t")
        wfout = csv.writer(fout, delimiter = "\t")
        chr_n = 0
        s_pos = 0
        e_pos = 0
        dup = True
        for row in rfin:
            row_chr_n = row[0]
            row_s_pos = row[3]
            row_e_pos = row[4]
            if (row_chr_n == chr_n) and (row_s_pos == s_pos) and (row_e_pos == e_pos):
                dup = True
            else:
                dup = False
            if dup == False:
                wfout.writerow(row)         
            chr_n = row_chr_n
            s_pos = row_s_pos
            e_pos = row_e_pos
    
            