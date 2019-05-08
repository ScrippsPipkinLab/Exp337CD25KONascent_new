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

###----- Remove duplicates
in_file = "GRCm38_exon.gff"
out_file = "GRCm38_exon_rmdup.gff"
out_file2 = "GRCm38_exon_rmdup.csv"
out2_header = ["Transcript stable ID", "Exon stable ID", "Chromosome/scaffold name", 
               "Strand", "Transcript start (bp)", "Transcript end (bp)", 
               "Transcription start site (TSS)", 
               "Transcript length (including UTRs and CDS)", "Gene name"]

with open(in_file, "r") as fin:
    with open(out_file, "w") as fout:
        with open(out_file2, "w") as fout2:
            rfin = csv.reader(fin, delimiter = "\t")
            wfout = csv.writer(fout, delimiter = "\t")
            wfout2 = csv.writer(fout2, delimiter = ",")
            wfout2.writerow(out2_header)
            chr_n = 0
            s_pos = 0
            e_pos = 0
            dup = True
            next(rfin)
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
                    wfout2.writerow(row)
                chr_n = row_chr_n
                s_pos = row_s_pos
                e_pos = row_e_pos
    
###----- Combine exons
in_file = "GRCm38_exon_rmdup_srt.csv"
out_file = "GRCm38_exon_rmdup_srt_cb.csv"
out_file2 = "GRCm38_exon_rmdup_srt_cb.gff"

def overlap(listx):
    #listx = [[1,3], [5,10], [3,6]]
    listx.sort()
    list_out = [listx[0]]
    for x in listx:
        x_len = x[1] - x[0]
        overlap = False
        for out_idx in range(0, len(list_out)):
            out = list_out[out_idx]
            out_len = out[1] - out[0]
            max_dist = max(x + out) - min(x + out)
            if max_dist <= x_len + out_len:
                overlap = True
                list_out[out_idx] = [min(x + out), max(x + out)]
        if overlap == False:
            list_out.append(x)
    return(list_out)


with open(in_file, "r") as fin:
    with open(out_file, "w") as fout:
        with open(out_file2, "w") as fout2:
            rfin = csv.reader(fin, delimiter = ",")
            wfout = csv.writer(fout, delimiter = ",")
            wfout2 = csv.writer(fout2, delimiter = "\t")
            header = next(rfin)
            wfout.writerow(header)
            exonID = ""
            st_ed = []
            last_row = []
            for row in rfin:
                if row[8] != exonID:
                    if len(st_ed) >= 1:
                        st_ed = overlap(st_ed)
                        if len(st_ed) > 1:
                            print("Error, exon non-overlap: %s"%exonID)
                            print(st_ed)
                        else:
                            new_row = last_row[:3] + st_ed[0] + last_row[5:]
                            wfout.writerow(new_row)
                            wfout2.writerow(new_row)
                    else:
                        print("Error, no start or end site: %s"%exonID)
                    # Initiate new exon
                    exonID = row[8]
                    st_ed = []
                    st_ed.append([int(row[3]), int(row[4])])
                    last_row = row
                else:
                    st_ed.append([int(row[3]), int(row[4])])
                    last_row = row
            if len(st_ed) >= 1:
                st_ed = overlap(st_ed)
                if len(st_ed) > 1:
                    print("Error, exon non-overlap: %s"%exonID)
                    print(st_ed)
                else:
                    new_row = last_row[:3] + st_ed[0] + last_row[5:]
                    wfout.writerow(new_row)
                    wfout2.writerow(new_row)
                    
                    
                
















