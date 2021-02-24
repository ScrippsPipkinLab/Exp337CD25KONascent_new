#!/use/bin/env python

import csv
import sys

print("----------------------------------------")
print("Useage:")
print("$python Ensembl_bed_to_usce_bed.py ensembl_bed_file")

infile = sys.argv[1]
outfile = infile.replace(".bed","") + ".chr.bed"
outfile = infile.replace(".bdg","") + ".chr.bdg"

chr_list = [str(x) for x in list(range(1,21))] + ["X", "Y"]

with open(infile, "r") as fin:
	with open(outfile, "w") as fout:
		rfin = csv.reader(fin, delimiter="\t")
		wfout = csv.writer(fout, delimiter="\t")
		for row in rfin:
			if str(row[0]) in chr_list:
				newrow = ["chr%s"%row[0].replace("MT","M")] + row[1:]
				wfout.writerow(newrow)
