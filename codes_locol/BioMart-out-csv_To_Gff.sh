#!/bin/bash

cd /Volumes/EXP337/Exp337CD25KONascent/References

awk -F "," \
'{if ($4==-1) \
{print $3 "\t" "." "\t" "exon" "\t" $5 "\t" $6 "\t" "." "\t" "-" "\t" "." "\t" "ID=" $2} \
else \
{print $3 "\t" "." "\t" "exon" "\t" $5 "\t" $6 "\t" "." "\t" "+" "\t" "." "\t" "ID=" $2} }' \
mart_export.csv > GRCm38_exon.gff

