#!/bin/bash

cd /Volumes/EXP337/Exp337CD25KONascent/References

awk -F "," \
'{if ($4==-1) \
{print $19 "\t" "." "\t" $14 "\t" $4 "\t" $5 "\t" "." "\t" "-" "\t" "." "\t" "ID=" $16} \
else \
{print $19 "\t" "." "\t" $14 "\t" $4 "\t" $5 "\t" "." "\t" "+" "\t" "." "\t" "ID=" $16} }' \
mart_export_20190509.csv > GRCm38.gff

