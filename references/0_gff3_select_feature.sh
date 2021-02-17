#!/bin/bash

awk -F "\t" '{ if ($3 == "exon") {print $0} }' Mus_musculus.GRCm38.102.gff3 > Mus_musculus.GRCm38.102.exon.gff3
awk -F "\t" '{ if ($3 == "mRNA") {print $0} }' Mus_musculus.GRCm38.102.gff3 > Mus_musculus.GRCm38.102.mRNA.gff3