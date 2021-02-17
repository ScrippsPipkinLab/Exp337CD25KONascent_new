#!/bin/bash

sort -k1,1 -k4n Mus_musculus.GRCm38.102.exon.simp.rmdup.gff3 > Mus_musculus.GRCm38.102.exon.simp.rmdup.srt.gff3
sort -k1,1 -k4n Mus_musculus.GRCm38.102.mRNA.simp.rmdup.gff3 > Mus_musculus.GRCm38.102.mRNA.simp.rmdup.srt.gff3