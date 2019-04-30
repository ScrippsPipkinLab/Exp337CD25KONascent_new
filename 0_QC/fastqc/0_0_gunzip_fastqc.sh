#!/bin/bash

wk_dir=/gpfs/group/pipkin/Exp337
cd $wk_dir

module load fastqc

for i in *.gz
do
  gunzip $i
done

for i in *.fastq
do
  fastqc $i
doen
