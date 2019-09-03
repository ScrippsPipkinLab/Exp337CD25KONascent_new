#!/bin/bash

in_dir=/gpfs/group/pipkin/Exp337/trimmed_cb_fq
cd $in_dir

module load fastqc

for i in *.fq
do
  fastqc $i
done
