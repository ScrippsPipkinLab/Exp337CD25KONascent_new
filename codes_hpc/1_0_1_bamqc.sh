#!/bin/bash

cd /gpfs/group/pipkin/Exp337/srt_flt_bam/dupr
for i in *srt_dupr_flt.bam
do
  /gpfs/home/hdiao/qualimap_v2.2.1/qualimap bamqc -bam $i
done

cd /gpfs/group/pipkin/Exp337/srt_flt_bam/non_dupr
for i in *srt_flt.bam
do
  /gpfs/home/hdiao/qualimap_v2.2.1/qualimap bamqc -bam $i
done
