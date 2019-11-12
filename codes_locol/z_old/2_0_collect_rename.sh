#!/bin/bash

wk_dir=/Volumes/EXP337/Exp337CD25KONascent/2_count/1_dupr_count_new
cd $wk_dir

sample_names=(WT_0h_rep1 WT_6h_rep1 WT_24h_rep1 WT_48h_rep1 KO_0h_rep1 KO_6h_rep1 KO_24h_rep1 KO_48h_rep1 WT_0h_rep2 WT_6h_rep2 WT_24h_rep2 WT_48h_rep2 KO_0h_rep2 KO_6h_rep2 KO_24h_rep2 KO_48h_rep2 WT_0h_rep3 WT_6h_rep3 WT_24h_rep3 WT_48h_rep3 KO_0h_rep3 KO_6h_rep3 KO_24h_rep3 KO_48h_rep3)

for i in $(seq 0 23)
do
  z=`expr $i + 1`
  mv 337-${z}_dupr_neg.count ${sample_names[$i]}_dupr_F.count
  mv 337-${z}_dupr_pos.count ${sample_names[$i]}_dupr_R.count
done