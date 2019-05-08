########## Calculate TPM for nascent transcripts ##########
# Author: Huitian (Yolanda) Diao
# May 8th, 2019

######################################## Imports ########################################
library(dplyr)
library(tidyverse)

######################################## Main ########################################
if (FALSE) {
  wk.dir <- "/Volumes/EXP337/Exp337CD25KONascent/3_tpm"
  setwd(wk.dir)
  
  # Read files
  count.file <- "/Volumes/EXP337/Exp337CD25KONascent/2_count/2_Compiled_csv/Exp337_dupr_all_count_c5_GN.csv"
  ref.file <- "/Volumes/EXP337/Exp337CD25KONascent/References/GRCm38_exon_rmdup.csv"
  count.tb <- read_csv(count.file)
  ref.tb <- read_csv(ref.file)
  
  # Transform tibbles
  ref.tb <- ref.tb %>%
    mutate(size = abs(start-end)) %>%
    separate(attribute, sep="=", into=c("type", "exonID")) %>%
    select("exonID", "size")
  nrow(ref.tb)
  
  # Combine tibble
  count.tb <- inner_join(count.tb, ref.tb, by="exonID")
  colnames(ref.tb)
  col.names <- colnames(count.tb)
  col.names
  
  # Calculate rpk
  for (col.name in col.names) {
    new.col.name <- paste(col.name, "rpk", sep="_")
    count.tb <- count.tb %>%
      mutate(!!sym(new.col.name) := !!sym(col.name) / )
  }
  
  rpk.tb <- count.tb %>%
  
  
  
  
  
  
}