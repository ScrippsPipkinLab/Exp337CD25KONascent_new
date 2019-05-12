########## Calculate TPM for nascent transcripts ##########
# Author: Huitian (Yolanda) Diao
# May 8th, 2019

######################################## Imports ########################################
library(dplyr)
library(tidyverse)

######################################## Main ########################################
if (FALSE) {
  wk.dir <- "/Volumes/EXP337/Exp337CD25KONascent/2_count/3_tpm"
  setwd(wk.dir)
  
  # Read files
  count.file <- "/Volumes/EXP337/Exp337CD25KONascent/2_count/2_compiled_csv/Exp337_all_c5.csv"
  ref.file <- "/Volumes/EXP337/Exp337CD25KONascent/References/GRCm38_exon_rmdup_srt_cb_srt_dupr.csv"
  count.tb <- read_csv(count.file)
  ref.tb <- read_csv(ref.file)
  
  # Transform tibbles
  ref.tb <- ref.tb %>%
    mutate(size = abs(start-end)) %>%
    separate(exon_id, sep="=", into=c("type", "name")) %>%
    select("name", "size")
  
  # Combine tibble
  count.tb <- inner_join(count.tb, ref.tb, by="name")
  colnames(ref.tb)
  col.names <- colnames(count.tb)
  col.numeric <- col.names[2:24]
  
  rpk.names <- c()
  # Calculate rpk
  for (col.name in col.numeric) {
    new.col.name <- paste(col.name, "rpk", sep="_")
    rpk.names <- c(rpk.names, c(new.col.name))
    count.tb <- count.tb %>%
      mutate(!!sym(new.col.name) := !!sym(col.name) / size * 1000)
  }

  # Calculate tpm
  tpm.names <- c()
  for (col.name in rpk.names) {
    new.col.name <- gsub("rpk", "tpm", col.name)
    tpm.names <- c(tpm.names, c(new.col.name))
    count.tb <- count.tb %>%
      mutate(!!sym(new.col.name) := !!sym(col.name) / (sum(!!sym(col.name)) / 1000000) )
  }
  
  # Select tpm out
  tpm.names <- c(c("name"), tpm.names)
  tpm.tb <- count.tb %>% select(one_of(tpm.names))
  
  # Calculate average tpm
  tpm.names.types <- c()
  for (name in tpm.names) {
    new.name <- gsub("_dupr_tpm", "", name)
    new.name <- gsub("_rep1", "", new.name)
    new.name <- gsub("_rep2", "", new.name)
    new.name <- gsub("_rep3", "", new.name)
    tpm.names.types <- c(tpm.names.types, c(new.name))
  }
  tpm.names.types <- tail(tpm.names.types, length(tpm.names.types) - 1)
  tpm.names.types <- unique(tpm.names.types)
  
  tpm.avg.tb <- tpm.tb %>% select(name)
  for (type in tpm.names.types) {
    # find columns of the same replicate type
    type.cols <- c()
    for (col in tpm.names) {
      if (grepl(type, col)) {
        type.cols <- c(type.cols, c(col))
      }
    }
    # Create a new tibble for the same replicate type...
    new.tb <- tpm.tb %>% select(one_of(type.cols))
    # Add type mean to new tibble
    tpm.avg.tb <- tpm.avg.tb %>% add_column(!!sym(type) := rowMeans(new.tb))
  }
  
  
  ###----- Write output
  write_csv(tpm.tb, "Exp337_dupr_tpm.csv")
  write_csv(tpm.avg.tb, "Exp337_dupr_tpm_avg.csv")
  
  
  
  
}