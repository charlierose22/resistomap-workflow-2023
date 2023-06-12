# load necessary packages and source function script
if (!exists("flag", mode = "function")) source("scripts/functions.R")
library(plyr)
library(tidyverse)

# import raw dataset
rawdata <- readxl::read_excel("raw-data/raw_results_drying_study_resistomap.xlsx",
                              sheet = "1A-21A", 
                              col_types = c("skip", 
                                            "skip", 
                                            "text", 
                                            "text", 
                                            "skip", 
                                            "numeric", 
                                            "numeric", 
                                            "numeric", 
                                            "text")) %>% 
  janitor::clean_names()

# import assay information
assayinformation <- readxl::read_excel(
  "~/GitHub/resistomap/raw-data/raw_results_drying_study_resistomap.xlsx",
  sheet = "ARG selection")

# import sample information
samples <- read_csv("~/GitHub/resistomap/raw-data/samples.csv")

# remove unsatisfactory flags
flags_removed <- flag(rawdata)
flags_removed$flags = NULL

# add individual sample locations for removing "solo" results
# "solo" results are when out of the three repeats, there was only amplification in one sample
flags_removed$id = NA
flags_removed$replicate = NA

# fill in replicate IDs
flags_removed <- mutate(flags_removed,
                        replicate = case_when(
                          str_detect(sample, "rep1") ~ "1",
                          str_detect(sample, "rep2") ~ "2",
                          str_detect(sample, "rep3") ~ "3" ))

# fill in sample IDs.
flags_removed <- mutate(flags_removed, id = case_when(
  str_starts(sample, "1A") ~ "A",
  str_starts(sample, "2A") ~ "B",
  str_starts(sample, "3A") ~ "C" ,
  str_starts(sample, "4A") ~ "D",
  str_starts(sample, "5A") ~ "E",
  str_starts(sample, "6A") ~ "F",
  str_starts(sample, "7A") ~ "G",
  str_starts(sample, "8A") ~ "H",
  str_starts(sample, "9A") ~ "I",
  str_starts(sample, "10A") ~ "J",
  str_starts(sample, "11A") ~ "K",
  str_starts(sample, "12A") ~ "L",
  str_starts(sample, "13A") ~ "M",
  str_starts(sample, "14A") ~ "N",
  str_starts(sample, "15A") ~ "O",
  str_starts(sample, "16A") ~ "P",
  str_starts(sample, "17A") ~ "Q",
  str_starts(sample, "18A") ~ "R",
  str_starts(sample, "19A") ~ "S",
  str_starts(sample, "20A") ~ "T",
  str_starts(sample, "21A") ~ "U"))

# remove "solo" results
solo_removed <- ddply(flags_removed, c("assay", "id"),
                      function(d) {if (nrow(d) > 1) d else NULL})

# remove Tm and Efficiency columns as we're focusing on cycle threshold
solo_removed$tm = NULL
solo_removed$efficiency = NULL

# mutate the sample column to not begin with a number for easier coding and recognition.
mutated_id <- mutate(solo_removed, sample = case_when(
  str_starts(sample, "1A-rep1") ~ "A-rep1",
  str_starts(sample, "1A-rep2") ~ "A-rep2",
  str_starts(sample, "1A-rep3") ~ "A-rep3",
  str_starts(sample, "2A-rep1") ~ "B-rep1",
  str_starts(sample, "2A-rep2") ~ "B-rep2",
  str_starts(sample, "2A-rep3") ~ "B-rep3",
  str_starts(sample, "3A-rep1") ~ "C-rep1",
  str_starts(sample, "3A-rep2") ~ "C-rep2",
  str_starts(sample, "3A-rep3") ~ "C-rep3",
  str_starts(sample, "4A-rep1") ~ "D-rep1",
  str_starts(sample, "4A-rep2") ~ "D-rep2",
  str_starts(sample, "4A-rep3") ~ "D-rep3",
  str_starts(sample, "5A-rep1") ~ "E-rep1",
  str_starts(sample, "5A-rep2") ~ "E-rep2",
  str_starts(sample, "5A-rep3") ~ "E-rep3",
  str_starts(sample, "6A-rep1") ~ "F-rep1",
  str_starts(sample, "6A-rep2") ~ "F-rep2",
  str_starts(sample, "6A-rep3") ~ "F-rep3",
  str_starts(sample, "7A-rep1") ~ "G-rep1",
  str_starts(sample, "7A-rep2") ~ "G-rep2",
  str_starts(sample, "7A-rep3") ~ "G-rep3",
  str_starts(sample, "8A-rep1") ~ "H-rep1",
  str_starts(sample, "8A-rep2") ~ "H-rep2",
  str_starts(sample, "8A-rep3") ~ "H-rep3",
  str_starts(sample, "9A-rep1") ~ "I-rep1",
  str_starts(sample, "9A-rep2") ~ "I-rep2",
  str_starts(sample, "9A-rep3") ~ "I-rep3",
  str_starts(sample, "10A-rep1") ~ "J-rep1",
  str_starts(sample, "10A-rep2") ~ "J-rep2",
  str_starts(sample, "10A-rep3") ~ "J-rep3",
  str_starts(sample, "11A-rep1") ~ "K-rep1",
  str_starts(sample, "11A-rep2") ~ "K-rep2",
  str_starts(sample, "11A-rep3") ~ "K-rep3",
  str_starts(sample, "12A-rep1") ~ "L-rep1",
  str_starts(sample, "12A-rep2") ~ "L-rep2",
  str_starts(sample, "12A-rep3") ~ "L-rep3",
  str_starts(sample, "13A-rep1") ~ "M-rep1",
  str_starts(sample, "13A-rep2") ~ "M-rep2",
  str_starts(sample, "13A-rep3") ~ "M-rep3",
  str_starts(sample, "14A-rep1") ~ "N-rep1",
  str_starts(sample, "14A-rep2") ~ "N-rep2",
  str_starts(sample, "14A-rep3") ~ "N-rep3",
  str_starts(sample, "15A-rep1") ~ "O-rep1",
  str_starts(sample, "15A-rep2") ~ "O-rep2",
  str_starts(sample, "15A-rep3") ~ "O-rep3",
  str_starts(sample, "16A-rep1") ~ "P-rep1",
  str_starts(sample, "16A-rep2") ~ "P-rep2",
  str_starts(sample, "16A-rep3") ~ "P-rep3",
  str_starts(sample, "17A-rep1") ~ "Q-rep1",
  str_starts(sample, "17A-rep2") ~ "Q-rep2",
  str_starts(sample, "17A-rep3") ~ "Q-rep3",
  str_starts(sample, "18A-rep1") ~ "R-rep1",
  str_starts(sample, "18A-rep2") ~ "R-rep2",
  str_starts(sample, "18A-rep3") ~ "R-rep3",
  str_starts(sample, "19A-rep1") ~ "S-rep1",
  str_starts(sample, "19A-rep2") ~ "S-rep2",
  str_starts(sample, "19A-rep3") ~ "S-rep3",
  str_starts(sample, "20A-rep1") ~ "T-rep1",
  str_starts(sample, "20A-rep2") ~ "T-rep2",
  str_starts(sample, "20A-rep3") ~ "T-rep3",
  str_starts(sample, "21A-rep1") ~ "U-rep1",
  str_starts(sample, "21A-rep2") ~ "U-rep2",
  str_starts(sample, "21A-rep3") ~ "U-rep3"))

# pivot the table, so assay is across the top
mutated_wide <- mutated_id %>%
  pivot_wider(names_from = assay, values_from = ct)

mutated_wide$id = NULL
mutated_wide$replicate = NULL

# pivot the table back to longer.
mutated_long <- mutated_wide %>%
  pivot_longer(cols = AY1:AY601, 
               names_to = "assay", 
               values_to = "ct")

# pivot the table again, so sample is across the top, rather than assay.
mutated_wide_2 <- mutated_long %>%
  pivot_wider(
    names_from = sample, values_from = ct)

# transpose the table.
transposed_ct <- t(mutated_wide_2[2:64])

# Make sure the top row is the Assay codes, by first extracting as a list.
assay_names <- t(mutated_wide_2[1])
colnames(transposed_ct) <- as.character(assay_names[1,])

# Calculate delta ct and make sure the output is as a data frame.
df_transposed_ct <- as.data.frame(transposed_ct)
delta_ct <- df_transposed_ct[ , 2:70] - df_transposed_ct[ , "AY1"]

# create delta ct csv
write.csv(delta_ct, "delta_ct_nostats.csv", row.names = TRUE)
