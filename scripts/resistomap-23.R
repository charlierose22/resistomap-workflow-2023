# load necessary packages and source function script
if (!exists("flag", mode = "function")) source("scripts/functions.R")
library(plyr)
library(tidyverse)
library(viridis)
library(hrbrthemes)

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
  sheet = "ARG selection") %>% 
  janitor::clean_names()
assayinformation$forward_primer = NULL
assayinformation$reverse_primer = NULL

# import moisture content
moisture_content <- read_csv("raw-data/moisture_content.csv", 
                             col_types = cols(Date = col_date(format = "%d-%m-%Y"), 
                                              Boat_Weight_g = col_skip(), 
                                              Start_Weight_g = col_skip(), 
                                              End_Weight_g = col_skip(), 
                                              Net_Start_Weight_g = col_double(), 
                                              Net_End_Weight_g = col_double())) %>% 
  janitor::clean_names()

# remove % from moisture_content_pc
moisture_content$moisture_content_pc = substr(
  moisture_content$moisture_content_pc, 1, 
  nchar(moisture_content$moisture_content_pc) - 1)

# import sample information
samples <- read_csv("~/GitHub/resistomap/raw-data/samples.csv") %>% 
  janitor::clean_names()

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
write.csv(delta_ct, "data/delta_ct_nostats.csv", row.names = TRUE)

# Turn the rownames into the first column to preserve them.
delta_ct_rownames <- rownames_to_column(delta_ct, "sample")

# pivot longer
delta_ct_long <- delta_ct_rownames %>%
  pivot_longer(cols = AY111:AY601, 
               names_to = "assay", 
               values_to = "delta_ct")

# add replicate columns again.
delta_ct_long$replicate = NA
delta_ct_long$id = NA

# fill in replicate IDs
delta_ct_long <- mutate(delta_ct_long,
                        replicate = case_when(
                          str_detect(sample, "rep1") ~ "1",
                          str_detect(sample, "rep2") ~ "2",
                          str_detect(sample, "rep3") ~ "3" ))

# fill in sample IDs.
delta_ct_long <- mutate(delta_ct_long, id = case_when(
  str_starts(sample, "A") ~ "A",
  str_starts(sample, "B") ~ "B",
  str_starts(sample, "C") ~ "C" ,
  str_starts(sample, "D") ~ "D",
  str_starts(sample, "E") ~ "E",
  str_starts(sample, "F") ~ "F",
  str_starts(sample, "G") ~ "G",
  str_starts(sample, "H") ~ "H",
  str_starts(sample, "I") ~ "I",
  str_starts(sample, "J") ~ "J",
  str_starts(sample, "K") ~ "K",
  str_starts(sample, "L") ~ "L",
  str_starts(sample, "M") ~ "M",
  str_starts(sample, "N") ~ "N",
  str_starts(sample, "O") ~ "O",
  str_starts(sample, "P") ~ "P",
  str_starts(sample, "Q") ~ "Q",
  str_starts(sample, "R") ~ "R",
  str_starts(sample, "S") ~ "S",
  str_starts(sample, "T") ~ "T",
  str_starts(sample, "U") ~ "U"))

# remove nas
nona_delta_ct <- delta_ct_long %>% drop_na()

# calculate means
means <- nona_delta_ct %>%
  group_by(pick(id, assay)) %>%
  summarise(mean = mean(delta_ct),
            std = sd(delta_ct),
            n = length(delta_ct),
            se = std/sqrt(n))

# left join for assay information and sample information
assay_means = means %>% left_join(assayinformation, by = "assay")
assay_samples_means = assay_means %>% left_join(samples, by = "id")

# swap pile length measurements for locations X and Y
assay_samples_means$slice = NA
assay_samples_means <- mutate(assay_samples_means,
                        slice = case_when(
                          str_detect(length, "half") ~ "Y",
                          str_detect(length, "quarter") ~ "Z"))
assay_samples_means$length = NULL

# Rearrange columns.
assay_samples_means <- subset(assay_samples_means, select = c(
  assay, gene, target_antibiotics_major, id, day, height, slice,
  mean, std, n, se))

# create separate data sets for location study and time series.
location_study <- assay_samples_means[grepl('\\<1\\>|\\<29\\>',
                                            assay_samples_means$day),]
time_study <- assay_samples_means[grepl('bottom',
                                        assay_samples_means$height),]
time_study <- time_study[grepl('Y', time_study$slice),]

# create csvs of annotated data
write.csv(assay_samples_means, "data/annotated_delta_ct_means.csv")
write.csv(time_study, "data/annotated_time_study.csv")
write.csv(location_study, "data/annotated_location_study.csv")

# location study plots
# Create a list of target antibiotics.
target_antibiotics <- unique(assay_samples_means$target_antibiotics_major)

# Create a loop for each target antibiotic.
for (target_antibiotic in target_antibiotics) {
  
  # Create a subset of the data for the current target antibiotic.
  data <- location_study[location_study$target_antibiotics_major == 
                                target_antibiotic, ]
  
  # Create a heatmap of the data.
  ggplot(data, aes(x = day, y = height, fill = mean)) +
    geom_tile() +
    scale_y_discrete(limits = rev) +
    scale_fill_gradient2(low = "turquoise3", high = "orange", mid = "yellow", 
                         midpoint = 11) +
    labs(x = "day", y = "height", colour = "normalised delta ct") +
    theme_ipsum(base_size = 10)
  
  # Save the linegraph to a file.
  ggsave(paste0("figures/heatmaps/location-heatmap-", target_antibiotic, ".png"))
}

# Create a loop for each target antibiotic.
for (target_antibiotic in target_antibiotics) {
  
  # Create a subset of the data for the current target antibiotic.
  data2 <- location_study[location_study$target_antibiotics_major == 
                                 target_antibiotic, ]
  
  # Create line graphs of the data.
  data2 %>% 
    group_by(slice, day) %>% 
  ggplot(aes(x = height, y = log(mean), fill = gene)) +
      geom_bar(position = "stack", stat = "identity") +
      scale_fill_viridis(discrete = T) +
      labs(x = "height", y = "normalised delta ct") +
      facet_grid(slice ~ day) +
      theme_ipsum(base_size = 10)
      
  # Save the linegraph to a file.
  ggsave(paste0("figures/bargraph/location-bargraph-", target_antibiotic, ".png"))
}

# time series plots

# Create a loop for each target antibiotic.
for (target_antibiotic in target_antibiotics) {
  
  # Create a subset of the data for the current target antibiotic.
  data3 <- time_study[time_study$target_antibiotics_major == 
                                target_antibiotic, ]
  
  # Create a heatmap of the data.
  ggplot(data3, aes(x = day, y = gene, fill = mean)) +
    geom_tile() +
    scale_y_discrete(limits = rev) +
    scale_fill_viridis(discrete = F) +
    labs(x = "day", y = "gene", colour = "normalised delta ct") +
    theme_ipsum(base_size = 10)
  
  # Save the heatmap to a file.
  ggsave(paste0("figures/heatmaps/time-heatmap-", target_antibiotic, ".png"))
}

# Create a loop for each target antibiotic.
for (target_antibiotic in target_antibiotics) {
  
  # Create a subset of the data for the current target antibiotic.
  data4 <- time_study[time_study$target_antibiotics_major == 
                                 target_antibiotic, ]
  
  # Create line graphs of the data.
  ggplot(data4, aes(x = day, 
                    y = log(mean), 
                    colour = gene)) +
    geom_jitter(width = 0.3) +
    geom_line(aes(color =  gene)) +
    labs(x = "day", y = "normalised delta ct") +
    scale_color_viridis(discrete = TRUE) +
    theme_ipsum(base_size = 10)
    
  # Save the linegraph to a file.
  ggsave(paste0("figures/linegraph/time-linegraph-", target_antibiotic, ".png"))
}

# split based on location study data
split <- split(location_study, location_study$target_antibiotics_major)
loc_aminoglycoside <- split$Aminoglycoside
loc_beta_lactam <- split$'Beta Lactam'
loc_integrons <- split$Integrons
loc_multidrug <- split$MDR
loc_mobile <- split$MGE
loc_mlsb <- split$MLSB
loc_other <- split$Other
loc_phenicol <- split$Phenicol
loc_quinolone <- split$Quinolone
loc_sulfonamide <- split$Sulfonamide
loc_tetracycline <- split$Tetracycline
loc_trimethoprim <- split$Trimethoprim
loc_vancomycin <- split$Vancomycin

# split based on time series data
split2 <- split(time_study, time_study$target_antibiotics_major)
time_aminoglycoside <- split2$Aminoglycoside
time_beta_lactam <- split2$'Beta Lactam'
time_integrons <- split2$Integrons
time_multidrug <- split2$MDR
time_mobile <- split2$MGE
time_mlsb <- split2$MLSB
time_other <- split2$Other
time_phenicol <- split2$Phenicol
time_quinolone <- split2$Quinolone
time_sulfonamide <- split2$Sulfonamide
time_tetracycline <- split2$Tetracycline
time_trimethoprim <- split2$Trimethoprim
time_vancomycin <- split2$Vancomycin

# graphs focus on tetracycline?

# can we do a double y graph to compare with moisture content?

