# location study plots
# Create a list of target antibiotics.
target_antibiotics <- unique(assay_samples_means$target_antibiotics_major)

# Create a loop for each target antibiotic.
for (target_antibiotic in target_antibiotics) {
  
  # Create a subset of the data for the current target antibiotic.
  data <- location_study[location_study$target_antibiotics_major == 
                           target_antibiotic, ]
  
  # Create a heatmap of the data.
  ggplot(data, aes(x = day, y = gene, fill = mean)) +
    geom_tile() +
    scale_y_discrete(limits = rev) +
    scale_fill_viridis(discrete = F) +
    labs(x = "day", y = "gene", colour = "abundance") +
    facet_grid(length ~ height) +
    theme_ipsum(base_size = 10)
  
  # Save the linegraph to a file.
  ggsave(paste0("figures/heatmaps/location-heatmap-", 
                target_antibiotic, ".png"), width = 7, height = 7)
}

# Create a loop for each target antibiotic.
for (target_antibiotic in target_antibiotics) {
  
  # Create a subset of the data for the current target antibiotic.
  data2 <- location_study[location_study$target_antibiotics_major == 
                            target_antibiotic, ]
  
  # Create bar graphs of the data.
  data2 %>% 
    group_by(length, height) %>% 
    ggplot(aes(x = day, y = mean, fill = gene)) +
    geom_bar(position = "stack", stat = "identity") +
    scale_fill_viridis(discrete = T) +
    labs(x = "day", y = "abundance", fill = "gene") +
    facet_grid(length ~ height) +
    theme_ipsum(base_size = 10)
  
  # Save the linegraph to a file.
  ggsave(paste0("figures/bargraph/location-bargraph-", 
                target_antibiotic, ".png"), width = 7, height = 7)
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
    labs(x = "day", y = "gene", colour = "abundance") +
    theme_ipsum(base_size = 10)
  
  # Save the heatmap to a file.
  ggsave(paste0("figures/heatmaps/time-heatmap-", 
                target_antibiotic, ".png"), width = 6, height = 5)
}

# Create a loop for each target antibiotic.
for (target_antibiotic in target_antibiotics) {
  
  # Create a subset of the data for the current target antibiotic.
  data4 <- time_study[time_study$target_antibiotics_major == 
                        target_antibiotic, ]
  
  # Create line graphs of the data.
  ggplot(data4, aes(x = day, 
                    y = mean, 
                    color = gene)) +
    geom_point() +
    geom_errorbar(aes(x = day,
                      ymin = mean - se,
                      ymax = mean + se),
                  width = .6) +
    geom_line(aes(color =  gene)) +
    labs(x = "day", y = "abundance", color = "gene") +
    scale_color_viridis(discrete = T) +
    theme_ipsum(base_size = 10)
  
  # Save the linegraph to a file.
  ggsave(paste0("figures/linegraph/time-error-linegraph-", 
                target_antibiotic, ".png"), width = 6, height = 5)
}

write.csv(samples_16S_plot, 
          "data/16S_abundance.csv", 
          row.names = FALSE)

# create a graph for total abundance of 16S data across samples.
ggplot(location_16S, aes(x = day, y = mean, fill = length)) +
  geom_col(position = 'dodge') +
  labs(x = "day", y = "ct", fill = "location") +
  scale_color_viridis(discrete = T) +
  facet_wrap(~ height) +
  theme_ipsum(base_size = 10)
ggsave(paste0("figures/16S-location.png"), width = 8, height = 5)
ggplot(time_16S, aes(x = day, y = mean, fill = length)) +
  geom_col(position = 'dodge') +
  labs(x = "day", y = "ct", fill = "location") +
  scale_color_viridis(discrete = T) +
  facet_wrap(~ height) +
  theme_ipsum(base_size = 10)
ggsave(paste0("figures/16S-time.png"), width = 8, height = 5)

