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
    group_by(slice, day) %>% 
    ggplot(aes(x = height, y = mean, fill = gene)) +
    geom_bar(position = "stack", stat = "identity") +
    scale_fill_viridis(discrete = T) +
    labs(x = "height", y = "normalised delta ct") +
    facet_grid(slice ~ day) +
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
    labs(x = "day", y = "gene", colour = "normalised delta ct") +
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
                    colour = gene)) +
    geom_point() +
    geom_errorbar(aes(x = day,
                      ymin = mean - se,
                      ymax = mean + se),
                  width = .6) +
    geom_line(aes(color =  gene)) +
    labs(x = "day", y = "normalised delta ct") +
    scale_color_viridis(discrete = TRUE) +
    theme_ipsum(base_size = 10)
  
  # Save the linegraph to a file.
  ggsave(paste0("figures/linegraph/time-error-linegraph-", 
                target_antibiotic, ".png"), width = 6, height = 5)
}