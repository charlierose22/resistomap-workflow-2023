# define function to remove results with certain flags.
flag <- function(data) 
  {data[!grepl('CurveFitFail|MultipleMeltPeak|NoAmplification',
                                     data$flags),]}

# define function to strip a prefix from a string.
prefix_strip <- function(str, prefix) 
{substr(str, 1 + nchar(prefix), nchar(str))}

# define function for heatmaps
heatmap <- function(data, x, y, midpoint, xlabs, ylabs) {
ggplot(data, aes({{x}}, {{y}}, fill = mean)) +
  geom_tile() +
  scale_y_discrete(limits = rev) +
  scale_fill_gradient2(low = "turquoise3", high = "orange", mid = "yellow", 
                       midpoint) +
  labs(x = "xlabs", y = "ylabs", colour = "normalised delta ct") +
  theme_bw(base_size = 10) +
  theme(panel.grid.major = element_line(colour = "gray80"),
        panel.grid.minor = element_line(colour = "gray80"),
        axis.text.x = element_text(angle = 90),
        legend.text = element_text(family = "serif", 
                                   size = 10), 
        axis.text = element_text(family = "serif", 
                                 size = 10),
        axis.title = element_text(family = "serif",
                                  size = 10, face = "bold", colour = "gray20"),
        legend.title = element_text(size = 10,
                                    family = "serif"),
        plot.background = element_rect(colour = NA,
                                       linetype = "solid"), 
        legend.key = element_rect(fill = NA)) + labs(fill = "intensity")}