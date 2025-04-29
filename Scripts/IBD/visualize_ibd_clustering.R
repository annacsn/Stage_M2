library(tidygraph)
library(ggraph)
library(ggplot2)
library(readr)
library(RColorBrewer)
library(cowplot)
library(dplyr)
library(gridExtra)
library(grid)

# Set working directory
setwd("/Volumes/ANNA_DD/Focus_R539T/hmmIBD/")

# Load data
ibd_data <- read_tsv("IBD_fract.txt")
metadata <- read_tsv("metadata_R539T.txt")

# Clean sample names (lowercase and trim whitespace)
metadata$Sample <- tolower(trimws(metadata$Sample))
ibd_data$sample1 <- tolower(trimws(ibd_data$sample1))
ibd_data$sample2 <- tolower(trimws(ibd_data$sample2))

# Filter IBD pairs above the threshold and present in metadata
ibd_threshold <- 0.49
filtered_ibd <- ibd_data %>%
  filter(IBD_value >= ibd_threshold) %>%
  filter(sample1 %in% metadata$Sample & sample2 %in% metadata$Sample)

# Build graph: nodes = metadata, edges = filtered IBD pairs
graph_tbl <- tbl_graph(
  nodes = metadata,
  edges = filtered_ibd[, c("sample1", "sample2")],
  directed = FALSE
)

# Define colors by Year and shapes by Location
year_colors <- brewer.pal(length(unique(metadata$Year)), "Set1")
names(year_colors) <- sort(unique(metadata$Year))

location_shapes <- setNames(
  0:(length(unique(metadata$Location)) - 1),
  sort(unique(metadata$Location))
)

# Add cluster ID to each node (connected component)
graph_tbl <- graph_tbl %>%
  mutate(ClusterID = group_components())

# Extract cluster information
metadata_clusters <- as_tibble(graph_tbl) %>%
  select(Sample, Year, Location, ClusterID)

# Calculate clustering statistics per year
clustering_stats <- metadata_clusters %>%
  group_by(Year) %>%
  summarise(
    nb_samples = n(),
    nb_clusters = n_distinct(ClusterID),
    clustering_ratio = nb_clusters / nb_samples,
    .groups = "drop"
  ) %>%
  arrange(Year)

# Create cascade bar plot for clustering ratio
cascade_plot <- clustering_stats %>%
  mutate(clustering_ratio = clustering_ratio * 100) %>%
  ggplot(aes(x = interaction(Year), y = clustering_ratio)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = paste0(round(clustering_ratio), "%")), vjust = -0.5, size = 5) +
  labs(
    title = "IBD Clustering Rate by Year",
    subtitle = "(Infectious Diversity)",
    x = "Year - Location",
    y = "Clusters / Samples (%)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.background = element_rect(fill = "grey95", color = NA),
    panel.grid.major.y = element_line(color = "grey80", linetype = "dashed"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12, margin = margin(b = 10))
  )

# Prepare data for the summary table
table_data <- clustering_stats %>%
  mutate(clustering_ratio = round(clustering_ratio, 2))

# Create a tableGrob
stat_table <- tableGrob(
  table_data,
  rows = NULL,
  theme = ttheme_default(
    base_size = 12,
    core = list(fg_params = list(hjust = 0.5)),
    colhead = list(fg_params = list(fontface = "bold", hjust = 0.5))
  )
)

# Add title to the table
title_table <- ggdraw() +
  draw_label(
    "Summary of Statistics by Year and Location",
    fontface = "bold",
    size = 14,
    hjust = -0.1,
    x = 0
  )

# Combine title and table
stat_table_block <- cowplot::plot_grid(
  title_table,
  stat_table,
  ncol = 1,
  rel_heights = c(0.12, 1)
)

# Create cluster graph plot (color = Year, shape = Location)
graph_plot <- ggraph(graph_tbl, layout = "fr") +
  geom_edge_link(color = "black", alpha = 0.6, width = 0.6) +
  geom_node_point(aes(
    color = as.factor(Year),
    shape = as.factor(Location)
  ), size = 3, stroke = 1.2, fill = "white") +
  scale_color_manual(values = year_colors, guide = "none") +  # No color legend
  scale_shape_manual(values = location_shapes, guide = "none") +  # No shape legend
  theme_void(base_size = 14) +
  theme(
    plot.background = element_rect(fill = "grey85", linewidth = 0.5),
    plot.title = element_text(hjust = 0.2, size = 16, face = "bold")
  ) +
  ggtitle("IBD Clusters (Color = Year, Shape = Location)")

# Create separate legend plot
legend_plot <- ggplot(metadata, aes(x = Year, y = Location, color = as.factor(Year), shape = as.factor(Location))) +
  geom_point(size = 3, stroke = 1.2, fill = "white") +
  scale_color_manual(values = year_colors, name = "Year") +
  scale_shape_manual(values = location_shapes, name = "Location") +
  theme_void() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 26, face = "bold"),
    legend.text = element_text(size = 21)
  )

# Extract the legend
the_legend <- cowplot::get_legend(legend_plot)

# Assemble left column: cascade plot + table
left_col <- cowplot::plot_grid(
  cascade_plot,
  stat_table_block,
  ncol = 1,
  rel_heights = c(0.6, 0.4),
  labels = c("A", "C"),
  label_size = 20,
  label_fontface = "bold",
  label_x = 0,
  label_y = c(1, 1)
)

# Center column: graph
middle_col <- cowplot::plot_grid(
  graph_plot,
  labels = "B",
  label_size = 20,
  label_fontface = "bold",
  label_x = 0,
  label_y = 1
)

# Right column: legend
right_col <- cowplot::plot_grid(
  the_legend,
  ncol = 1
)

# Final assembly of all plots
final_plot <- cowplot::plot_grid(
  left_col,
  middle_col,
  right_col,
  ncol = 3,
  rel_widths = c(0.35, 0.5, 0.15)
)

# Add global title
final_plot_titled <- cowplot::plot_grid(
  ggdraw() + draw_label(
    paste0(
      "Clustering of R539T Samples based on Identity by Descent (IBD) Percentage,\nThreshold >= ", ibd_threshold
    ),
    fontface = "bold",
    size = 16,
    hjust = 0.5,
    lineheight = 1.1
  ),
  final_plot,
  ncol = 1,
  rel_heights = c(0.12, 1)
)

# Display the final plot
print(final_plot_titled)
