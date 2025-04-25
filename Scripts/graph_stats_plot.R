library(tidygraph)
library(ggraph)
library(ggplot2)
library(readr)
library(RColorBrewer)
library(cowplot)
library(dplyr)
library(gridExtra)
library(grid)

setwd("/Volumes/ANNA_DD/New_data_Dak_Nong/hmmIBD")

# Lire les données
ibd_data <- read_tsv("IBD_fract.txt")
metadata <- read_tsv("metadata_mutations.txt")

# Nettoyer les noms d'échantillons
metadata$Sample <- tolower(trimws(metadata$Sample))
ibd_data$sample1 <- tolower(trimws(ibd_data$sample1))
ibd_data$sample2 <- tolower(trimws(ibd_data$sample2))

# Filtrage IBD
ibd_threshold <- 0.55
filtered_ibd <- ibd_data %>% filter(IBD_value >= ibd_threshold)

# Sécurité : ne garder que les arêtes valides
filtered_ibd <- filtered_ibd %>%
  filter(sample1 %in% metadata$Sample & sample2 %in% metadata$Sample)

# Graphe tidygraph avec tous les nœuds
graph_tbl <- tbl_graph(
  nodes = metadata,
  edges = filtered_ibd[, c("sample1", "sample2")],
  directed = FALSE
)

# Couleurs & formes
year_colors <- brewer.pal(length(unique(metadata$Year)), "Set1")
names(year_colors) <- sort(unique(metadata$Year))
mutation_shapes <- c(
  "WT" = 21, "C580Y" = 22, "R539T" = 23,
  "C469F" = 24, "Y511H" = 25, "P553L" = 8
)

# Ajouter l'ID de cluster
graph_tbl <- graph_tbl %>%
  mutate(ClusterID = group_components())

# Extraire infos pour stats
metadata_clusters <- as_tibble(graph_tbl) %>%
  select(Sample, Year, Mutation, ClusterID)

# Statistiques de clustérisation
clustering_stats <- metadata_clusters %>%
  group_by(Year) %>%
  summarise(
    nb_samples = n(),
    nb_clusters = n_distinct(ClusterID),
    clustering_ratio = nb_clusters / nb_samples
  ) %>%
  arrange(Year)

# Cascade plot
cascade_plot <- clustering_stats %>%
  mutate(clustering_ratio = clustering_ratio * 100) %>%
  ggplot(aes(x = factor(Year), y = clustering_ratio)) +
  geom_col(fill = "#69b3a2", width = 0.7) +
  geom_text(aes(label = paste0(round(clustering_ratio), "%")), vjust = -0.5, size = 5) +
  labs(
    title = "Taux de clustérisation IBD par année",
    subtitle = "(Diversité infectieuse)",
    x = "Année",
    y = "Clusters / Échantillons (%)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.background = element_rect(fill = "grey95", color = NA),
    panel.grid.major.y = element_line(color = "grey80", linetype = "dashed"),
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12, margin = margin(b = 10))
  )

# Mise en forme du tableau
table_data <- clustering_stats %>%
  mutate(clustering_ratio = round(clustering_ratio * 100, 2))

stat_table <- tableGrob(
  table_data,
  rows = NULL,
  theme = ttheme_default(
    base_size = 15,
    core = list(fg_params = list(hjust = 0.5)),
    colhead = list(fg_params = list(fontface = "bold", hjust = 0.5))
  )
)

title_table <- ggdraw() +
  draw_label(
    "Résumé des statistiques par année",
    fontface = "bold",
    size = 16,
    hjust = -0.2,
    x = 0
  )

stat_table_block <- cowplot::plot_grid(
  title_table,
  stat_table,
  ncol = 1,
  rel_heights = c(0.12, 1)
)

# Graphe sans légende
graph_plot <- ggraph(graph_tbl, layout = "fr") +
  geom_edge_link(color = "black", alpha = 0.6, width = 0.6) +
  geom_node_point(aes(
    shape = Mutation,
    color = as.factor(Year)
  ), size = 3, stroke = 1.2, fill = "white") +
  scale_color_manual(values = year_colors, guide = "none") +
  scale_shape_manual(values = mutation_shapes, guide = "none") +
  theme_void(base_size = 14) +
  theme(
    plot.background = element_rect(fill = "grey85", linewidth = 0.5),
    plot.title = element_text(hjust = 0.07, size = 16, face = "bold")
  ) +
  ggtitle("Réprésentation des clusters")

# Légende
legend_plot <- ggplot(metadata, aes(x = Year, y = Mutation, color = as.factor(Year), shape = Mutation)) +
  geom_point(size = 3, fill = "white", stroke = 1.2) +
  scale_color_manual(values = year_colors, name = "Année") +
  scale_shape_manual(values = mutation_shapes, name = "Mutation") +
  theme_void() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 26, face = "bold"),
    legend.text = element_text(size = 21)
  )

the_legend <- cowplot::get_legend(legend_plot)

# Colonne gauche
left_col <- cowplot::plot_grid(
  cascade_plot,
  stat_table_block,
  ncol = 1,
  rel_heights = c(0.6, 0.4),
  labels = c("A", "C"),
  label_size = 20,
  label_fontface = "bold"
)

# Centre
middle_col <- cowplot::plot_grid(
  graph_plot,
  labels = "B",
  label_size = 20,
  label_fontface = "bold"
)

# Droite
right_col <- cowplot::plot_grid(
  the_legend,
  ncol = 1
)

# Assemblage final
final_plot <- cowplot::plot_grid(
  left_col,
  middle_col,
  right_col,
  ncol = 3,
  rel_widths = c(0.35, 0.5, 0.15)
)

final_plot_titled <- cowplot::plot_grid(
  ggdraw() + draw_label(
    "Clustérisation des échantillons basée sur leur pourcentage d'Identity by Descent (IBD)\ndans la région de Dak Nong (Vietnam), Seuil >= 0.55",
    fontface = "bold",
    size = 16,
    hjust = 0.5,
    lineheight = 1.1
  ),
  final_plot,
  ncol = 1,
  rel_heights = c(0.12, 1)
)

print(final_plot_titled)