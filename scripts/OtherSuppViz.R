packages <- c("sf", "terra", "here", "rgrass", "dplyr", "data.table", 
              "rredlist", "tidyr", "Matrix", "tidyverse", "pROC",
              "kableExtra", "ggplot2", "reshape2", "parallel")
lapply(packages, library, character.only = TRUE)

setwd('/')
setwd('Users/wenxinyang/Desktop/GitHub/colander')
source('scripts/refineBiomodelos/funcs.R')
source('scripts/refineBiomodelos/info.R')

pos_data <- read.csv('results/glm_btst_2012_keephab/gen7_glm_btst_2012_keephab/btst_ideam_randomCI_gen7__0615_2026_pos.csv')
good_col_names <- gsub('\\.', ' ', colnames(pos_data))
good_col_names <- gsub("Wetlands  inland", "Wetland (inland)", good_col_names)
colnames(pos_data) <- good_col_names

if("X" %in% colnames(pos_data)){
  pos_data$X <- NULL
}

# Get habitat columns (exclude land_cover, auc, n_samples)
habitat_cols <- colnames(pos_data)[!colnames(pos_data) %in% c("", "land_cover", "auc", "n_samples")]

# Convert to long format for plotting
pos_long <- pos_data %>%
  select(land_cover, all_of(habitat_cols)) %>%
  tidyr::gather(habitat, count_value, -land_cover)

# Create color categories for habitat columns
pos_long$color_category <- cut(pos_long$count_value, 
                               breaks = c(-Inf, 850, Inf),
                               labels = c("Not a pair (<850)", "AOH (≥850)"),
                               include.lowest = TRUE)

# Add AUC and n_samples data
auc_data <- data.frame(
  land_cover = pos_data$land_cover,
  habitat = "AUC",
  count_value = pos_data$auc,
  color_category = "AUC"
)

n_samples_data <- data.frame(
  land_cover = pos_data$land_cover,
  habitat = "n_samples",
  count_value = round(pos_data$n_samples),
  color_category = "n_samples"
)

# Combine all data
plot_data <- rbind(pos_long, auc_data, n_samples_data)

# Reorder color categories for legend: Low at top, Not a pair at bottom
plot_data$color_category <- factor(plot_data$color_category,
                                   levels = c("AOH (≥850)", "Not a pair (<850)", "AUC", "n_samples"))

# Order land cover types alphabetically in descending order (Z to A)
plot_data$land_cover <- factor(plot_data$land_cover, 
                               levels = sort(unique(plot_data$land_cover), decreasing = TRUE))

# Set factor levels for habitat to include AUC and n_samples at the end
habitat_levels <- c(habitat_cols, "AUC", "n_samples")
plot_data$habitat <- factor(plot_data$habitat, levels = habitat_levels)

# Create heatmap
p <- ggplot2::ggplot(plot_data, aes(x = habitat, y = land_cover, fill = color_category)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(color_category == "AUC", 
                               sprintf("%.2f", count_value),
                               ifelse(color_category == "n_samples",
                                      as.character(count_value),
                                      as.character(count_value)))), 
            size = 5, fontface = "bold", color = "black", family = "Arial") +
  scale_fill_manual(
    values = c("Not a pair (<850)" = "#f7f7f7", "AOH (≥850)" = "#659c6a",
               "AUC" = "white", "n_samples" = "white"),
    name = "Certainty level (counts)"
  ) +
  scale_x_discrete(position = "bottom") +
  theme_minimal() +
  theme(
    text = element_text(family = "Arial"),
    axis.text.x = element_text(size = 12, family = "Arial", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 15, family = "Arial"),
    axis.title.x = element_text(size = 16, face = 'bold', family = "Arial"),
    axis.title.y = element_text(size = 16, face = 'bold', family = "Arial"),
    legend.title = element_text(size = 16, family = "Arial"),
    legend.text = element_text(size = 14, family = "Arial"),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold", family = "Arial"),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  labs(
    # title = "Habitat-land cover translation matrix",
    x = "Habitat Classes",
    y = "Land Cover Classes"
  ) +
  coord_fixed(ratio = 0.6)
p


p2 <- ggplot2::ggplot(plot_data, aes(x = habitat, y = land_cover, fill = count_value)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_tile(data = subset(plot_data, color_category == "n_samples"),
            fill = "white", color = "white", linewidth = 0.5) +
  geom_tile(data = subset(plot_data, count_value >= 850 & color_category != "n_samples"),
            color = "#006d2c", linewidth = 2, fill = NA) +
  geom_text(aes(label = ifelse(color_category == "AUC", 
                               sprintf("%.2f", count_value),
                               ifelse(color_category == "n_samples",
                                      as.character(count_value),
                                      as.character(count_value)))), 
            size = 5, fontface = "bold", color = "black", family = "Arial") +
  scale_fill_gradient(
    low = "white",
    high = "#66c2a4",
    na.value = "grey90",
    name = "Counts",
    limits = c(0, 1000),
    breaks = seq(0, 1000, 200)
  ) +
  guides(
    fill = guide_colorbar(
      title = "Counts",
      title.position = "top",
      title.hjust = 0.5,
      barwidth = 1.5,
      barheight = 15,
      title.theme = element_text(
        size = 16, 
        family = "Arial",
        margin = margin(b = 15)
      ),
      label.position = "right",  # Move labels to the right of the color bar
      label.hjust = 0.5,         # Horizontal justification
      label.vjust = 0.5          # Vertical justification
    )
  ) +
  scale_x_discrete(position = "bottom") +
  theme_minimal() +
  theme(
    text = element_text(family = "Arial"),
    axis.text.x = element_text(size = 12, family = "Arial", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 15, family = "Arial"),
    axis.title.x = element_text(size = 16, face = 'bold', family = "Arial"),
    axis.title.y = element_text(size = 16, face = 'bold', family = "Arial"),
    legend.title = element_text(
      size = 16, 
      family = "Arial",
      face = "bold",
      margin = margin(b = 10)
    ),
    legend.text = element_text(
      size = 14, 
      family = "Arial",
      margin = margin(l = 5, t = 5, b = 5),  # Added left margin to push text right
      hjust = 0  # Left-align the text
    ),
    legend.position = "right",
    legend.spacing.y = unit(0.5, "cm"),
    legend.box.spacing = unit(0.5, "cm"),
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold", family = "Arial"),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  labs(
    x = "Habitat Classes",
    y = "Land Cover Classes"
  ) +
  coord_fixed(ratio = 0.6)
p2