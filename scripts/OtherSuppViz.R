packages <- c("sf", "terra", "here", "rgrass", "dplyr", "data.table", 
              "rredlist", "tidyr", "Matrix", "tidyverse", "pROC",
              "kableExtra", "ggplot2", "reshape2", "parallel", "gridExtra")
lapply(packages, library, character.only = TRUE)

setwd('/')
setwd('Users/wenxinyang/Desktop/GitHub/colander')
source('scripts/refineBiomodelos/funcs.R')
source('scripts/refineBiomodelos/info.R')


## -------------------- Fig. S1 ---------------------
gen_number <- 5
df_tertile <- read.csv(sprintf('results/glm_btst_2012_keephab/gen%d_glm_btst_2012_keephab/btst_ideam_randomCI_gen%d__0615_2026_tertiles.csv',
                              gen_number, gen_number))


# Calculate maximum frequency across both histograms to set common y-axis range
h1 <- hist(df_tertile$one.third, breaks = 25, plot = FALSE)
h2 <- hist(df_tertile$two.third, breaks = 25, plot = FALSE)
max_freq <- max(c(max(h1$counts), max(h2$counts)))
# Add some padding (10% above max)
y_max <- ceiling(max_freq * 1.1)

# Create histogram for one.third
p_one_third <- ggplot(df_tertile, aes(x = one.third)) +
  geom_histogram(fill = "#4a90a4", color = "#2c5f6b", alpha = 0.8, bins = 30) +
  scale_y_continuous(limits = c(0, y_max)) +
  theme_minimal() +
  theme(
    text = element_text(family = "Arial"),
    axis.text = element_text(size = 12, family = "Arial"),
    axis.title = element_text(size = 14, face = "bold", family = "Arial"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, family = "Arial"),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  labs(
    #title = "Distribution of One-Third Tertile Values",
    x = "One-Third Tertile Values",
    y = "Frequency"
  )

# Create histogram for two.third
p_two_third <- ggplot(df_tertile, aes(x = two.third)) +
  geom_histogram(fill = "#4a90a4", color = "#2c5f6b", alpha = 0.8, bins = 30) +
  scale_y_continuous(limits = c(0, y_max)) +
  theme_minimal() +
  theme(
    text = element_text(family = "Arial"),
    axis.text = element_text(size = 12, family = "Arial"),
    axis.title = element_text(size = 14, face = "bold", family = "Arial"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, family = "Arial"),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  labs(
    #title = "Distribution of Two-Third Tertile Values",
    x = "Two-Third Tertile Values",
    y = "Frequency"
  )

# Display both histograms side by side
p_all <- grid.arrange(p_one_third, p_two_third, ncol = 2)
ggsave(sprintf('viz/gen%d_glm_btst_keephab_terthist.png', gen_number), p_all,
       height=6, width=10)

## -------------------- Fig. S2 ---------------------
gen_number <- 7

pos_data <- read.csv(sprintf('results/glm_btst_2012_rmhab/gen%d_glm_btst_2012_rmhab/btst_ideam_randomCI_gen%d_rmhab_0615_2026_pos.csv', 
                             gen_number, gen_number))
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

p2 <- ggplot2::ggplot(plot_data, aes(x = habitat, y = land_cover, fill = count_value)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_tile(data = subset(plot_data, color_category == "n_samples"),
            fill = "white", color = "white", linewidth = 0.5) +
  #geom_tile(data = subset(plot_data, count_value >= 850 & color_category != "n_samples"),
   #         color = "#006d2c", linewidth = 2, fill = NA) +
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

ggsave(sprintf('viz/gen%d_glm_btst_2012_keephab_mat.png', gen_number), p2,
       height = 10, width = 12)

## ------------------------ Fig. S3 & S5 ---------------------
df_btst7 <- read.csv('results/glm_btst_2012_keephab/gen7_glm_btst_2012_keephab/btst_ideam_randomCI_gen7__0615_2026_pos.csv')
#df_rmhab7 <- read.csv('results/glm_btst_2012_rmhab/gen7_glm_btst_2012_rmhab/btst_ideam_randomCI_gen7_rmhab_0615_2026_pos.csv')


df_subsamp7 <- read.csv('results/glm_subsamp_2012_keephab/gen7_glm_subsampling_2012_keephab/btst_ideam_randomCI_gen7__0615_2026_pos.csv')
pos_data = df_subsamp7

df_nobal7 <- read.csv('results/glm_btst_2012_keephab/gen7_glm_bootstrap_2012_keephab_imb/btst_ideam_randomCI_gen7__nobal_0701_2026_pos.csv')
#pos_data = df_nobal7

df_btst7 <- df_btst7 %>% select(-all_of(c('X', 'land_cover', 'auc', 'n_samples')))
#df_rmhab7 <- df_rmhab7 %>% select(-all_of(c('X', 'land_cover', 'auc', 'n_samples')))
#df_btst71 <- df_btst7 %>% select(all_of(colnames(df_rmhab7)))

df_subsamp7 <- df_subsamp7 %>% select(-all_of(c('X', 'land_cover', 'auc', 'n_samples')))

matbtst7 <- as.vector(as.matrix((df_btst7)))
#matrmhab7 <- as.vector(as.matrix((df_rmhab7)))
#matbtst71 <- as.vector(as.matrix(df_btst71))
matsubsamp7 <- as.vector(as.matrix(df_subsamp7))

pearson_cor <- cor(matbtst7, matsubsamp7, method='pearson')
#cor(matbtst71, matrmhab7, method='pearson')
pearson_cor

spearman_cor <- cor(matbtst7, matsubsamp7, method='spearman')
#cor(matbtst71, matrmhab7, method='spearman')
spearman_cor

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


p2 <- ggplot2::ggplot(plot_data, aes(x = habitat, y = land_cover, fill = count_value)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_tile(data = subset(plot_data, color_category == "n_samples"),
            fill = "white", color = "white", linewidth = 0.5) +
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

ggsave('viz/nobal7_matrix.png', p2, height=10, width=12)


## --------------------- boxplot response -----------------
df <- read.csv('results/validation/aoh_validation_0626.csv') %>% filter(threshold==850)


a <- df %>% select(species, taxa) %>% unique()
table(a$taxa)

df <- df %>% filter(!is.na(point_prevalence) & !is.na(model_prevalence))

library(tidyr)
library(ggplot2)

# Reshape data to long format
df_long <- df %>%
  pivot_longer(
    cols = c(point_prevalence, model_prevalence),
    names_to = "prevalence_type",
    values_to = "prevalence_value"
  )

# Create side-by-side boxplots
ggplot(df_long, aes(x = prevalence_type, y = prevalence_value, fill = prevalence_type)) +
  geom_boxplot() +
  scale_fill_manual(
    values = c("point_prevalence" = "#5996FF", "model_prevalence" = "#FF9E20"),
    labels = c("point_prevalence" = "Point Prevalence", 
               "model_prevalence" = "Model Prevalence")
  ) +
  facet_wrap(~taxa, scales = "free") +
  labs(
    x = "Prevalence Type",
    y = "Prevalence Value",
    fill = "Prevalence"
  ) +
  theme_minimal()


## -------------------- Fig. S6.1. ---------------------
gen_number <- 7

pos_data <- read.csv('results/btst_2012_lv1/gen7_glm_bootstrap_2012_keephab/btst_ideam_randomCI_gen7__0702_2026_pos.csv')
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

p2 <- ggplot2::ggplot(plot_data, aes(x = habitat, y = land_cover, fill = count_value)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_tile(data = subset(plot_data, color_category == "n_samples"),
            fill = "white", color = "white", linewidth = 0.5) +
  #geom_tile(data = subset(plot_data, count_value >= 850 & color_category != "n_samples"),
  #         color = "#006d2c", linewidth = 2, fill = NA) +
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

ggsave('viz/lv1_mat.png', p2, height = 8, width = 12)


## -------------------- Fig. S6.2. ---------------------
df_all <- read.csv('data/occ_pts/allinfo_ideam_coords_2012_0605.csv')
colnames(df_all)
median(table(df_all$nvl_3_n))
write.csv(as.data.frame(table(df_all$nvl_3_n)), 'viz/tables6.csv')
?write.csv()

