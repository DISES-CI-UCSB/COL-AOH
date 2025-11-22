packages <- c("sf", "terra", "here", "rgrass", "dplyr", "data.table", 
              "rredlist", "tidyr", "Matrix", "tidyverse", "pROC",
              "kableExtra", "ggplot2", "reshape2", "parallel", "gridExtra")
lapply(packages, library, character.only = TRUE)

setwd('/')
setwd('Users/wenxinyang/Desktop/GitHub/colander')
source('scripts/refineBiomodelos/funcs.R')
source('scripts/refineBiomodelos/info.R')
# ============== 1. sensitivity of tertiles =============
dft_folder <- 'results/aoh_results_randomCI_final_2022/'
tertile_file <- list.files(dft_folder, pattern = "_tertiles.csv", full.names = TRUE)
df_tertile <- read.csv(tertile_file)

# Calculate maximum frequency across both histograms to set common y-axis range
h1 <- hist(df_tertile$one.third, breaks = 30, plot = FALSE)
h2 <- hist(df_tertile$two.third, breaks = 30, plot = FALSE)
max_freq <- max(c(max(h1$counts), max(h2$counts)))
# Add some padding (10% above max)
y_max <- ceiling(max_freq * 1.1)

# Create histogram for one.third
p_one_third <- ggplot(df_tertile, aes(x = one.third)) +
  geom_histogram(fill = "#4a90a4", color = "#2c5f6b", alpha = 0.8, bins = 30) +
  scale_y_continuous(limits = c(0, y_max)) +
  theme_minimal() +
  theme(
    text = element_text(family = "Times New Roman"),
    axis.text = element_text(size = 12, family = "Times New Roman"),
    axis.title = element_text(size = 14, face = "bold", family = "Times New Roman"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, family = "Times New Roman"),
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
    text = element_text(family = "Times New Roman"),
    axis.text = element_text(size = 12, family = "Times New Roman"),
    axis.title = element_text(size = 14, face = "bold", family = "Times New Roman"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, family = "Times New Roman"),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  labs(
    #title = "Distribution of Two-Third Tertile Values",
    x = "Two-Third Tertile Values",
    y = "Frequency"
  )

# Display both histograms side by side
grid.arrange(p_one_third, p_two_third, ncol = 2)

# ============== 2. sensitivity of generalist species ============
# Function to create heatmap using the entire _pos.csv file
create_count_heatmap <- function(count_matrix, dft_folder) {
  # Find and read the _pos.csv file
  pos_file <- list.files(dft_folder, pattern = "_pos.csv", full.names = TRUE)
  if (length(pos_file) == 0) {
    stop("No _pos.csv file found in the specified folder")
  }
  
  pos_data <- read.csv(pos_file[1])
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
                                 breaks = c(-Inf, 300, 600, 900, Inf),
                                 labels = c("Not a pair (≤300)", "Low (300-600)", "Medium (600-900)", "High (>900)"),
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
                                     levels = c("High (>900)", "Medium (600-900)", "Low (300-600)", "Not a pair (≤300)", "AUC", "n_samples"))
  
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
      values = c("Not a pair (≤300)" = "#f7f7f7", "Low (300-600)" = "#00b3ff", "Medium (600-900)" = "#dbb300", "High (>900)" = "#659c6a",
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
  
  return(p)
}



dft_folder <- 'results/aoh_results_randomCI_final_7gen_2022/'
count_matrix <- read.csv(paste0(dft_folder, "/count_above_1_matrix.csv"))
p1 <- count_heatmap <- create_count_heatmap(count_matrix, dft_folder)
p1

