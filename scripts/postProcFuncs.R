# ================== Analysis and Visualization Functions =====================
# Create matrices and plots from the results

# Function to extract habitat associations and raw odds ratios from iteration results
extract_habitat_data <- function(iteration_results) {
  # Combine all habitat association data
  all_hab_data <- do.call(rbind, lapply(iteration_results, function(x) x$hab_assoc))
  
  # Combine all raw odds ratios data
  all_raw_odds <- do.call(rbind, lapply(iteration_results, function(x) x$raw_odds))
  
  # Extract land cover and habitat pairs
  land_covers <- unique(all_hab_data$land_cover)
  habitats <- unique(all_raw_odds$habitat)
  
  return(list(all_hab_data = all_hab_data, all_raw_odds = all_raw_odds, 
              land_covers = land_covers, habitats = habitats))
}

# Function to create count matrix (iterations above 1)
create_count_matrix <- function(all_hab_data, land_covers, habitats) {
  count_matrix <- matrix(0, nrow = length(land_covers), ncol = length(habitats))
  rownames(count_matrix) <- land_covers
  colnames(count_matrix) <- habitats
  
  for (lc in land_covers) {
    lc_data <- all_hab_data[all_hab_data$land_cover == lc, ]
    for (hab in habitats) {
      # Count how many iterations have values > 1 for this pair
      count_above_1 <- sum(grepl(hab, lc_data$pos_odds_habitats, fixed = TRUE), na.rm = TRUE)
      count_matrix[lc, hab] <- count_above_1
    }
  }
  
  return(count_matrix)
}


# Function to create 95% CI matrix (width)
create_ci_matrix <- function(all_raw_odds, land_covers, habitats) {
  ci_matrix <- matrix(NA, nrow = length(land_covers), ncol = length(habitats))
  rownames(ci_matrix) <- land_covers
  colnames(ci_matrix) <- habitats
  
  for (lc in land_covers) {
    for (hab in habitats) {
      # Get all odds ratios for this land cover-habitat pair
      pair_data <- all_raw_odds[all_raw_odds$land_cover == lc & all_raw_odds$habitat == hab, ]
      values <- pair_data$odds_ratio
      values <- values[!is.na(values)]
      
      if (length(values) > 0) {
        # Calculate 95% confidence interval
        ci <- quantile(values, c(0.025, 0.975))
        ci_matrix[lc, hab] <- ci[2] - ci[1]  # CI width
      }
    }
  }
  
  return(ci_matrix)
}




# Function to create 95% CI matrix (width)
create_ci_matrix <- function(all_raw_odds, land_covers, habitats) {
  ci_matrix <- matrix(NA, nrow = length(land_covers), ncol = length(habitats))
  rownames(ci_matrix) <- land_covers
  colnames(ci_matrix) <- habitats
  
  for (lc in land_covers) {
    for (hab in habitats) {
      # Get all odds ratios for this land cover-habitat pair
      pair_data <- all_raw_odds[all_raw_odds$land_cover == lc & all_raw_odds$habitat == hab, ]
      values <- pair_data$odds_ratio
      values <- values[!is.na(values)]
      
      if (length(values) > 0) {
        # Calculate 95% confidence interval
        ci <- quantile(values, c(0.025, 0.975))
        ci_matrix[lc, hab] <- ci[2] - ci[1]  # CI width
      }
    }
  }
  
  return(ci_matrix)
}

# Function to create 95% CI bounds table
create_ci_bounds_table <- function(all_raw_odds, land_covers, habitats) {
  ci_bounds_matrix <- matrix(NA, nrow = length(land_covers), ncol = length(habitats))
  rownames(ci_bounds_matrix) <- land_covers
  colnames(ci_bounds_matrix) <- habitats
  
  for (lc in land_covers) {
    for (hab in habitats) {
      # Get all odds ratios for this land cover-habitat pair
      pair_data <- all_raw_odds[all_raw_odds$land_cover == lc & all_raw_odds$habitat == hab, ]
      values <- pair_data$odds_ratio
      values <- values[!is.na(values)]
      
      if (length(values) > 0) {
        # Calculate 95% confidence interval bounds
        ci <- quantile(values, c(0.025, 0.975))
        lower_bound <- round(ci[1], 3)
        upper_bound <- round(ci[2], 3)
        ci_bounds_matrix[lc, hab] <- paste0("(", lower_bound, ", ", upper_bound, ")")
      } else {
        ci_bounds_matrix[lc, hab] <- "(NA, NA)"
      }
    }
  }
  
  return(ci_bounds_matrix)
}

# Function to create barplot of CI lengths
create_ci_barplot <- function(ci_matrix, dft_folder) {
  # Convert matrix to long format for plotting
  ci_df <- as.data.frame(ci_matrix)
  ci_df$land_cover <- rownames(ci_df)
  ci_long <- tidyr::gather(ci_df, habitat, ci_width, -land_cover)
  ci_long <- ci_long[!is.na(ci_long$ci_width), ]
  
  # Create barplot
  p <- ggplot2::ggplot(ci_long, aes(x = reorder(paste(land_cover, habitat), ci_width), y = ci_width)) +
    geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
    labs(title = "95% Confidence Interval Widths for Land Cover-Habitat Pairs",
         x = "Land Cover-Habitat Pair",
         y = "95% CI Width") +
    coord_flip()
  
  # Save plot
  ggplot2::ggsave(paste0(dft_folder, "/ci_width_barplot.png"), p, width = 12, height = 8, dpi = 300)
  
  return(p)
}

# Function to create heatmap of CI bounds with color coding
create_ci_bounds_heatmap <- function(all_raw_odds, land_covers, habitats, dft_folder) {
  # Create data frame for plotting
  plot_data <- data.frame()
  
  for (lc in land_covers) {
    for (hab in habitats) {
      # Get all odds ratios for this land cover-habitat pair
      pair_data <- all_raw_odds[all_raw_odds$land_cover == lc & all_raw_odds$habitat == hab, ]
      values <- pair_data$odds_ratio
      values <- values[!is.na(values)]
      
      if (length(values) > 0) {
        # Calculate 95% confidence interval bounds
        ci <- quantile(values, c(0.025, 0.975))
        lower_bound <- ci[1]
        upper_bound <- ci[2]
        mean_value <- mean(values)
        
        plot_data <- rbind(plot_data, data.frame(
          land_cover = lc,
          habitat = hab,
          lower_bound = lower_bound,
          upper_bound = upper_bound,
          mean_value = mean_value,
          ci_text = paste0("(", round(lower_bound, 3), ", ", round(upper_bound, 3), ")"),
          color_group = ifelse(lower_bound > 1, "Significant", "Not Significant")
        ))
      } else {
        plot_data <- rbind(plot_data, data.frame(
          land_cover = lc,
          habitat = hab,
          lower_bound = NA,
          upper_bound = NA,
          mean_value = NA,
          ci_text = "(NA, NA)",
          color_group = "No Data"
        ))
      }
    }
  }
  
  # Create heatmap
  p <- ggplot2::ggplot(plot_data, aes(x = habitat, y = land_cover)) +
    geom_tile(fill = "white", color = "grey80", linewidth = 0.5) +
    geom_text(aes(label = ci_text, color = color_group), 
              size = 3.5, fontface = "bold") +
    scale_color_manual(
      values = c("Significant" = "red", "Not Significant" = "black", "No Data" = "grey50"),
      name = "Significance"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 11),
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", color = NA)
    ) +
    labs(
      title = "95% Confidence Intervals for Land Cover-Habitat Pairs",
      x = "Habitat Type",
      y = "Land Cover Type"
    ) +
    coord_fixed(ratio = 0.6)  # Make cells even wider by reducing the ratio further
  
  # Save plot
  ggplot2::ggsave(paste0(dft_folder, "/ci_bounds_heatmap.png"), p, width = 16, height = 14, dpi = 300)
  
  return(p)
}

# Function to create heatmap using the entire _pos.csv file
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
