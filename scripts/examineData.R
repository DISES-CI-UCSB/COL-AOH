# Examine data structure and create heatmap of land cover x habitat counts
# Author: Wenxin Yang
# Date: January, 2025

# ================== Load libraries =====================
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(viridis)
library(kableExtra)
library(reshape2)

# ================== Set up working directory =====================
setwd('/')
setwd('Users/wenxinyang/Desktop/GitHub/colander')
source('scripts/refineBiomodelos/funcs.R')
source('scripts/refineBiomodelos/info.R')

# ================== Read and prepare data =====================
# Read the main data file
df_all <- read.csv('data/occ_pts/allinfo_ideam_cgls_coords.csv')

# Read species preference data
df_pref <- read.csv('data/occ_pts/col_animal_pref_cleaned.csv')

# Add generalist info (using default 5 generalists)
df_basic_info <- addGeneralistInfo(df_pref, num_generalist = 5)

# Remove habitat codes not included in the analysis
df_all_info <- df_all %>% select(-all_of(drop_cols))

# Remove those near boundary (using 0 as default)
d_near <- 0
df_all_info <- df_all_info %>% filter(dst_t_b > d_near) %>% select(-dst_t_b)

# Merge several artificial habitat types
df_all_info <- df_all_info %>% mutate(
  hab_14.12 = ifelse(hab_14.1+hab_14.2 >0, 1, 0),
  hab_14.36 = ifelse(hab_14.3+hab_14.6 >0, 1, 0),
  hab_14.45 = ifelse(hab_14.4+hab_14.5 >0, 1, 0)
)

# Remove the original ones
df_all_info <- df_all_info %>% select(-all_of(c('hab_14.1', 'hab_14.2', 'hab_14.3', 'hab_14.4', 'hab_14.5', 'hab_14.6')))

# Remove rows with NA values
df_all_info <- df_all_info[complete.cases(df_all_info), ]
df_all_info$X <- NULL

# Sum all columns that start with "hab_"
df_all_info$sum <- rowSums(df_all_info[,grep("^hab_", colnames(df_all_info))])
df_all_info <- df_all_info %>% filter(sum>0) %>% select(-sum)

# Remove generalist species
li_generalist <- unique((df_basic_info %>% filter(type == 'generalist'))$name)
df_all_info <- df_all_info %>% filter(!species %in% li_generalist)

# Merge with basic info
new_rec <- merge(df_all_info, df_basic_info, by.x='species', by.y='name')

cat("Total records:", nrow(new_rec), "\n")
cat("Total species:", length(unique(new_rec$species)), "\n")
cat("Total land cover types:", length(unique(new_rec$nvl_2_n)), "\n")

# ================== Create land cover x habitat count matrix =====================
# Get habitat columns
habitat_cols <- colnames(new_rec)[grep("^hab_", colnames(new_rec))]

# Create count matrix
count_matrix <- matrix(0, nrow = length(unique(new_rec$nvl_2_n)), ncol = length(habitat_cols))
rownames(count_matrix) <- sort(unique(new_rec$nvl_2_n))
colnames(count_matrix) <- habitat_cols

# Fill the matrix with counts
for (lc in rownames(count_matrix)) {
  for (hab in colnames(count_matrix)) {
    # Count records where this land cover has this habitat
    count <- sum(new_rec$nvl_2_n == lc & new_rec[[hab]] == 1)
    count_matrix[lc, hab] <- count
  }
}

# Get land cover names for better labels
lc_names <- sapply(rownames(count_matrix), function(x) {
  ideam_lc_info$ideam_lc_name[match(paste0('lc_',x), ideam_lc_info$ideam_lc_code)]
})

# Get habitat names for better labels
hab_names <- sapply(colnames(count_matrix), function(x) {
  habitat_info1$habitat_name[match(x, habitat_info1$habitat_code)]
})

# Create labeled matrix for display
labeled_matrix <- count_matrix
rownames(labeled_matrix) <- lc_names  # Use only names, not codes
colnames(labeled_matrix) <- hab_names  # Use only names, not codes

# ================== Create heatmap similar to the image =====================
# Create heatmap with counts displayed
heatmap_plot <- pheatmap(
  labeled_matrix,
  display_numbers = TRUE,  # Show counts in cells
  number_format = "%.0f",  # No decimal places
  fontsize_number = 9,    # Size of numbers
  fontsize_row = 10,       # Size of row labels
  fontsize_col = 10,       # Size of column labels
  color = colorRampPalette(c("white", "lightgreen", "darkgreen"))(100),  # Green gradient like the image
  cluster_rows = FALSE,    # Don't cluster rows to maintain order
  cluster_cols = FALSE,    # Don't cluster columns to maintain order
  main = "Land Cover x Habitat Count Matrix",
  angle_col = 45,          # Angle column labels
  width = 12,              # Width in inches
  height = 10,             # Height in inches
  border_color = "white", # Border color for cells
  cellwidth = 20,          # Cell width
  cellheight = 15          # Cell height
)




# Normalize the count matrix by row to get percentages
row_sums <- rowSums(count_matrix)
# Avoid division by zero
row_sums[row_sums == 0] <- 1
percent_matrix <- sweep(count_matrix, 1, row_sums, FUN = "/") * 100

# Round percentages for display
percent_matrix_rounded <- round(percent_matrix, 0)

# Create labeled matrix for display (percentages)
labeled_percent_matrix <- percent_matrix_rounded
rownames(labeled_percent_matrix) <- lc_names
colnames(labeled_percent_matrix) <- hab_names

# Create heatmap with percentages displayed
heatmap_percent_plot <- pheatmap(
  labeled_percent_matrix,
  display_numbers = TRUE,  # Show percentages in cells
  number_format = "%.0f",  # No decimal places
  fontsize_number = 9,
  fontsize_row = 10,
  fontsize_col = 10,
  color = colorRampPalette(c("white", "lightblue", "blue"))(100),  # Blue gradient for percent
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  main = "Land Cover x Habitat Percentage Matrix",
  angle_col = 45,
  width = 12,
  height = 10,
  border_color = "white",
  cellwidth = 20,
  cellheight = 15
)

# ================= test if rows are independent =================
# Test if rows of the count matrix are independent (chi-squared test of independence)
cat("\n=== TESTING ROW INDEPENDENCE (Chi-squared test) ===\n")
chisq_test <- chisq.test(count_matrix)

cat("Chi-squared test statistic:", chisq_test$statistic, "\n")
cat("Degrees of freedom:", chisq_test$parameter, "\n")
cat("P-value:", chisq_test$p.value, "\n")

if (chisq_test$p.value < 0.05) {
  cat("Result: Rows and columns are NOT independent (p < 0.05)\n")
} else {
  cat("Result: No evidence against independence (p >= 0.05)\n")
}

# ================== Find habitats without specialist =====================
colnames(new_rec)

spp <- df_basic_info %>% filter(type == 'specialist')
spp <- unique(spp$name) # 1048 specialist species
df_spp <- new_rec %>% filter(species %in% spp)
hab_cols <- grep("^hab_", colnames(df_spp), value = TRUE)

df_spp$specialist_habitat <- apply(df_spp[, hab_cols], 1, function(x) {
  habs <- hab_cols[which(x == 1)]
  if (length(habs) == 0) {
    return(NA)
  } else if (length(habs) == 1) {
    return(habs)
  } else {
    return(paste(habs, collapse = ";"))
  }
})

unique(df_spp$specialist_habitat)


df_allsp <- new_rec
df_allsp$specialist_habitat <- apply(df_allsp[, hab_cols], 1, function(x) {
  habs <- hab_cols[which(x == 1)]
  if (length(habs) == 0) {
    return(NA)
  } else if (length(habs) == 1) {
    return(habs)
  } else {
    return(paste(habs, collapse = ";"))
  }
})

li_allcombo <- unique(df_allsp$specialist_habitat)
hab_15_or_1412 <- li_allcombo[
  grepl("hab_15", li_allcombo) | grepl("hab_14.12", li_allcombo)]
cat("Habitats containing 'hab_15' or 'hab_14.12':\n")
print(hab_15_or_1412)


count_habs <- as.data.frame(table(df_allsp$specialist_habitat))
count_habs$num_habs <- sapply(as.character(count_habs$Var1), function(x) {
  if (is.na(x) || x == "NA") {
    return(0)
  } else {
    return(length(strsplit(x, ";")[[1]]))
  }
})

## ------ find habitat combo without 1, 3, 4, 5 ---------
no_sp <- count_habs %>% filter(!grepl('hab_1|hab_3|hab_4|hab_5', Var1))



# ================== Create summary statistics =====================
cat("Creating summary statistics...\n")

# Summary by land cover
lc_summary <- data.frame(
  land_cover_code = rownames(count_matrix),
  land_cover_name = lc_names,
  total_records = rowSums(count_matrix),
  habitats_with_records = rowSums(count_matrix > 0),
  max_habitat_count = apply(count_matrix, 1, max),
  min_habitat_count = apply(count_matrix, 1, function(x) min(x[x > 0]))
)

# Summary by habitat
hab_summary <- data.frame(
  habitat_code = colnames(count_matrix),
  habitat_name = hab_names,
  total_records = colSums(count_matrix),
  land_covers_with_records = colSums(count_matrix > 0),
  max_land_cover_count = apply(count_matrix, 2, max),
  min_land_cover_count = apply(count_matrix, 2, function(x) min(x[x > 0]))
)

# Print summaries
cat("\n=== LAND COVER SUMMARY ===\n")
print(lc_summary)

cat("\n=== HABITAT SUMMARY ===\n")
print(hab_summary)

# ================== Create data quality checks =====================
cat("\n=== DATA QUALITY CHECKS ===\n")

# Check for land covers with very few records
low_count_lc <- lc_summary[lc_summary$total_records < 50, ]
if (nrow(low_count_lc) > 0) {
  cat("Land covers with < 50 records:\n")
  print(low_count_lc)
} else {
  cat("All land covers have >= 50 records\n")
}

# Check for habitats with very few records
low_count_hab <- hab_summary[hab_summary$total_records < 50, ]
if (nrow(low_count_hab) > 0) {
  cat("Habitats with < 50 records:\n")
  print(low_count_hab)
} else {
  cat("All habitats have >= 50 records\n")
}

# Check for zero-variance habitats (all 0s or all 1s)
zero_var_habitats <- c()
for (hab in habitat_cols) {
  if (length(unique(new_rec[[hab]])) == 1) {
    zero_var_habitats <- c(zero_var_habitats, hab)
  }
}

if (length(zero_var_habitats) > 0) {
  cat("Zero-variance habitats (all same value):\n")
  print(zero_var_habitats)
} else {
  cat("All habitats have variation\n")
}

# ================== Save results =====================
# Save count matrix
write.csv(count_matrix, "data/land_cover_habitat_counts.csv")

# Save labeled matrix
write.csv(labeled_matrix, "data/land_cover_habitat_counts_labeled.csv")

# Save summaries
write.csv(lc_summary, "data/land_cover_summary.csv")
write.csv(hab_summary, "data/habitat_summary.csv")

cat("\n=== FILES SAVED ===\n")
cat("land_cover_habitat_counts.csv - Raw count matrix\n")
cat("land_cover_habitat_counts_labeled.csv - Labeled count matrix\n")
cat("land_cover_summary.csv - Land cover summary statistics\n")
cat("habitat_summary.csv - Habitat summary statistics\n")
cat("land_cover_habitat_heatmap.png - Heatmap visualization\n")

# ================== Create additional visualizations =====================
cat("\nCreating additional visualizations...\n")

# 1. Bar plot of total records by land cover
lc_bar <- ggplot(lc_summary, aes(x = reorder(land_cover_name, total_records), y = total_records)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Total Records by Land Cover Type",
       x = "Land Cover Type",
       y = "Number of Records") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))

ggsave("figures/land_cover_records_bar.png", lc_bar, width = 10, height = 8, dpi = 300)

# 2. Bar plot of total records by habitat
hab_bar <- ggplot(hab_summary, aes(x = reorder(habitat_name, total_records), y = total_records)) +
  geom_bar(stat = "identity", fill = "darkgreen") +
  coord_flip() +
  labs(title = "Total Records by Habitat Type",
       x = "Habitat Type",
       y = "Number of Records") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))

ggsave("figures/habitat_records_bar.png", hab_bar, width = 10, height = 8, dpi = 300)

# 3. Correlation heatmap of habitats
habitat_cor <- cor(new_rec[, habitat_cols])
habitat_cor_plot <- pheatmap(
  habitat_cor,
  display_numbers = TRUE,
  number_format = "%.2f",
  fontsize_number = 6,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "Habitat Correlation Matrix",
  width = 12,
  height = 10
)

ggsave("figures/habitat_correlation_heatmap.png", 
       plot = habitat_cor_plot, 
       width = 12, height = 10, 
       dpi = 300)

cat("Additional visualizations saved to figures/ directory\n")

# ================== Final summary =====================
cat("\n=== FINAL SUMMARY ===\n")
cat("Data examination complete!\n")
cat("Total records analyzed:", nrow(new_rec), "\n")
cat("Land cover types:", nrow(count_matrix), "\n")
cat("Habitat types:", ncol(count_matrix), "\n")
cat("Non-zero entries in matrix:", sum(count_matrix > 0), "\n")
cat("Matrix sparsity:", round((1 - sum(count_matrix > 0) / (nrow(count_matrix) * ncol(count_matrix))) * 100, 1), "%\n")
cat("Maximum count in matrix:", max(count_matrix), "\n")
cat("Minimum non-zero count in matrix:", min(count_matrix[count_matrix > 0]), "\n")
