# This script visualizes the habitat land cover matrix
# Author: Wenxin Yang
# Date: July, 2025

# read in packages
library(reshape2)
library(dplyr)
library(ggplot2)

dft_folder <- "aoh_results_balanced/"
# read in all_results.RData
load(paste0(dft_folder, "/all_results.RData"))

test <- all_results[[1]]$pos

# if all_results is not saved
#all_results <- list()
#for(i in 1:10){
#  path_pos <- file.path(dft_folder, paste0('btst_ideam_', as.character(i), '_gen5_boot1000_0708_2025.csv'))
#  ideam <- read.csv(path_pos)
#  pos <- selectPairs(ideam, 'positive') 
#  path_tert <- file.path(dft_folder, paste0('btst_ideam_', as.character(i), '_gen5_boot1000_0708_2025_tertiles.csv'))
#  tert <- read.csv(path_tert)

#  all_results[[i]] <- list(pos = pos, tert = tert)

#}

# ================== Process results from all_results =====================
## ======= output 1 heatmaps =======
# for all elements in all_results, get pos
# output 1: for all elements in the list, get pos, and generate a data frame
# columns: land cover, habitat, number of runs where the value is > 500, number of runs where the value is > 800

# extract positive associations from all iterations
output1_data <- data.frame()

for(iter in seq_along(all_results)) {
  # get the positive associations for this iteration
  pos_data <- all_results[[iter]]$pos
  
  # for each land cover type
  for(lc in unique(pos_data$land_cover)) {
    lc_data <- pos_data[pos_data$land_cover == lc, ]
    
    # for each habitat column (excluding land_cover, auc, n_samples)
    habitat_cols <- colnames(lc_data)[!colnames(lc_data) %in% c("land_cover", "auc", "n_samples")]
    
    for(hab in habitat_cols) {
      # get the count for this habitat-land cover pair
      count_val <- lc_data[[hab]][1]  # Should be the same for all rows of same land cover
      
      # add to output1_data
      output1_data <- rbind(output1_data, data.frame(
        land_cover = lc,
        habitat = hab,
        iteration = iter,
        positive_runs = count_val
      ))
    }
  }
}

# create the new output1 with threshold counts
# get unique land covers and habitats
unique_land_covers <- unique(output1_data$land_cover)
unique_habitats <- unique(output1_data$habitat)

# create output1 with threshold counts
output1_thresholds <- data.frame()

for(lc in unique_land_covers) {
  for(hab in unique_habitats) {
    # get all iterations for this habitat-land cover pair
    pair_data <- output1_data[output1_data$habitat == hab & output1_data$land_cover == lc, ]
    
    if(nrow(pair_data) > 0) {
      # count runs above thresholds
      runs_above_500 <- sum(pair_data$positive_runs > 500)
      runs_above_800 <- sum(pair_data$positive_runs > 800)
      
      # add to output1_thresholds
      output1_thresholds <- rbind(output1_thresholds, data.frame(
        land_cover = lc,
        habitat = hab,
        runs_above_500 = runs_above_500,
        runs_above_800 = runs_above_800
      ))
    }
  }
}

#write.csv(output1_data, paste0(dft_folder, "/output1_detailed_positive_associations.csv"), row.names = FALSE)
#write.csv(output1_thresholds, paste0(dft_folder, "/output1_threshold_counts.csv"), row.names = FALSE)

# create two matrices to show land cover, habitat, # of runs for 800 and 500 respectively


# Define custom order for habitats
custom_habitat_order <- c(
  "Artificial-arable and pasture",
  "Artificial-degraded forest and plantation", 
  "Artificial-urban areas and rural gardens",
  "Desert",
  "Grassland",
  "Savanna",
  "Shrubland", 
  "Wetlands (inland)",
  "Artificial-aquatic",
  "Forest",
  "Rocky areas"
)

# Sort unique habitats by custom order and land covers alphabetically

# alphabetical order
unique_habitats_sorted <- sort(unique_habitats)
unique_land_covers_sorted <- sort(unique_land_covers)
# custom order
#unique_habitats_sorted <- unique_habitats[order(match(unique_habitats, custom_habitat_order))]



# for runs above 500 threshold
# colnames: land cover
# rownames: habitat
# cell values: number of runs where positive_runs > 500
matrix_500 <- matrix(0, nrow = length(unique_habitats_sorted), ncol = length(unique_land_covers_sorted))
rownames(matrix_500) <- unique_habitats_sorted
colnames(matrix_500) <- unique_land_covers_sorted

# for runs above 800 threshold
# colnames: land cover
# rownames: habitat
# cell values: number of runs where positive_runs > 800
matrix_800 <- matrix(0, nrow = length(unique_habitats_sorted), ncol = length(unique_land_covers_sorted))
rownames(matrix_800) <- unique_habitats_sorted
colnames(matrix_800) <- unique_land_covers_sorted

# Fill the matrices with threshold counts
for(hab in unique_habitats_sorted) {
  for(lc in unique_land_covers_sorted) {
    # Get all iterations for this habitat-land cover pair
    pair_data <- output1_data[output1_data$habitat == hab & output1_data$land_cover == lc, ]
    
    if(nrow(pair_data) > 0) {
      # Count runs above each threshold
      matrix_500[hab, lc] <- sum(pair_data$positive_runs > 500)
      matrix_800[hab, lc] <- sum(pair_data$positive_runs > 800)
    }
  }
}

# Save the threshold matrices
#write.csv(matrix_500, paste0(dft_folder, "/matrix_runs_above_500.csv"))
#write.csv(matrix_800, paste0(dft_folder, "/matrix_runs_above_800.csv"))


# heatmap for runs above 500
matrix_500_long <- melt(matrix_500)
colnames(matrix_500_long) <- c("Habitat", "Land_Cover", "Runs_Above_500")

# Get unique habitats and land covers in custom/alphabetical order
unique_habitats_ordered <- unique_habitats_sorted
unique_land_covers_ordered <- sort(unique(matrix_500_long$Land_Cover))

# Set factor levels to ensure custom ordering for habitats and alphabetical for land covers
matrix_500_long$Habitat <- factor(matrix_500_long$Habitat, levels = unique_habitats_ordered)
matrix_500_long$Land_Cover <- factor(matrix_500_long$Land_Cover, levels = unique_land_covers_ordered)

heatmap_500 <- ggplot(matrix_500_long, aes(y = Land_Cover, x = Habitat, fill = Runs_Above_500)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "orange", name = "Runs > 500") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Number of Runs Above 500 Threshold",
       subtitle = paste("Based on", length(all_results), "iterations"),
       x = "Habitat Type",
       y = "Land Cover Type")
heatmap_500

# heatmap for runs above 800
matrix_800_long <- melt(matrix_800)
colnames(matrix_800_long) <- c("Habitat", "Land_Cover", "Runs_Above_800")

# Set factor levels to ensure alphabetical ordering (using the same ordered levels as matrix_500)
matrix_800_long$Habitat <- factor(matrix_800_long$Habitat, levels = unique_habitats_ordered)
matrix_800_long$Land_Cover <- factor(matrix_800_long$Land_Cover, levels = unique_land_covers_ordered)

heatmap_800 <- ggplot(matrix_800_long, aes(x = Habitat, y = Land_Cover, fill = Runs_Above_800)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red", name = "Runs > 800") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Number of Runs Above 800 Threshold",
       subtitle = paste("Based on", length(all_results), "iterations"),
       x = "Habitat Type",
       y = "Land Cover Type")

heatmap_800

# save heatmaps
#ggsave(paste0(dft_folder, "/viz/heatmap_runs_above_500.png"), heatmap_500, width = 12, height = 8, dpi = 300)
#ggsave(paste0(dft_folder, "/viz/heatmap_runs_above_800.png"), heatmap_800, width = 12, height = 8, dpi = 300)


## ====== output 1.5 =======
# get number of positive iterations for each habitat - land cover pair out of 1000 for 100 iterations
df_hab_lc_counts <- data.frame()
for(i in seq_along(all_results)){
  pos <- all_results[[i]]$pos
  for(lc in unique_land_covers_ordered){
    pos_lc <- pos %>% filter(land_cover==lc)
    for(hab in unique_habitats_ordered){
      df_hab_lc_counts <- rbind(df_hab_lc_counts, data.frame(
        pair = paste0(hab, '_', lc),
        iteration = i,
        counts = pos_lc[hab][[1]]
      ))
    }
  }
}

# ggplot create histograms of df_hab_lc_counts, facet_wrap by pair, alphabetical order
ggplot(df_hab_lc_counts, aes(x = counts)) +
  geom_histogram(alpha = 0.7, position = "identity", bins = 30) +
  facet_wrap(~pair) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank()
  )


library(DescTools)
library(e1071)
# get mean, math mode, sd, skewness, and variation for each pair
# also get 2.5% and 97.5% quantiles
summary_stats_hab_lc <- df_hab_lc_counts %>%
  group_by(pair) %>%
  summarise(mean = mean(counts),
            # mode = Mode(counts),
            sd = sd(counts),
            skewness = skewness(counts),
            variation = var(counts),
            q2.5 = quantile(counts, 0.025),
            q97.5 = quantile(counts, 0.975))
hist(summary_stats_hab_lc$mean)
hist(summary_stats_hab_lc$sd)

group1 <- summary_stats_hab_lc %>% filter(q97.5 < 500)
group2 <- summary_stats_hab_lc %>% filter(q2.5 > 500)
group3 <- summary_stats_hab_lc %>% filter(q2.5 <= 500 & q97.5 >=500)

# order summary_stats_hab_lc by high mean, low variance, and left skewed
summary_stats_hab_lc <- summary_stats_hab_lc %>%
  arrange(desc(mean), sd, skewness)


## ======= test cluster-based separation =======
# Cluster-based separation to identify high mean, low variance, left skewed distributions
# prepare data for clustering (scale the features)
clustering_data <- summary_stats_hab_lc %>%
  mutate(skewness = ifelse(is.na(skewness), 0, skewness),
         CI_width = q97.5 - q2.5) %>%
  select(mean, skewness, CI_width) %>%
  scale()
### ====== test for normality ======
# shapiro-Wilk test for normality
shapiro_test <- shapiro.test(clustering_data[,2])
print(paste("Shapiro-Wilk test p-value:", round(shapiro_test$p.value, 4)))
print(paste("Shapiro-Wilk test statistic:", round(shapiro_test$statistic, 4)))

# Kolmogorov-Smirnov test for normality
ks_test <- ks.test(clustering_data[,2], "pnorm", 
                   mean = mean(clustering_data[,2]), 
                   sd = sd(clustering_data[,2]))
print(paste("Kolmogorov-Smirnov test p-value:", round(ks_test$p.value, 4)))

# not normal at all

#### ====== determine optimal k using elbow plot/ silhouette analysis
# determine optimal number of clusters using elbow method
wss <- numeric(10)
for(i in 1:10) {
  tryCatch({
    kmeans_result <- kmeans(clustering_data, centers = i, nstart = nrow(clustering_data))
    wss[i] <- kmeans_result$tot.withinss
  }, error = function(e) {
    print(paste("Error with k =", i, ":", e$message))
    wss[i] <- NA
  })
}

# plot elbow curve to determine optimal k
elbow_plot <- ggplot(data.frame(k = 1:10, wss = wss), aes(x = k, y = wss)) +
  geom_line() + geom_point() +
  theme_minimal() +
  labs(title = "Elbow Method for Optimal k",
       x = "Number of Clusters (k)",
       y = "Total Within Sum of Squares")
elbow_plot
# a bit ambiguous

# use silhouette analysis to determine optimal k
# calculate silhouette scores for different k values
sil_scores <- numeric(9)
sil_count_negative <- numeric(9)
for(k in 2:10) {
  km <- kmeans(clustering_data, centers=k, nstart=25)
  ss <- silhouette(km$cluster, dist(clustering_data))
  sil_scores[k-1] <- mean(ss[,3])
  sil_count_negative[k-1] <- length(ss[ss[,3]<0,])
}

sil_results <- data.frame(k=2:10, score=sil_scores, count=sil_count_negative)
# plot silhouette scores & counts side by side two panels
p1 <- ggplot(sil_results, aes(x=k, y=score)) +
  geom_line() + geom_point() +
  theme_minimal() +
  labs(title="Silhouette Scores",
       x="Number of Clusters (k)", 
       y="Average Silhouette Score")

p2 <- ggplot(sil_results, aes(x=k, y=count)) +
  geom_line(color="red") + geom_point(color="red") +
  theme_minimal() +
  labs(title="Negative Silhouette Counts",
       x="Number of Clusters (k)",
       y="Count of Negative Silhouettes")

gridExtra::grid.arrange(p1, p2, ncol=2)


# find all k values with minimum negative silhouette count
min_count <- min(sil_results$count)
vals1 <- sil_results$k[sil_results$count == min_count]

max_scores <- max(sil_results$score)
vals2 <- sil_results$k[sil_results$score == max_scores]
# find intersection of vals1 and vals_2
optimal_k <- intersect(vals1, vals2)

##====== perform k-means clustering with optimal k
kclust <- kmeans(clustering_data, centers = 4, nstart = nrow(clustering_data))

# Add cluster assignments to the original data
# create a mapping for the clean data
clean_indices <- which(complete.cases(clustering_data) & 
                      !is.infinite(rowSums(as.matrix(clustering_data))))

# create cluster assignments for all data (NA for excluded rows)
all_clusters <- rep(NA, nrow(summary_stats_hab_lc))
all_clusters[clean_indices] <- kclust$cluster

summary_stats_hab_lc_clustered <- summary_stats_hab_lc %>%
  mutate(
    skewness = ifelse(is.na(skewness), 0, skewness),
    cluster = as.factor(all_clusters),
    left_skew = as.factor(ifelse(skewness<0, 1, ifelse(skewness==0, 0, -1))))


write.csv(summary_stats_hab_lc_clustered, 'aoh_results/matrix_4clusters_0730.csv',
          row.names = FALSE)

cluster_summary <- summary_stats_hab_lc_clustered %>%
  filter(!is.na(cluster)) %>%  # remove rows that couldn't be clustered
  group_by(cluster) %>%
  summarise(
    n_pairs = n(),
    avg_mean = mean(mean),
    avg_sd = mean(sd),
    avg_variation = mean(variation),
    avg_skewness = mean(skewness),
    avg_q2.5 = mean(q2.5),
    avg_q97.5 = mean(q97.5),
    cv = avg_sd / avg_mean,  # coefficient of variation (lower = more consistent)
    left_skew_score = avg_skewness<0  # TRUE if left skewed
  ) %>%
  arrange(desc(avg_mean), cv, avg_skewness) # sort by high mean, low CV, then left skew

print(cluster_summary)
cluster_order <- cluster_summary$cluster


# Identify the cluster with highest mean, lowest variance, and left skew
high_mean_low_var_left_skew_cluster <- cluster_summary %>%
  filter(avg_mean == max(avg_mean)) %>%
  filter(cv == min(cv)) %>%
  filter(avg_skewness == min(avg_skewness)) %>%
  pull(cluster)

# Create groups based on clustering
group_high_mean_low_var_left_skew <- summary_stats_hab_lc_clustered %>%
  filter(cluster == high_mean_low_var_left_skew_cluster)

group_high_mean_low_var_left_skew$pair
# create custom group
a_group <- summary_stats_hab_lc_clustered %>% filter(cluster == 4)
a_group$pair

# Plot the clustered data
cluster_plot <- ggplot(summary_stats_hab_lc_clustered %>% filter(!is.na(cluster)), 
                       aes(x = mean, y = sd, color = cluster, 
                           size = q97.5 - q2.5, shape = left_skew)) +
  geom_point(alpha = 0.7) +
  scale_color_discrete(name = "Cluster") +
  scale_size_continuous(name = "95% CI Width") +
  scale_shape_discrete(name = "Skewness", labels = c("Right", "None", "Left")) +
  theme_minimal() +
  labs(title = "K-means Clustering of Habitat-Land Cover Pairs",
       x = "Mean Count",
       y = "Standard Deviation")

cluster_plot

# plot a group
ggplot(a_group %>% 
         arrange(desc(mean), sd, skewness), 
       aes(x = mean, y = reorder(pair, mean, decreasing = FALSE))) +
  geom_pointrange(
    size = 0.3,
    aes(xmin = q2.5, xmax = q97.5, color = sd),
    position = position_dodge(width = 2)) +
  scale_color_gradient(low = "green", high = "red", name = "SD") +
  scale_x_continuous(limits = c(0, 1000)) +
  labs(title = "Distributions",
       subtitle = paste("Cluster", unique(a_group$cluster)[[1]], "-", nrow(a_group), "pairs"),
       y = "Habitat-Land Cover Pair",
       x = "Mean Count")

## ====== silouhette analysis ======
# Silhouette analysis for cluster validation
library(cluster)

# Calculate silhouette for the current clustering
silhouette_obj <- silhouette(as.numeric(as.character(kclust$cluster)), 
                            dist(clustering_data))

# Calculate average silhouette width
avg_silhouette <- mean(silhouette_obj[,3])
print(paste("Average silhouette width:", round(avg_silhouette, 3)))

# Create silhouette plot using base R and ggplot2
silhouette_df <- data.frame(
  cluster = silhouette_obj[,1],
  mean = clustering_data[, 1],
  sd = clustering_data[, 2],
  neighbor = silhouette_obj[,2],
  sil_width = silhouette_obj[,3],
  point_id = 1:nrow(silhouette_obj)
)

# Order by cluster and silhouette width for better visualization
silhouette_df <- silhouette_df %>%
  mutate(cluster = factor(cluster, levels = cluster_order)) %>%
  arrange(cluster, desc(mean), sd) %>%
  mutate(ordered_id = 1:n())

# create custom silhouette plot
# order rows by cluster order in cluster summary
silhouette_plot <- ggplot(silhouette_df, 
                          aes(x = sil_width, y = ordered_id)) +
  geom_segment(aes(x = 0, xend = sil_width, 
                   y = nrow(silhouette_df)-ordered_id, 
                   yend = nrow(silhouette_df)-ordered_id, 
                   color = cluster), 
               size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  scale_color_discrete(name = "Cluster") +
  theme_minimal() +
  labs(title = "Silhouette Analysis",
       subtitle = paste("Average silhouette width:", round(avg_silhouette, 3)),
       x = "Silhouette width",
       y = "Points") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
silhouette_plot
## ======= plot interactive heatmap of summary_stats_hab_lc =======
# x is habitat, y is land cover, color is cluster group 
# show mean, sd, skewness when mouse hovers over the point
# show as an interactive visualization
library(plotly)

# Create interactive heatmap of summary_stats_hab_lc with cluster colors
# Parse pair names to extract habitat and land cover
summary_stats_hab_lc_clustered_parsed <- summary_stats_hab_lc_clustered %>%
  filter(!is.na(cluster)) %>%
  mutate(
    # Split pair name by underscore and take first part as habitat, rest as land cover
    habitat = sapply(strsplit(pair, "_"), function(x) x[1]),
    land_cover = sapply(strsplit(pair, "_"), function(x) paste(x[-1], collapse = "_"))
  )

# Create interactive heatmap with habitat on x-axis and land cover on y-axis
# order land cover by alphabetical order ascending
# color of cluster goes from red to blue, by the ordering in cluster_summary

n_clusters <- length(cluster_order)
cluster_colors <- colorRampPalette(c("red", "orange", "white", "lightblue", "blue"))(n_clusters)
names(cluster_colors) <- cluster_order

# Order land cover alphabetically and create factor levels
summary_stats_hab_lc_clustered_parsed <- summary_stats_hab_lc_clustered_parsed %>%
  mutate(
    land_cover = factor(land_cover, levels = sort(unique(land_cover))),
    cluster = factor(cluster, levels = cluster_order)
  )



summary_stats_hab_lc_clustered_parsed <- summary_stats_hab_lc_clustered_parsed %>%
  mutate(CI_width = round(q97.5 - q2.5, 0))

interactive_heatmap <- ggplot(summary_stats_hab_lc_clustered_parsed, 
                              aes(x = habitat, 
                                  y = land_cover, 
                                  fill = cluster,
                                  text = paste(
                                    #"Habitat:", habitat,
                                     #         "\nLand Cover:", land_cover,
                                              "Cluster:", cluster,
                                              "\nMean:", round(mean, 2),
                                              "\nSD:", round(sd, 2),
                                              "\nSkewness:", round(skewness, 3),
                                              "\nCI width:", CI_width
                                              ,"\nQ2.5:", round(q2.5, 2)
                                              ,"\nQ97.5:", round(q97.5, 2)
                                  )))+
  geom_tile() +
  scale_fill_manual(values = cluster_colors, name = "Cluster") +
  geom_text(data = summary_stats_hab_lc_clustered_parsed,
            aes(label = CI_width), color = "white",
            size = 3) +
  labs(title = "Clusters",
       x = "Habitat",
       y = "Land Cover") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        panel.grid = element_blank())
interactive_heatmap
# Convert to interactive plot
interactive_heatmap_plotly <- ggplotly(interactive_heatmap, tooltip = "text")
interactive_heatmap_plotly



##======== output 2: AUC values matrix =======
# colnames: land cover
# rownames: iteration 1, 2, ..., 100
# cell values: auc values

# get unique land covers from the first iteration
first_iter_pos <- all_results[[1]]$pos
land_covers <- unique(first_iter_pos$land_cover)

# create matrix for AUC values
output2_matrix <- matrix(nrow = length(all_results), ncol = length(land_covers))
colnames(output2_matrix) <- land_covers
rownames(output2_matrix) <- paste0("iteration_", 1:length(all_results))

# fill the matrix with AUC values
for(iter in seq_along(all_results)) {
  pos_data <- all_results[[iter]]$pos
  
  for(lc in land_covers) {
    lc_data <- pos_data[pos_data$land_cover == lc, ]
    if(nrow(lc_data) > 0) {
      output2_matrix[iter, lc] <- lc_data$auc[1]
    } else {
      output2_matrix[iter, lc] <- NA
    }
  }
}

# summary statistics for each column in output 2
summary_stats_output2 <- data.frame(
  land_cover = colnames(output2_matrix),
  mean_auc = colMeans(output2_matrix, na.rm = TRUE),
  sd_auc = apply(output2_matrix, 2, sd, na.rm = TRUE),
  min_auc = apply(output2_matrix, 2, min, na.rm = TRUE),
  max_auc = apply(output2_matrix, 2, max, na.rm = TRUE)
)

# format summary_stats_output2 to have 2 digits
summary_stats_output2 <- summary_stats_output2 %>%
  mutate(across(where(is.numeric), ~ round(., 2)))

# save output 2
write.csv(output2_matrix, paste0(dft_folder, "/output2_auc_matrix.csv"))

## ===== output 3 tertile sensitivity ======
# output 3: tertile sensitivity
tmp <- all_results[[1]]$tert
summary(tmp$two.third)
summary(tmp$one.third)
hist(tmp$one.third)
hist(tmp$two.third)

# for each iteration, record min mean max
one_third_data <- data.frame()
two_third_data <- data.frame()

for(iter in seq_along(all_results)) {
  tert_data <- all_results[[iter]]$tert
  
  one_third_data <- rbind(one_third_data, 
                          data.frame(min_val = min(tert_data$one.third),
                                     mean_val = mean(tert_data$one.third),
                                     max_val = max(tert_data$one.third)))
  
  two_third_data <- rbind(two_third_data, 
                          data.frame(min_val = min(tert_data$two.third),
                                     mean_val = mean(tert_data$two.third),
                                     max_val = max(tert_data$two.third)))
}

one_third_data$iteration <- seq(1, nrow(one_third_data), 1)
two_third_data$iteration <- seq(1, nrow(two_third_data), 1)

# convert to long format for plotting
one_third_long <- melt(one_third_data, id.vars = "iteration", 
                       variable.name = "statistic", value.name = "value")

# create histogram for one.third with different colors and varying x-axis ranges
one_third_hist <- ggplot(one_third_long, aes(x = value, fill = statistic)) +
  geom_histogram(alpha = 0.7, position = "identity", bins = 30) +
  scale_fill_manual(values = c("min_val" = "red", "mean_val" = "blue", "max_val" = "green"),
                    labels = c("min_val" = "Minimum", "mean_val" = "Mean", "max_val" = "Maximum"),
                    name = "Statistic") +
  theme_minimal() +
  labs(title = "Distribution of One-Third Tertile Values",
       subtitle = paste("Based on", length(all_results), "iterations"),
       x = "Value",
       y = "Frequency") +
  facet_wrap(~statistic, scales = "free_x", labeller = labeller(
    statistic = c("min_val" = "Minimum", "mean_val" = "Mean", "max_val" = "Maximum")
  )) +
  # add individual x-axis scales for each facet
  scale_x_continuous(breaks = function(x) pretty(x, n = 5))

one_third_hist

# convert to long format for plotting
two_third_long <- melt(two_third_data, id.vars = "iteration", 
                       variable.name = "statistic", value.name = "value")

# create histogram for two.third with different colors and varying x-axis ranges
two_third_hist <- ggplot(two_third_long, aes(x = value, fill = statistic)) +
  geom_histogram(alpha = 0.7, position = "identity", bins = 30) +
  scale_fill_manual(values = c("min_val" = "red", "mean_val" = "blue", "max_val" = "green"),
                    labels = c("min_val" = "Minimum", "mean_val" = "Mean", "max_val" = "Maximum"),
                    name = "Statistic") +
  theme_minimal() +
  labs(title = "Distribution of Two-Third Tertile Values",
       subtitle = paste("Based on", length(all_results), "iterations"),
       x = "Value",
       y = "Frequency") +
  facet_wrap(~statistic, scales = "free_x", labeller = labeller(
    statistic = c("min_val" = "Minimum", "mean_val" = "Mean", "max_val" = "Maximum")
  )) +
  # add individual x-axis scales for each facet
  scale_x_continuous(breaks = function(x) pretty(x, n = 5))

two_third_hist