# Compute null models for regressions
# Author: Wenxin Yang
# Date: May, 2025

# ================== Prep =====================
# Load libraries
packages <- c("sf", "terra", "here", "rgrass", "dplyr", "data.table", 
              "rredlist", "tidyr", "Matrix", "tidyverse", "pROC",
              "kableExtra", "patchwork")
lapply(packages, library, character.only = TRUE)

setwd('/')
setwd('Users/wenxinyang/Desktop/GitHub/colander')
source('scripts/refineBiomodelos/funcs.R')
source('scripts/refineBiomodelos/info.R')
df_seed = 2025

# Create viz directory if it doesn't exist
if (!dir.exists("viz")) {
  dir.create("viz")
}

# ================== Read in data =====================
d_near <- 0.0001
df <- read_csv('data/occ_pts/allinfo_ideam_cgls.csv') %>% 
  filter(dst_t_b > d_near) %>% select(-dst_t_b)
df <- df %>% select(-all_of(c('hab_2', 'hab_8')))

df <- df %>% mutate(
  hab_14.12 = ifelse(hab_14.1+hab_14.2 >0, 1, 0),
  hab_14.36 = ifelse(hab_14.3+hab_14.6 >0, 1, 0),
  hab_14.45 = ifelse(hab_14.4+hab_14.5 >0, 1, 0)
)

# remove the original ones
df <- df %>% select(-all_of(c('hab_14.1', 'hab_14.2', 'hab_14.3', 'hab_14.4', 'hab_14.5', 'hab_14.6', 'hab_')))
df <- df %>% select(-all_of(drop_cols))

lc_codes <- unique(df$nvl_2_n)
lc_codes

# for each land cover type, randomize 1 and 0 values for 1000 times and fit regressions

getRegInfo <- function(mod, lc_code){
  # Get model coefficients summary and AUC
  model_summary <- summary(mod)$coefficients
  
  # Get variable names from the model
  var_names <- names(coef(mod))[-1]  # exclude intercept
  
  # Initialize vectors for odds ratios and p-values
  odds_ratios <- rep(0, length(var_names))
  p_values <- rep(1, length(var_names))  # default to 1 (not significant)
  names(odds_ratios) <- var_names
  names(p_values) <- var_names
  
  # Fill in values for non-NA coefficients
  valid_coefs <- which(!is.na(coef(mod)[-1]))
  if(length(valid_coefs) > 0) {
    odds_ratios[valid_coefs] <- exp(coef(mod)[-1][valid_coefs])
    p_values[valid_coefs] <- model_summary[-1, "Pr(>|z|)"][valid_coefs]
  }
  
  # Round values
  odds_ratios <- round(odds_ratios, 2)
  p_values <- round(p_values, 3)
  
  # Create table with odds ratios and p-values
  odds_ratios_table <- data.frame(
    odds_ratios = odds_ratios,
    p_values = p_values,
    row.names = var_names
  )
  
  # Set column names and get land cover names
  colnames(odds_ratios_table) <- c(paste0("lc_", lc_code), "p_values")
  lc_names <- rownames(odds_ratios_table)
  
  return(odds_ratios_table)
}

get_null_info <- function(df, lc_code, seed = 226, lc_dataset){
  print(lc_code)
  splits_lc <- create_balanced_splits(lc_dataset, df, lc_code, seed)
  
  # get the train set
  train0 <- splits_lc$train
  
  # fit a model for the observed pattern
  obs_mod <- glm(presence ~ ., data = train0, family = 'binomial')
  obs_info <- getRegInfo(obs_mod, lc_code)
  obs_info$seed <- 0
  obs_info$hab_code <- rownames(obs_info)
  rownames(obs_info) <- seq(1:nrow(obs_info))
  
  # fit models for random patterns
  for(i in 1:1000){
    this_seed <- i
    # randomly assign 1 and 0 to the presence while keeping the sum the same
    this_train <- train0
    this_train$presence <- sample(this_train$presence)
    # fit a model for the random pattern
    this_mod <- glm(presence ~ ., data = this_train, family = 'binomial')
    this_info <- getRegInfo(this_mod, lc_code)
    this_info$seed <- i
    this_info$hab_code <- rownames(this_info)
    rownames(this_info) <- seq(1:nrow(this_info))
    
    obs_info <- rbind(obs_info, this_info)
  }
 
   return(obs_info)
}

plot_all_habitats_facet <- function(data, lc_code, lc_name) {
  # Calculate quantiles and actual values for all habitats
  data %>%
    group_by(hab_code) %>%
    ggplot(aes(x = .data[[lc_code]])) +
    geom_histogram(fill = "#B0C4B1", color = "#4A5759") +
    geom_vline(data = . %>% filter(seed == 0),
               aes(xintercept = .data[[lc_code]], color = "Observed Value"),
               linetype = "dashed") +
    # Add 95% CI lines
    geom_vline(data = . %>% 
                 filter(seed > 0) %>%
                 summarise(q025 = quantile(.data[[lc_code]], 0.025),
                           q975 = quantile(.data[[lc_code]], 0.975)),
               aes(xintercept = q025, color = "95% CI"),
               linetype = "dotted") +
    geom_vline(data = . %>%
                 filter(seed > 0) %>% 
                 summarise(q025 = quantile(.data[[lc_code]], 0.025),
                           q975 = quantile(.data[[lc_code]], 0.975)),
               aes(xintercept = q975, color = "95% CI"),
               linetype = "dotted") +
    scale_color_manual(name = "Lines",
                      values = c("Observed Value" = "red", 
                               "95% CI" = "blue")) +
    facet_wrap(~hab_code, scales = "free") +
    labs(x = "Coefficient Value",
         y = "Count", 
         title = paste("Null Distribution for", lc_name)) +
    theme_minimal(base_family = "Times New Roman")
}


# for grass/shrub land cover
info_32 <- get_null_info(df, "32", df_seed, "ideam")
unique(info_32$hab_code)
# plot all
plot_all_habitats_facet(info_32, "lc_32", "Grass/Shrub")
# plot only those with significant coefficients
info_32_sig_rand <- info_32 %>% filter(seed > 0 & p_values < 0.05)
info_32_obs <- info_32 %>% filter(seed==0)
info_32_sig <- rbind(info_32_sig_rand, info_32_obs)
plot_all_habitats_facet(info_32_sig, "lc_32", "Grass/Shrub")