# Parameterized version of 02_computeMatrix.R for single matrix generation
# with random sampling and statistics
# Author: Wenxin Yang
# Date: July, 2025
# THIS IS THE SCRIPT I ENDED UP USING FOR THE MANUSCRIPT

# ================== Function to run single matrix analysis with parameters =====================
run_analysis <- function(
    # core parameters
  num_generalist = 7,
  d_near = 0,
  random_seed = 2025,
  # data filtering parameters
  balance_specialist_generalist = 1, # balance specialist vs non specialist by taxa
  remove_desert_rocky_aa = TRUE, # remove desert (hab_6), rocky (hab_8), and artificial-aquatic (hab_15) habitats
  # output parameters
  output_suffix = paste0("_", format(Sys.Date(), "%m%d"), "_", format(Sys.Date(), "%Y")),
  save_results = TRUE,
  dft_folder = "aoh_results_randomCI",
  
  # parallel processing parameters
  n_cores = NULL
) {
  
  # ================== Prep =====================
  # load libraries
  packages <- c("sf", "terra", "here", "rgrass", "dplyr", "data.table", 
                "rredlist", "tidyr", "Matrix", "tidyverse", "pROC",
                "kableExtra", "ggplot2", "reshape2", "parallel")
  lapply(packages, library, character.only = TRUE)
  
  setwd('/')
  setwd('Users/wenxinyang/Desktop/GitHub/colander')
  source('scripts/refineBiomodelos/funcs.R')
  source('scripts/refineBiomodelos/info.R')
  
  # create aoh_results directory structure if it doesn't exist
  if (!dir.exists(dft_folder)) dir.create(dft_folder)
  if (!dir.exists(paste0(dft_folder, "/viz"))) dir.create(paste0(dft_folder, "/viz"))
  
  
  # ======= 2. Read in & prep data =======
  df_all <- read.csv('data/occ_pts/allinfo_ideam_cgls_coords_2022.csv')
  # df_all <- read.csv('data/occ_pts/allinfo_ideam_cgls_coords.csv') # 2018
  
  if(!'nvl_2_n' %in% colnames(df_all)){
    df_all$nvl_2_n <- df_all$nivel_2_num
  }
  
  if(!'31' %in% unique(df_all$nvl_2_n)){
    df_all$nvl_2_n <- unlist(lapply(df_all$nvl_2_n, function(x) gsub("\\.", "", x)))
  }
  
  # add species generalist vs specialist info
  df_pref <- read.csv('data/occ_pts/col_animal_pref_cleaned.csv')
  df_basic_info <- addGeneralistInfo(df_pref, num_generalist)
  rm(df_pref)
  
  # remove habitat codes not included in the analysis per the original paper
  df_all_info <- df_all %>% select(-all_of(drop_cols))
  
  # remove those near boundary
  if('dst_t_b' %in% colnames(df_all_info)){
    summary(df_all_info$dst_t_b)
    df_all_info <- df_all_info %>% filter(dst_t_b > d_near) %>% select(-dst_t_b)
  }

  
  # merge several artificial habitat types
  df_all_info <- df_all_info %>% mutate(
    hab_14.12 = ifelse(hab_14.1+hab_14.2 >0, 1, 0),
    hab_14.36 = ifelse(hab_14.3+hab_14.6 >0, 1, 0),
    hab_14.45 = ifelse(hab_14.4+hab_14.5 >0, 1, 0)
  )
  
  # remove the original ones
  df_all_info <- df_all_info %>% select(-all_of(c('hab_14.1', 'hab_14.2', 'hab_14.3', 'hab_14.4', 'hab_14.5', 'hab_14.6')))
  
  # ======= 2.3 Remove some habitats if requested =======
  if (remove_desert_rocky_aa) {
    df_all_info <- df_all_info %>% select(-all_of(c('hab_6', 'hab_8', 'hab_15')))
    # Update colname_dataset to exclude these habitats
    colname_dataset <- c("lc_code", "n_samples", habitat_info1$habitat_code[!habitat_info1$habitat_code %in% c('hab_6', 'hab_8', 'hab_15')], "auc")
  } else {
    colname_dataset <- c("lc_code", "n_samples", habitat_info1$habitat_code, "auc")
  }
  
  # remove rows with NA values
  df_all_info <- df_all_info[complete.cases(df_all_info), ]
  sum(is.na(df_all_info))
  df_all_info$X <- NULL
  
  # sum all columns that start with "hab_"
  df_all_info$sum <- rowSums(df_all_info[,grep("^hab_", colnames(df_all_info))])
  df_all_info <- df_all_info %>% filter(sum>0) %>% select(-sum)
  
  # colSums(df_all_info[,grep("^hab_", colnames(df_all_info))])
  
  ## ======= 2.1 remove generalist species =======
  li_generalist <- unique((df_basic_info %>% filter(type == 'generalist'))$name)
  df_all_info <- df_all_info %>% filter(!species %in% li_generalist)
  
  ## ========= 4.1 Generate 1000 matrices in parallel =======
  cat("Setting up parallel processing for 1000 iterations...\n")
  
  # Set up parallel processing
  if (is.null(n_cores)) {
    n_cores <- detectCores() - 3  # Leave three cores free
  }
  n_cores <- min(n_cores, detectCores() - 1, 1000)  # Don't use more cores than available or needed
  
  cat("Using", n_cores, "cores for parallel processing\n")
  # Create cluster
  cl <- makeCluster(n_cores)
  
  # Export necessary data and functions to cluster
  clusterExport(cl, c("habitat_info1", "ideam_lc_info", "cgls_lc_info", 
                      "get_a_single_row", "create_balanced_splits", 
                      "build_evaluate_model", "get_odds_ratios_row"))
  
  # Load required packages and define functions on cluster
  clusterEvalQ(cl, {
    library(dplyr)
    library(tidyr)
    library(pROC)
    library(stringr)
    
    # Define generate_matrix function on cluster
  generate_matrix <- function(dat, landcover_colname, landcover_dataset, seed, cols){
    # dat, 'nvl_2_n', 'ideam', seed, colname_dataset
    tryCatch({
        # Debug: check the data structure
        if (!landcover_colname %in% colnames(dat)) {
          stop("Land cover column '", landcover_colname, "' not found in data")
        }
        
        # Get unique land cover values and ensure they are character
        landcover_values <- dat[[landcover_colname]]
        if (!is.character(landcover_values)) {
          landcover_values <- as.character(landcover_values)
        }
        unique_landcovers <- unique(landcover_values)
        
        # Remove any NA or empty values
        unique_landcovers <- unique_landcovers[!is.na(unique_landcovers) & unique_landcovers != ""]
        
        if (length(unique_landcovers) == 0) {
          stop("No valid land cover values found")
        }
        
        info_list <- lapply(unique_landcovers,
                          function(x) {
                            tryCatch({
                              get_a_single_row(dat, x, seed = seed, lc_dataset = landcover_dataset)
                            }, error = function(e) {
                              cat("Error processing land cover", x, ":", e$message, "\n")
                              # Return empty row with correct structure
                              empty_row <- data.frame(matrix(0, nrow = 1, ncol = length(cols)))
                              colnames(empty_row) <- cols
                              empty_row$lc_code <- paste0("lc_", x)
                              empty_row$auc <- 0.5
                              empty_row$n_samples <- 0
                              return(empty_row)
                            })
                          })
      
      # Filter out NULL results and combine
      info_list <- info_list[!sapply(info_list, is.null)]
      
      if (length(info_list) == 0) {
        # If no valid results, return empty matrix
        empty_matrix <- data.frame(matrix(0, nrow = 0, ncol = length(cols)))
        colnames(empty_matrix) <- cols
        return(empty_matrix)
      }
      
      info_list <- do.call(rbind, info_list) %>% select(all_of(cols))
      # reset row index
      rownames(info_list) <- NULL
      
      # change column names that start with "hab_" to the habitat names in habitat_info
      colnames(info_list) <- sapply(colnames(info_list), function(x){
        if (x %in% habitat_info1$habitat_code) {
          return(habitat_info1$habitat_name[match(x, habitat_info1$habitat_code)])
        } else {
          return(x)
        }
      })
      
      # change lc_code values to the land cover names in ideam_lc_info
      info_list$lc_code <- sapply(info_list$lc_code, function(x) {
        if(landcover_dataset=='ideam'){
          return(ideam_lc_info$ideam_lc_name[match(x, ideam_lc_info$ideam_lc_code)])
        } else if(landcover_dataset=='cgls'){
          return(cgls_lc_info$cgls_lc_name[match(x, cgls_lc_info$cgls_lc_code)])
        }
      })
      
      # rename lc_code to "land cover"
      colnames(info_list)[colnames(info_list)=='lc_code'] <- 'Land Cover'
      
      return(info_list)
    }, error = function(e) {
      cat("Error in generate_matrix:", e$message, "\n")
      # Return empty matrix with correct structure
      empty_matrix <- data.frame(matrix(0, nrow = 0, ncol = length(cols)))
      colnames(empty_matrix) <- cols
      return(empty_matrix)
    })
  }
  
    # Define getHabAssoc function on cluster
  getHabAssoc <- function(df, seed, landcover_dataset){
    ## determine what is a good land cover - habitat pair to keep
    ## using the Lumbierres et al. thresholds
    habitat_associations <- data.frame(
      land_cover = df$`Land Cover`,
      pos_odds_habitats = sapply(1:nrow(df), function(i) {
        # Get habitat columns (exclude irrelevant columns)
        habitat_cols <- colnames(df)[!colnames(df) %in% c("Land Cover", "n_samples", "auc")]
        # Get odds ratios for this land cover type
        odds_ratios <- as.numeric(df[i, habitat_cols])
        # Get habitats with positive odds ratios
        pos_odds <- habitat_cols[odds_ratios > 1]
        paste(pos_odds, collapse=", ")
      })
    )
    habitat_associations$seed_val <- seed
    habitat_associations$dataset <- landcover_dataset
    # add auc
    habitat_associations$auc <- df$auc
    # add n_samples
    habitat_associations$n_samples <- df$n_samples
    
    return(habitat_associations)
  }
  
  })
  
  # Function to run a single iteration with a specific seed
  run_single_iteration <- function(seed, df_all_info, df_basic_info, colname_dataset, balance_specialist_generalist) {
    # seed, df_all_info, df_basic_info, colname_dataset, balance_specialist_generalist
    tryCatch({
      # Balance specialist vs non-specialist within each taxa for this iteration
      if (balance_specialist_generalist) {
        set.seed(seed)  # Use the iteration seed
        records <- merge(df_all_info, df_basic_info, by.x='species', by.y='name')
        balanced_rec <- data.frame(matrix(nrow=0, ncol=ncol(records)))
        colnames(balanced_rec) <- colnames(records)
        
        # Balance specialist vs non-specialist within each taxa
        for (t in unique(records$taxa)) {
          sub_records <- records %>% filter(taxa == t)
          set.seed(seed)  # Use the iteration seed for consistency
          numsamp <- balance_specialist_generalist*sum(sub_records$type=='specialist') # keep more info from non-specialist species if param >1
          sub_records_other <- sub_records %>% 
            filter(type == 'other') %>% 
            sample_n(sum(sub_records$type == 'specialist'))
          sub_records <- rbind(sub_records %>% filter(type != 'other'), sub_records_other)
          balanced_rec <- rbind(balanced_rec, sub_records)
        }
        
        # Remove duplicates that might occur from overlapping strata
        balanced_rec <- balanced_rec %>% distinct()
      } else {
        balanced_rec <- merge(df_all_info, df_basic_info, by.x='species', by.y='name')
      }
      
      # Generate matrix for this iteration
      dat <- balanced_rec
      
      # Debug: check data structure
      if (nrow(dat) == 0) {
        stop("No data available for iteration with seed ", seed)
      }
      
      if (!'nvl_2_n' %in% colnames(dat)) {
        stop("Land cover column 'nvl_2_n' not found in data for iteration with seed ", seed)
      }
      
      this_ideam <- generate_matrix(dat, 'nvl_2_n', 'ideam', seed, colname_dataset)
      
      all_values <- as.vector(unlist(this_ideam %>% select(-`Land Cover`, auc, n_samples)))
      all_values <- all_values[all_values > 1] # keep positive values
      
      if (length(all_values) == 0) {
        # If no positive values, create empty results
        this_tert <- data.frame(
          min = NA,`one-third` = NA,`two-third` = NA, max = NA,seed = seed
        )
        tmp <- data.frame(
          land_cover = character(),
          pos_odds_habitats = character(),
          seed_val = seed,
          dataset = 'ideam',
          auc = numeric(),
          n_samples = numeric()
        )
        
        # Create empty raw odds ratios dataframe
        raw_odds <- data.frame(
          land_cover = character(),
          habitat = character(),
          odds_ratio = numeric(),
          seed = numeric()
        )
      } else {
        this_tert <- data.frame(
          min = min(all_values),
          `one-third` = quantile(all_values, 1/3),
          `two-third` = quantile(all_values, 2/3), 
          max = max(all_values),
          seed = seed
        )
        
        tmp <- getHabAssoc(this_ideam, seed, 'ideam')
        
        # Extract raw odds ratios for this iteration
        raw_odds <- data.frame()
        for (i in 1:nrow(this_ideam)) {
          land_cover <- this_ideam$`Land Cover`[i]
          # Get habitat columns (exclude metadata columns)
          habitat_cols <- colnames(this_ideam)[!colnames(this_ideam) %in% c("Land Cover", "n_samples", "auc")]
          
          for (hab in habitat_cols) {
            odds_ratio <- this_ideam[i, hab]
            if (!is.na(odds_ratio)) {
              raw_odds <- rbind(raw_odds, data.frame(
                land_cover = land_cover,
                habitat = hab,
                odds_ratio = odds_ratio,
                seed = seed
              ))
            }
          }
        }
      }
      
      gc()
      
      return(list(tert = this_tert, hab_assoc = tmp, raw_odds = raw_odds))
    }, error = function(e) {
      # Return empty results on error
      cat("Error in iteration with seed", seed, ":", e$message, "\n")
      
      this_tert <- data.frame(
        min = NA,
        `one-third` = NA,
        `two-third` = NA, 
        max = NA,
        seed = seed
      )
      
      tmp <- data.frame(
        land_cover = character(),
        pos_odds_habitats = character(),
        seed_val = seed,
        dataset = 'ideam',
        auc = numeric(),
        n_samples = numeric()
      )
      
      # Create empty raw odds ratios dataframe
      raw_odds <- data.frame(
        land_cover = character(),
        habitat = character(),
        odds_ratio = numeric(),
        seed = numeric()
      )
      
      return(list(tert = this_tert, hab_assoc = tmp, raw_odds = raw_odds))
    })
  }
  
  # Define a simplified selectPairs function that only works with positive odds
  selectPairs <- function(df, posOrHigh){
    if(posOrHigh=='positive'){
      df_sel <- df %>%
        select(land_cover, pos_odds_habitats) %>%
        separate_rows(pos_odds_habitats, sep=", ") %>%
        group_by(land_cover, pos_odds_habitats) %>%
        summarise(count=n(), .groups='drop') %>%
        right_join(
          expand.grid(
            land_cover=unique(df$land_cover),
            pos_odds_habitats=unique(unlist(strsplit(df$pos_odds_habitats, ", ")))
          ),
          by=c("land_cover", "pos_odds_habitats")
        ) %>%
        mutate(count = replace_na(count, 0)) %>%
        pivot_wider(names_from=pos_odds_habitats, values_from=count, values_fill=0)
      
      # add missing habitats
      all_habitats <- unique(unlist(strsplit(df$pos_odds_habitats, ", ")))
      
      for (hab in all_habitats){
        if(!hab %in% colnames(df_sel)){
          df_sel[hab] = 0
        }
      }
      
      mean_auc_n_samples <- df %>%
        group_by(land_cover) %>%
        summarise(auc = mean(auc), n_samples = mean(n_samples), .groups='drop')
      
      df_sel <- merge(df_sel, mean_auc_n_samples, by='land_cover')
      
      return(df_sel)
    } else {
      stop("Only 'positive' is supported for posOrHigh parameter")
    }
  }
  
  # Run iterations in parallel (start with fewer for testing)
  n_iterations <- 1000  # Start with 10 for testing, can increase to 1000 later
  cat("Running", n_iterations, "iterations in parallel...\n")
  
  iteration_results <- parLapply(cl, 1:n_iterations, function(seed) {
    run_single_iteration(seed, df_all_info, df_basic_info, colname_dataset, balance_specialist_generalist)
  })
  
  # Stop cluster
  stopCluster(cl)
  
  # Combine all results
  btst_ideam <- do.call(rbind, lapply(iteration_results, function(x) x$hab_assoc))
  tert_ideam <- do.call(rbind, lapply(iteration_results, function(x) x$tert))
  raw_odds_all <- do.call(rbind, lapply(iteration_results, function(x) x$raw_odds))
  
  # ======= 6. Process and save results =======
  if (save_results) {
    # Create filename with parameters
    filename_base <- paste0(dft_folder, "/btst_ideam_randomCI_", 
                            "gen", num_generalist, "_",
                            ifelse(d_near > 0, paste0("near", d_near, "_"), ""),
                            ifelse(remove_desert_rocky_aa, "_noDesRockAA", ""),
                            ifelse(!balance_specialist_generalist, "_nobal", ""),
                            output_suffix)
    
    write.csv(btst_ideam, paste0(filename_base, ".csv"))
    write.csv(tert_ideam, paste0(filename_base, "_tertiles.csv"))
    
    # process results for visualization
    # count positive odds ratios pairs
    btst_ideam_pos <- selectPairs(btst_ideam, 'positive')
    write.csv(btst_ideam_pos, paste0(filename_base, "_pos.csv"))
    

    
    # Save all results as a single RDS file
    all_results <- list(
      iteration_results = iteration_results,  # All 1000 individual results
      combined_hab_assoc = btst_ideam,       # Combined habitat associations
      combined_tertiles = tert_ideam,        # Combined tertiles
      raw_odds_ratios = raw_odds_all,        # Raw odds ratios for CI calculations
      parameters = list(
        num_generalist = num_generalist,
        d_near = d_near,
        balance_specialist_generalist = balance_specialist_generalist,
        remove_desert_rocky_aa = remove_desert_rocky_aa,
        n_iterations = n_iterations,
        n_cores_used = n_cores
      )
    )
    
    saveRDS(all_results, paste0(dft_folder, "/allData.rds"))
    
    cat("\n=== SUMMARY OF RESULTS ===\n")
    cat("Total iterations completed:", length(iteration_results), "\n")
    cat("Results saved in:", dft_folder, "\n")
    cat("Main output file: allData.rds\n")
  }
  
  # return results
  return(list(
    ideam = btst_ideam,
    tert = tert_ideam,
    all_results = all_results
  ))
}

# ================== Run =====================
dft_folder <- "results/aoh_results_randomCI_final_7gen_2022_nobal"
if (!dir.exists(dft_folder)) {
  dir.create(dft_folder, recursive = TRUE)
}

# Example run with 1000 parallel iterations
example_result <- run_analysis(
  num_generalist = 7,
  d_near = 0,
  random_seed = 2025,  # Base seed (not used for iterations)
  balance_specialist_generalist = FALSE,
  remove_desert_rocky_aa = FALSE,  # Set to TRUE to remove desert and rocky habitats
  save_results = TRUE,
  dft_folder = dft_folder,
  n_cores = NULL  # Will use detectCores() - 3
)

# ================== Analysis and Visualization =====================
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

# Analyze the results
if (exists("example_result") && !is.null(example_result$all_results)) {
  cat("\n=== Analyzing Results ===\n")
  
  # Extract data
  extracted_data <- extract_habitat_data(example_result$all_results$iteration_results)
  
  # Create matrices
  count_matrix <- create_count_matrix(extracted_data$all_hab_data, 
                                     extracted_data$land_covers, 
                                     extracted_data$habitats)
  
  ci_matrix <- create_ci_matrix(extracted_data$all_raw_odds, 
                               extracted_data$land_covers, 
                               extracted_data$habitats)
  
  # Create CI bounds table
  ci_bounds_table <- create_ci_bounds_table(extracted_data$all_raw_odds, 
                                           extracted_data$land_covers, 
                                           extracted_data$habitats)
  
  # Save matrices
  write.csv(count_matrix, paste0(dft_folder, "/count_above_1_matrix.csv"))
  write.csv(ci_matrix, paste0(dft_folder, "/ci_width_matrix.csv"))
  write.csv(ci_bounds_table, paste0(dft_folder, "/ci_bounds_table.csv"))
  
  # Create and save barplot
  ci_plot <- create_ci_barplot(ci_matrix, dft_folder)
  
  # Create and save CI bounds heatmap
  ci_heatmap <- create_ci_bounds_heatmap(extracted_data$all_raw_odds, 
                                        extracted_data$land_covers, 
                                        extracted_data$habitats, 
                                        dft_folder)
  
  # Create and save count heatmap
  count_heatmap <- create_count_heatmap(count_matrix, dft_folder)
  
  # ci_heatmap
  # Print summary statistics
  cat("Count matrix dimensions:", dim(count_matrix), "\n")
  cat("CI matrix dimensions:", dim(ci_matrix), "\n")
  cat("Total pairs with values > 1:", sum(count_matrix > 0), "\n")
  cat("Average CI width:", mean(ci_matrix, na.rm = TRUE), "\n")
  
  cat("\nResults saved in:", dft_folder, "\n")
  cat("- count_above_1_matrix.csv: Number of iterations > 1 for each pair\n")
  cat("- ci_width_matrix.csv: 95% CI width for each pair\n")
  cat("- ci_bounds_table.csv: 95% CI bounds (lower, upper) for each pair\n")
  cat("- ci_width_barplot.png: Barplot of CI widths\n")
  cat("- ci_bounds_heatmap.png: Heatmap of CI bounds with significance coloring\n")
  cat("- count_above_1_heatmap.png: Heatmap of count of iterations > 1\n")
  
} else {
  cat("\nNo results found to analyze. Make sure to run the analysis first.\n")
}

gc()


tert_df <- read.csv(paste0(dft_folder, '/btst_ideam_randomCI_gen7__noDesRockAA_1104_2025_tertiles.csv'))
summary(tert_df$one.third)
summary(tert_df$two.third)
hist(tert_df$one.third)
hist(tert_df$two.third)
