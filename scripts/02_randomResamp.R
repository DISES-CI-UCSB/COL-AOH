# Matrix calculation
# with random resampling and statistics
# Author: Wenxin Yang
# Date: July, 2025
# Revised from 02_randomCI.R in: June, 2026
# THIS IS THE SCRIPT I ENDED UP USING FOR THE MANUSCRIPT (revision June)

# ================== Function to run single matrix analysis with parameters =====================
run_analysis <- function(
    # core parameters
  num_generalist = 7,
  random_seed = 2025,
  # data filtering parameters
  balance_specialist_generalist = 1, # balance specialist vs non specialist by taxa
  remove_habitats = NULL, # c("hab_6", "hab_8") for removing desert (hab_6), rocky (hab_8)
  # output parameters
  output_suffix = paste0("_", format(Sys.Date(), "%m%d"), "_", format(Sys.Date(), "%Y")),
  save_results = TRUE,
  resampling_approach = 'bootstrap', # bootstrap or subsampling
  dft_folder = "aoh_results_randomCI",
  modtype = "glm", # glm or firth
  # parallel processing parameters
  n_cores = NULL
) {
  
  # ================== Prep =====================
  # load libraries
  packages <- c("sf", "terra", "here", "rgrass", "dplyr", "data.table", 
                "rredlist", "tidyr", "Matrix", "tidyverse", "pROC",
                "kableExtra", "ggplot2", "reshape2", "parallel", "logistf")
  lapply(packages, library, character.only = TRUE)
  
  setwd('/')
  setwd('Users/wenxinyang/Desktop/GitHub/colander')
  source('scripts/refineBiomodelos/funcs.R')
  source('scripts/refineBiomodelos/postProcFuncs.R')
  source('scripts/refineBiomodelos/info.R')
  
  # create aoh_results directory structure if it doesn't exist
  if (!dir.exists(dft_folder)) dir.create(dft_folder)
  if (!dir.exists(paste0(dft_folder, "/viz"))) dir.create(paste0(dft_folder, "/viz"))
  
  
  # ======= 2. Read in & prep data =======
  df_all <- read.csv('data/occ_pts/allinfo_ideam_coords_2012_0605.csv') %>% select(-any_of(c('nvl_1_n', 'nivel_1_num')))
  if('scientificName' %in% colnames(df_all)){
    df_all <- df_all %>% mutate(species = scientificName) %>% select(-scientificName)
  }
  if('scntfcN' %in% colnames(df_all)){
    df_all <- df_all %>% mutate(species = scntfcN) %>% select(-scntfcN)
  }
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
  
  # merge several artificial habitat types
  df_all_info <- df_all_info %>% mutate(
    hab_14.12 = ifelse(hab_14.1+hab_14.2 >0, 1, 0),
    hab_14.36 = ifelse(hab_14.3+hab_14.6 >0, 1, 0),
    hab_14.45 = ifelse(hab_14.4+hab_14.5 >0, 1, 0)
  )
  
  # remove the original ones
  df_all_info <- df_all_info %>% select(-all_of(c('hab_14.1', 'hab_14.2', 'hab_14.3', 'hab_14.4', 'hab_14.5', 'hab_14.6')))
  
  # ======= 2.3 Remove some habitats if requested =======
  if (!is.null(remove_habitats)) {
    df_all_info <- df_all_info %>% select(-all_of(remove_habitats))
    # Update colname_dataset to exclude these habitats
    colname_dataset <- c("lc_code", "n_samples", habitat_info1$habitat_code[!habitat_info1$habitat_code %in% remove_habitats], "auc")
  } else {
    colname_dataset <- c("lc_code", "n_samples", habitat_info1$habitat_code, "auc")
  }
  
  # remove rows with NA values
  df_all_info <- df_all_info %>% select(-any_of(c('year', 'eventYr', 'finalYr','id')))
  # tmp <- df_all_info[!complete.cases(df_all_info),]
  df_all_info <- df_all_info[complete.cases(df_all_info), ]
  sum(is.na(df_all_info))
  df_all_info$X <- NULL
  
  # sum all columns that start with "hab_"
  df_all_info$sum <- rowSums(df_all_info[,grep("^hab_", colnames(df_all_info))])
  df_all_info <- df_all_info %>% filter(sum>0) %>% select(-sum)
  if('scntfcN' %in% colnames(df_all_info)){
    df_all_info$species <- df_all_info$scntfcN
    df_all_info$scntfcN <- NULL
  }
  
  # colSums(df_all_info[,grep("^hab_", colnames(df_all_info))])
  
  ## ======= 2.1 remove generalist species =======
  li_generalist <- unique((df_basic_info %>% filter(type == 'generalist'))$name)
  if('species' %in% colnames(df_all_info)){
    df_all_info <- df_all_info %>% filter(!species %in% li_generalist)
  } else if('scntfcN' %in% colnames(df_all_info)){
    df_all_info <- df_all_info %>% filter(!scntfcN %in% li_generalist)
  }
  
  
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
                      "create_splits_replace",
                      "build_evaluate_model", "get_odds_ratios_row",
                      "resampling_approach"))
  
  # Load required packages and define functions on cluster
  clusterEvalQ(cl, {
    library(dplyr)
    library(tidyr)
    library(pROC)
    library(stringr)
    
    # Define generate_matrix function on cluster
    generate_matrix <- function(dat, landcover_colname, landcover_dataset, seed, cols, modtype){
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
        #unique_landcovers <- as.numeric(unique_landcovers)
        
        if (length(unique_landcovers) == 0) {
          stop("No valid land cover values found")
        }
        
        info_list <- lapply(unique_landcovers,
                            function(x) {
                              tryCatch({
                                get_a_single_row(dat, x, seed = seed, 
                                                 lc_dataset = landcover_dataset, 
                                                 modtype=modtype, filterpval=0,
                                                 resampling_approach=resampling_approach)
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
      if (nrow(df) == 0) {
        return(data.frame(
          land_cover = character(0),
          pos_odds_habitats = character(0),
          seed_val = numeric(0),
          dataset = character(0),
          auc = numeric(0),
          n_samples = numeric(0),
          stringsAsFactors = FALSE
        ))
      }
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
  run_single_iteration <- function(seed, df_all_info, df_basic_info, colname_dataset, balance_specialist_generalist, modtype) {
    cat(seed)
    # seed, df_all_info, df_basic_info, colname_dataset, balance_specialist_generalist
    tryCatch({
      # Balance specialist vs non-specialist within each taxa for this iteration
      if (balance_specialist_generalist) {
        set.seed(seed)  # Use the iteration seed
        spp_field <- ifelse('species' %in% colnames(df_all_info), 'species', 'scntfcN')
        records <- merge(df_all_info, df_basic_info, by.x=spp_field, by.y='name')
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
        spp_field <- ifelse('species' %in% colnames(df_all_info), 'species', 'scntfcN')
        balanced_rec <- merge(df_all_info, df_basic_info, by.x=spp_field, by.y='name')
      }
      
      # Generate matrix for this iteration
      dat <- balanced_rec
      
      # Make sure data structure is okay
      if (nrow(dat) == 0) {
        stop("No data available for iteration with seed ", seed)
      }
      
      if (!'nvl_2_n' %in% colnames(dat)) {
        stop("Land cover column 'nvl_2_n' not found in data for iteration with seed ", seed)
      }
      
      this_ideam <- generate_matrix(dat, 'nvl_2_n', 'ideam', seed, colname_dataset, modtype)
      
      all_values <- as.vector(unlist(this_ideam %>% select(-`Land Cover`, auc, n_samples)))
      all_values <- all_values[all_values > 1] # keep positive values
      
      if (length(all_values) == 0) {
        # If no positive values, create empty results
        this_tert <- data.frame(
          min = NA,`one-third` = NA,`two-third` = NA, max = NA, seed = seed
        )
        tmp <- data.frame(
          land_cover = character(0),
          pos_odds_habitats = character(0),
          seed_val = numeric(0),
          dataset = character(0),
          auc = numeric(0),
          n_samples = numeric(0),
          stringsAsFactors = FALSE
        )
        
        # Create empty raw odds ratios dataframe
        raw_odds <- data.frame(
          land_cover = character(0),
          habitat = character(0),
          odds_ratio = numeric(0),
          seed = numeric(0),
          stringsAsFactors = FALSE
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
        land_cover = character(0),
        pos_odds_habitats = character(0),
        seed_val = numeric(0),
        dataset = character(0),
        auc = numeric(0),
        n_samples = numeric(0),
        stringsAsFactors = FALSE
      )
      
      # Create empty raw odds ratios dataframe
      raw_odds <- data.frame(
        land_cover = character(0),
        habitat = character(0),
        odds_ratio = numeric(0),
        seed = numeric(0),
        stringsAsFactors = FALSE
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
    run_single_iteration(seed, df_all_info, df_basic_info, colname_dataset, balance_specialist_generalist, modtype)
  })
  
  ### ---------------- debug codes ----------------
#  iteration_results <- lapply(1:1000, function(seed) {
#    run_single_iteration(seed, df_all_info, df_basic_info, colname_dataset, balance_specialist_generalist, modtype)
#  })
  ### --------------------- debug codes end -------------------  
  
  cat("Parallel iterations complete.\n")
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
                            ifelse(!is.null(remove_habitats), "rmhab", ""),
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
        balance_specialist_generalist = balance_specialist_generalist,
        remove_habitats = remove_habitats,
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


# ======================= Analyze & Visualize Results =========================
save_results <- function(example_result){
  if (!is.null(example_result$all_results)) {
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
    
    se_table <- create_stability_metrics_tables(extracted_data$all_raw_odds, 
                                                extracted_data$land_covers, 
                                                extracted_data$habitats, 
                                                'bootstrap_se')
    cv_table <- create_stability_metrics_tables(extracted_data$all_raw_odds, 
                                                extracted_data$land_covers, 
                                                extracted_data$habitats, 
                                                'cv')
    mean_table <- create_stability_metrics_tables(extracted_data$all_raw_odds, 
                                                  extracted_data$land_covers, 
                                                  extracted_data$habitats, 
                                                  'mean')
    med_table <- create_stability_metrics_tables(extracted_data$all_raw_odds, 
                                                 extracted_data$land_covers, 
                                                 extracted_data$habitats, 
                                                 'median')
    
    # Save matrices
    write.csv(count_matrix, paste0(dft_folder, "/count_above_1_matrix.csv"))
    write.csv(ci_matrix, paste0(dft_folder, "/ci_width_matrix.csv"))
    write.csv(ci_bounds_table, paste0(dft_folder, "/ci_bounds_table.csv"))
    write.csv(se_table, paste0(dft_folder, "/se_table.csv"))
    write.csv(cv_table, paste0(dft_folder, "/cv_table.csv"))
    write.csv(mean_table, paste0(dft_folder, "/mean_table.csv"))
    write.csv(med_table, paste0(dft_folder, "/median_table.csv"))
    
    # Create and save barplot
    #ci_plot <- create_ci_barplot(ci_matrix, dft_folder)
    
    # Create and save CI bounds heatmap
    #ci_heatmap <- create_ci_bounds_heatmap(extracted_data$all_raw_odds, 
    #                                       extracted_data$land_covers, 
    #                                       extracted_data$habitats, 
    #                                       dft_folder)
    
    # Create and save count heatmap
    count_heatmap <- create_count_heatmap(count_matrix, dft_folder)
    ggsave(paste0(dft_folder, "/viz/count_heatmap.png"), count_heatmap, width=12, height=8)
    
  } else {
    cat("\nNo results found to analyze. Make sure to run the analysis first.\n")
  }
}


# ================== Run =====================
overall_folder <- 'results/glm_btst_2012_rmhab'
if (!dir.exists(overall_folder)) {
  dir.create(overall_folder, recursive = TRUE)
}

for(i in c(5,6,7,8)){
  cat('Generalist spp has ', i, ' preferences\n')
  # genN modtype resampling_approach year_spp_filter remove_habitat
  dft_folder <- file.path(overall_folder, paste0('gen', i, '_glm_btst_2012_rmhab'))
  if (!dir.exists(dft_folder)) {
    dir.create(dft_folder, recursive = TRUE)
  }
  
  # Example run with 1000 parallel iterations
  example_result <- run_analysis(
    num_generalist = i,
    random_seed = 2025,  # Base seed (not used for iterations)
    balance_specialist_generalist = 1,
    remove_habitats = c('hab_6', 'hab_8'),  # Set to TRUE to remove desert and rocky habitats
    save_results = TRUE,
    modtype = 'glm',
    dft_folder = dft_folder,
    resampling_approach='bootstrap',
    n_cores = NULL  # Will use detectCores() - 3
  )
  
  save_results(example_result)
  
}




