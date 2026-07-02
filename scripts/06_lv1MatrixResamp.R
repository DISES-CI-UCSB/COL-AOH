# Script to compute translation matrix for level 1 land cover
# Author: Wenxin Yang
# Date: July, 2025
# Modified: June, 2026

# ================== Function to run single matrix analysis with parameters =====================
run_analysis <- function(
    # core parameters
  num_generalist = 7,
  random_seed = 2025,
  # data filtering parameters
  balance_specialist_generalist = 1, # balance specialist vs non specialist by taxa
  remove_habitats = NULL, # e.g., c("hab_6", "hab_8") for removing desert and rocky
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
  df_all <- read.csv('data/occ_pts/allinfo_ideam_coords_2012_0605.csv')
  
  if('scientificName' %in% colnames(df_all)){
    df_all <- df_all %>% mutate(species = scientificName) %>% select(-scientificName)
  }
  if('scntfcN' %in% colnames(df_all)){
    df_all <- df_all %>% mutate(species = scntfcN) %>% select(-scntfcN)
  }
  
  if(!'nvl_1_n' %in% colnames(df_all)){
    df_all$nvl_1_n <- df_all$nivel_1_num
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
  
  # remove rows with NA values and unneeded cols
  df_all_info <- df_all_info %>% select(-any_of(c('year', 'eventYr', 'finalYr','id')))
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
  n_cores <- min(n_cores, detectCores() - 1, 1000)
  
  cat("Using", n_cores, "cores for parallel processing\n")
  cl <- makeCluster(n_cores)
  
  # Export necessary data and functions to cluster
  # FIX: Added ideam_lc_info_lv1 back in so the parallel workers know the names
  clusterExport(cl, c("habitat_info1", "ideam_lc_info", "ideam_lc_info_lv1", "cgls_lc_info", 
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
    generate_matrix <- function(dat, landcover_colname, landcover_dataset, seed, cols, modtype, resampling_approach){
      tryCatch({
        if (!landcover_colname %in% colnames(dat)) {
          stop("Land cover column '", landcover_colname, "' not found in data")
        }
        
        landcover_values <- dat[[landcover_colname]]
        if (!is.character(landcover_values)) {
          landcover_values <- as.character(landcover_values)
        }
        unique_landcovers <- unique(landcover_values)
        unique_landcovers <- unique_landcovers[!is.na(unique_landcovers) & unique_landcovers != ""]
        
        if (length(unique_landcovers) == 0) {
          stop("No valid land cover values found")
        }
        
        info_list <- lapply(unique_landcovers,
                            function(x) {
                              tryCatch({
                                get_a_single_row(dat, x, seed = seed, 
                                                 lc_dataset = landcover_dataset, 
                                                 modtype = modtype, filterpval = 0,
                                                 resampling_approach = resampling_approach)
                              }, error = function(e) {
                                cat("Error processing land cover", x, ":", e$message, "\n")
                                empty_row <- data.frame(matrix(0, nrow = 1, ncol = length(cols)))
                                colnames(empty_row) <- cols
                                empty_row$lc_code <- paste0("lc_", x)
                                empty_row$auc <- 0.5
                                empty_row$n_samples <- 0
                                return(empty_row)
                              })
                            })
        
        info_list <- info_list[!sapply(info_list, is.null)]
        
        if (length(info_list) == 0) {
          empty_matrix <- data.frame(matrix(0, nrow = 0, ncol = length(cols)))
          colnames(empty_matrix) <- cols
          return(empty_matrix)
        }
        
        info_list <- do.call(rbind, info_list) %>% select(all_of(cols))
        rownames(info_list) <- NULL
        
        colnames(info_list) <- sapply(colnames(info_list), function(x){
          if (x %in% habitat_info1$habitat_code) {
            return(habitat_info1$habitat_name[match(x, habitat_info1$habitat_code)])
          } else {
            return(x)
          }
        })
        
        # FIX: Restored the level 1/level 2 distinction for naming
        info_list$lc_code <- sapply(info_list$lc_code, function(x) {
          if(landcover_dataset=='ideam' & landcover_colname == 'nvl_2_n'){
            return(ideam_lc_info$ideam_lc_name[match(x, ideam_lc_info$ideam_lc_code)])
          } else if(landcover_dataset=='ideam' & landcover_colname == 'nvl_1_n'){
            return(ideam_lc_info_lv1$ideam_lc_name[match(x, ideam_lc_info_lv1$ideam_lc_code)])
          } else if(landcover_dataset=='cgls'){
            return(cgls_lc_info$cgls_lc_name[match(x, cgls_lc_info$cgls_lc_code)])
          }
        })
        
        colnames(info_list)[colnames(info_list)=='lc_code'] <- 'Land Cover'
        return(info_list)
      }, error = function(e) {
        cat("Error in generate_matrix:", e$message, "\n")
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
      
      habitat_associations <- data.frame(
        land_cover = df$`Land Cover`,
        pos_odds_habitats = sapply(1:nrow(df), function(i) {
          habitat_cols <- colnames(df)[!colnames(df) %in% c("Land Cover", "n_samples", "auc")]
          odds_ratios <- as.numeric(df[i, habitat_cols])
          pos_odds <- habitat_cols[odds_ratios > 1]
          paste(pos_odds, collapse=", ")
        })
      )
      habitat_associations$seed_val <- seed
      habitat_associations$dataset <- landcover_dataset
      habitat_associations$auc <- df$auc
      habitat_associations$n_samples <- df$n_samples
      
      return(habitat_associations)
    }
  })
  
  # Function to run a single iteration with a specific seed
  run_single_iteration <- function(seed, df_all_info, df_basic_info, colname_dataset, balance_specialist_generalist, modtype, resampling_approach) {
    tryCatch({
      if (balance_specialist_generalist > 0) {
        set.seed(seed)
        spp_field <- ifelse('species' %in% colnames(df_all_info), 'species', 'scntfcN')
        records <- merge(df_all_info, df_basic_info, by.x=spp_field, by.y='name')
        balanced_rec <- data.frame(matrix(nrow=0, ncol=ncol(records)))
        colnames(balanced_rec) <- colnames(records)
        
        for (t in unique(records$taxa)) {
          sub_records <- records %>% filter(taxa == t)
          set.seed(seed)
          
          # FIX: Now properly using numsamp in sample_n
          numsamp <- balance_specialist_generalist * sum(sub_records$type == 'specialist')
          
          sub_records_other <- sub_records %>% 
            filter(type == 'other') %>% 
            sample_n(numsamp)
          
          sub_records <- rbind(sub_records %>% filter(type != 'other'), sub_records_other)
          balanced_rec <- rbind(balanced_rec, sub_records)
        }
        balanced_rec <- balanced_rec %>% distinct()
      } else {
        spp_field <- ifelse('species' %in% colnames(df_all_info), 'species', 'scntfcN')
        balanced_rec <- merge(df_all_info, df_basic_info, by.x=spp_field, by.y='name')
      }
      
      dat <- balanced_rec
      
      if (nrow(dat) == 0) {
        stop("No data available for iteration with seed ", seed)
      }
      
      if (!'nvl_1_n' %in% colnames(dat)) {
        stop("Land cover column 'nvl_1_n' not found in data for iteration with seed ", seed)
      }
      
      this_ideam <- generate_matrix(dat, 'nvl_1_n', 'ideam', seed, colname_dataset, modtype, resampling_approach)
      
      all_values <- as.vector(unlist(this_ideam %>% select(-`Land Cover`, auc, n_samples)))
      all_values <- all_values[all_values > 1]
      
      if (length(all_values) == 0) {
        this_tert <- data.frame(min = NA, `one-third` = NA, `two-third` = NA, max = NA, seed = seed)
        tmp <- data.frame(land_cover = character(0), pos_odds_habitats = character(0), seed_val = numeric(0), dataset = character(0), auc = numeric(0), n_samples = numeric(0), stringsAsFactors = FALSE)
        raw_odds <- data.frame(land_cover = character(0), habitat = character(0), odds_ratio = numeric(0), seed = numeric(0), stringsAsFactors = FALSE)
      } else {
        this_tert <- data.frame(
          min = min(all_values),
          `one-third` = quantile(all_values, 1/3),
          `two-third` = quantile(all_values, 2/3), 
          max = max(all_values),
          seed = seed
        )
        
        tmp <- getHabAssoc(this_ideam, seed, 'ideam')
        raw_odds <- data.frame()
        
        for (i in 1:nrow(this_ideam)) {
          land_cover <- this_ideam$`Land Cover`[i]
          habitat_cols <- colnames(this_ideam)[!colnames(this_ideam) %in% c("Land Cover", "n_samples", "auc")]
          for (hab in habitat_cols) {
            odds_ratio <- this_ideam[i, hab]
            if (!is.na(odds_ratio)) {
              raw_odds <- rbind(raw_odds, data.frame(land_cover = land_cover, habitat = hab, odds_ratio = odds_ratio, seed = seed))
            }
          }
        }
      }
      
      gc()
      return(list(tert = this_tert, hab_assoc = tmp, raw_odds = raw_odds))
    }, error = function(e) {
      cat("Error in iteration with seed", seed, ":", e$message, "\n")
      this_tert <- data.frame(min = NA, `one-third` = NA, `two-third` = NA, max = NA, seed = seed)
      tmp <- data.frame(land_cover = character(0), pos_odds_habitats = character(0), seed_val = numeric(0), dataset = character(0), auc = numeric(0), n_samples = numeric(0), stringsAsFactors = FALSE)
      raw_odds <- data.frame(land_cover = character(0), habitat = character(0), odds_ratio = numeric(0), seed = numeric(0), stringsAsFactors = FALSE)
      return(list(tert = this_tert, hab_assoc = tmp, raw_odds = raw_odds))
    })
  }
  
  selectPairs <- function(df, posOrHigh){
    if(posOrHigh=='positive'){
      df_sel <- df %>%
        select(land_cover, pos_odds_habitats) %>%
        separate_rows(pos_odds_habitats, sep=", ") %>%
        group_by(land_cover, pos_odds_habitats) %>%
        summarise(count=n(), .groups='drop') %>%
        right_join(
          expand.grid(land_cover=unique(df$land_cover), pos_odds_habitats=unique(unlist(strsplit(df$pos_odds_habitats, ", ")))),
          by=c("land_cover", "pos_odds_habitats")
        ) %>%
        mutate(count = replace_na(count, 0)) %>%
        pivot_wider(names_from=pos_odds_habitats, values_from=count, values_fill=0)
      
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
  
  n_iterations <- 1000
  cat("Running", n_iterations, "iterations in parallel...\n")
  
  iteration_results <- parLapply(cl, 1:n_iterations, function(seed) {
    run_single_iteration(seed, df_all_info, df_basic_info, colname_dataset, balance_specialist_generalist, modtype, resampling_approach)
  })
  
  cat("Parallel iterations complete.\n")
  stopCluster(cl)
  
  btst_ideam <- do.call(rbind, lapply(iteration_results, function(x) x$hab_assoc))
  tert_ideam <- do.call(rbind, lapply(iteration_results, function(x) x$tert))
  raw_odds_all <- do.call(rbind, lapply(iteration_results, function(x) x$raw_odds))
  
  # ======= 6. Process and save results =======
  if (save_results) {
    filename_base <- paste0(dft_folder, "/btst_ideam_randomCI_", 
                            "gen", num_generalist, "_",
                            ifelse(!is.null(remove_habitats), "rmhab", ""),
                            ifelse(!balance_specialist_generalist, "_nobal", ""),
                            output_suffix)
    
    write.csv(btst_ideam, paste0(filename_base, ".csv"))
    write.csv(tert_ideam, paste0(filename_base, "_tertiles.csv"))
    
    btst_ideam_pos <- selectPairs(btst_ideam, 'positive')
    write.csv(btst_ideam_pos, paste0(filename_base, "_pos.csv"))
    
    all_results <- list(
      iteration_results = iteration_results,
      combined_hab_assoc = btst_ideam,
      combined_tertiles = tert_ideam,
      raw_odds_ratios = raw_odds_all,
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
  
  return(list(
    ideam = btst_ideam,
    tert = tert_ideam,
    all_results = all_results
  ))
}

# ================== Analysis and Visualization Functions =====================
extract_habitat_data <- function(iteration_results) {
  all_hab_data <- do.call(rbind, lapply(iteration_results, function(x) x$hab_assoc))
  all_raw_odds <- do.call(rbind, lapply(iteration_results, function(x) x$raw_odds))
  land_covers <- unique(all_hab_data$land_cover)
  habitats <- unique(all_raw_odds$habitat)
  return(list(all_hab_data = all_hab_data, all_raw_odds = all_raw_odds, land_covers = land_covers, habitats = habitats))
}

create_count_matrix <- function(all_hab_data, land_covers, habitats) {
  count_matrix <- matrix(0, nrow = length(land_covers), ncol = length(habitats))
  rownames(count_matrix) <- land_covers
  colnames(count_matrix) <- habitats
  for (lc in land_covers) {
    lc_data <- all_hab_data[all_hab_data$land_cover == lc, ]
    for (hab in habitats) {
      count_above_1 <- sum(grepl(hab, lc_data$pos_odds_habitats, fixed = TRUE), na.rm = TRUE)
      count_matrix[lc, hab] <- count_above_1
    }
  }
  return(count_matrix)
}

create_ci_matrix <- function(all_raw_odds, land_covers, habitats) {
  ci_matrix <- matrix(NA, nrow = length(land_covers), ncol = length(habitats))
  rownames(ci_matrix) <- land_covers
  colnames(ci_matrix) <- habitats
  for (lc in land_covers) {
    for (hab in habitats) {
      pair_data <- all_raw_odds[all_raw_odds$land_cover == lc & all_raw_odds$habitat == hab, ]
      values <- pair_data$odds_ratio
      values <- values[!is.na(values)]
      if (length(values) > 0) {
        ci <- quantile(values, c(0.025, 0.975))
        ci_matrix[lc, hab] <- ci[2] - ci[1]
      }
    }
  }
  return(ci_matrix)
}

create_ci_bounds_table <- function(all_raw_odds, land_covers, habitats) {
  ci_bounds_matrix <- matrix(NA, nrow = length(land_covers), ncol = length(habitats))
  rownames(ci_bounds_matrix) <- land_covers
  colnames(ci_bounds_matrix) <- habitats
  for (lc in land_covers) {
    for (hab in habitats) {
      pair_data <- all_raw_odds[all_raw_odds$land_cover == lc & all_raw_odds$habitat == hab, ]
      values <- pair_data$odds_ratio
      values <- values[!is.na(values)]
      if (length(values) > 0) {
        ci <- quantile(values, c(0.025, 0.975))
        ci_bounds_matrix[lc, hab] <- paste0("(", round(ci[1], 3), ", ", round(ci[2], 3), ")")
      } else {
        ci_bounds_matrix[lc, hab] <- "(NA, NA)"
      }
    }
  }
  return(ci_bounds_matrix)
}

create_count_heatmap <- function(count_matrix, dft_folder) {
  pos_file <- list.files(dft_folder, pattern = "_pos.csv", full.names = TRUE)
  if (length(pos_file) == 0) stop("No _pos.csv file found in the specified folder")
  
  pos_data <- read.csv(pos_file[1])
  if("X" %in% colnames(pos_data)) pos_data$X <- NULL
  
  habitat_cols <- colnames(pos_data)[!colnames(pos_data) %in% c("", "land_cover", "auc", "n_samples")]
  
  pos_long <- pos_data %>%
    select(land_cover, all_of(habitat_cols)) %>%
    tidyr::gather(habitat, count_value, -land_cover)
  
  pos_long$color_category <- cut(pos_long$count_value, 
                                 breaks = c(-Inf, 300, 600, 900, Inf),
                                 labels = c("Not a pair (≤300)", "High (300-600)", "Medium (600-900)", "Low (>900)"),
                                 include.lowest = TRUE)
  
  auc_data <- data.frame(land_cover = pos_data$land_cover, habitat = "AUC", count_value = pos_data$auc, color_category = "AUC")
  n_samples_data <- data.frame(land_cover = pos_data$land_cover, habitat = "n_samples", count_value = round(pos_data$n_samples), color_category = "n_samples")
  
  plot_data <- rbind(pos_long, auc_data, n_samples_data)
  plot_data$color_category <- factor(plot_data$color_category, levels = c("Low (>900)", "Medium (600-900)", "High (300-600)", "Not a pair (≤300)", "AUC", "n_samples"))
  plot_data$land_cover <- factor(plot_data$land_cover, levels = sort(unique(plot_data$land_cover), decreasing = TRUE))
  plot_data$habitat <- factor(plot_data$habitat, levels = c(habitat_cols, "AUC", "n_samples"))
  
  p <- ggplot2::ggplot(plot_data, aes(x = habitat, y = land_cover, fill = color_category)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = ifelse(color_category == "AUC", sprintf("%.2f", count_value), as.character(count_value))), 
              size = 5, fontface = "bold", color = "black", family = "Times New Roman") +
    scale_fill_manual(values = c("Not a pair (≤300)" = "#f7f7f7", "High (300-600)" = "#fee8c8", "Medium (600-900)" = "#fdbb84", "Low (>900)" = "#e34a33", "AUC" = "white", "n_samples" = "white"), name = "Uncertainty level (counts)") +
    scale_x_discrete(position = "top", labels = function(x) { sapply(x, function(label) { if (nchar(label) > 15) { paste(strwrap(label, width=15), collapse="\n") } else { label } }) }) +
    theme_minimal() +
    theme(text = element_text(family = "Times New Roman"), axis.text.x.top = element_text(size = 12, hjust = 0.5), axis.text.y = element_text(size = 15), axis.title.x = element_text(size = 16, face = 'bold'), axis.title.y = element_text(size = 16, face = 'bold'), legend.title = element_text(size = 16), legend.text = element_text(size = 14), legend.position = "right", plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), panel.grid = element_blank(), panel.background = element_rect(fill = "white", color = NA)) +
    labs(x = "Habitat Classes", y = "Land Cover Classes") +
    coord_fixed(ratio = 0.6)
  
  ggplot2::ggsave(paste0(dft_folder, "/viz/count_heatmap.png"), p, width = 16, height = 14, dpi = 300)
  return(p)
}

# ================== Run =====================
overall_folder <- "results/btst_2012_lv1"
if (!dir.exists(overall_folder)) dir.create(overall_folder, recursive = TRUE)

resampling_approach <- 'bootstrap'
keephab <- 1
removehab <- if(keephab == 1) NULL else c('hab_6', 'hab_8')

i <- 7
dft_folder <- file.path(overall_folder, paste0('gen', i, '_glm_', resampling_approach, '_2012_keephab'))
if (!dir.exists(dft_folder)) dir.create(dft_folder, recursive = TRUE)

example_result <- run_analysis(
  num_generalist = i,
  random_seed = 2025,
  balance_specialist_generalist = 1,
  remove_habitats = removehab, 
  save_results = TRUE,
  modtype = 'glm',
  dft_folder = dft_folder,
  resampling_approach = resampling_approach,
  n_cores = NULL 
)

# ======================= Analyze & Visualize Results =========================
if (!is.null(example_result$all_results)) {
  cat("\n=== Analyzing Results ===\n")
  
  extracted_data <- extract_habitat_data(example_result$all_results$iteration_results)
  
  count_matrix <- create_count_matrix(extracted_data$all_hab_data, extracted_data$land_covers, extracted_data$habitats)
  ci_matrix <- create_ci_matrix(extracted_data$all_raw_odds, extracted_data$land_covers, extracted_data$habitats)
  ci_bounds_table <- create_ci_bounds_table(extracted_data$all_raw_odds, extracted_data$land_covers, extracted_data$habitats)
  
  se_table <- create_stability_metrics_tables(extracted_data$all_raw_odds, extracted_data$land_covers, extracted_data$habitats, 'bootstrap_se')
  cv_table <- create_stability_metrics_tables(extracted_data$all_raw_odds, extracted_data$land_covers, extracted_data$habitats, 'cv')
  mean_table <- create_stability_metrics_tables(extracted_data$all_raw_odds, extracted_data$land_covers, extracted_data$habitats, 'mean')
  med_table <- create_stability_metrics_tables(extracted_data$all_raw_odds, extracted_data$land_covers, extracted_data$habitats, 'median')
  
  write.csv(count_matrix, paste0(dft_folder, "/count_above_1_matrix.csv"))
  write.csv(ci_matrix, paste0(dft_folder, "/ci_width_matrix.csv"))
  write.csv(ci_bounds_table, paste0(dft_folder, "/ci_bounds_table.csv"))
  write.csv(se_table, paste0(dft_folder, "/se_table.csv"))
  write.csv(cv_table, paste0(dft_folder, "/cv_table.csv"))
  write.csv(mean_table, paste0(dft_folder, "/mean_table.csv"))
  write.csv(med_table, paste0(dft_folder, "/median_table.csv"))
  
  count_heatmap <- create_count_heatmap(count_matrix, dft_folder)
  
  cat("Count matrix dimensions:", dim(count_matrix), "\n")
  cat("CI matrix dimensions:", dim(ci_matrix), "\n")
  cat("Total pairs with values > 1:", sum(count_matrix > 0), "\n")
  cat("Average CI width:", mean(ci_matrix, na.rm = TRUE), "\n")
  cat("\nResults saved in:", dft_folder, "\n")
  
} else {
  cat("\nNo results found to analyze. Make sure to run the analysis first.\n")
}

gc()

mean(cv_table)
max(cv_table)
