# Parameterized version of 02_computeMatrix.R for bootstrapping 
# Author: Wenxin Yang
# Date: July, 2025

# ================== Function to run analysis with parameters =====================
run_analysis <- function(
    # core parameters
  num_generalist = 5,
  d_near = 0,
  df_seed = 2025,
  bal_seed = 0,
  n_bootstrap = 1000,
  
  # data filtering parameters
  remove_desert_savanna = FALSE,
  balance_specialist_generalist = TRUE, # balance specialist vs non specialist
  
  # output parameters
  output_suffix = paste0("_", format(Sys.Date(), "%m%d"), "_", format(Sys.Date(), "%Y")),
  save_results = TRUE,
  dft_folder = "aoh_results"
) {
  
  # ================== Prep =====================
  # load libraries
  packages <- c("sf", "terra", "here", "rgrass", "dplyr", "data.table", 
                "rredlist", "tidyr", "Matrix", "tidyverse", "pROC",
                "kableExtra", "ggplot2", "reshape2")
  lapply(packages, library, character.only = TRUE)
  
  setwd('/')
  setwd('Users/wenxinyang/Desktop/GitHub/colander')
  source('scripts/refineBiomodelos/funcs.R')
  source('scripts/refineBiomodelos/info.R')
  
  # create aoh_results directory structure if it doesn't exist
  if (!dir.exists(dft_folder)) dir.create(dft_folder)
  if (!dir.exists(paste0(dft_folder, "/viz"))) dir.create(paste0(dft_folder, "/viz"))
  
  cat("Running analysis with parameters:\n")
  cat("# of associations for a generalist:", num_generalist, "\n")
  cat("near distance:", d_near, "\n")
  cat("random seed:", df_seed, "\n")
  cat("balance seed:", bal_seed, "\n") # for balancing # of pts for specialist vs non specialist
  cat("number of bootstraps for computing the matrix:", n_bootstrap, "\n")
  cat("remove desert savanna:", remove_desert_savanna, "\n")
  cat("balance specialist generalist:", balance_specialist_generalist, "\n")
  cat("output suffix:", output_suffix, "\n")
  
  # ======= 2. Read in & prep data =======
  df_all <- read.csv('data/occ_pts/allinfo_ideam_cgls.csv')
  
  # add species generalist vs specialist info
  df_pref <- read.csv('data/occ_pts/col_animal_pref_cleaned.csv')
  df_basic_info <- addGeneralistInfo(df_pref, num_generalist)
  rm(df_pref)
  
  # remove habitat codes not included in the analysis per the original paper
  df_all_info <- df_all %>% select(-all_of(drop_cols))
  
  # remove those near boundary
  summary(df_all_info$dst_t_b)
  df_all_info <- df_all_info %>% filter(dst_t_b > d_near) %>% select(-dst_t_b)
  
  # merge several artificial habitat types
  df_all_info <- df_all_info %>% mutate(
    hab_14.12 = ifelse(hab_14.1+hab_14.2 >0, 1, 0),
    hab_14.36 = ifelse(hab_14.3+hab_14.6 >0, 1, 0),
    hab_14.45 = ifelse(hab_14.4+hab_14.5 >0, 1, 0)
  )
  
  # remove the original ones
  df_all_info <- df_all_info %>% select(-all_of(c('hab_14.1', 'hab_14.2', 'hab_14.3', 'hab_14.4', 'hab_14.5', 'hab_14.6')))
  
  # remove rows with NA values
  df_all_info <- df_all_info[complete.cases(df_all_info), ]
  sum(is.na(df_all_info))
  df_all_info$X <- NULL
  
  # sum all columns that start with "hab_"
  df_all_info$sum <- rowSums(df_all_info[,grep("^hab_", colnames(df_all_info))])
  df_all_info <- df_all_info %>% filter(sum>0) %>% select(-sum)
  
  ## ======= 2.1 remove generalist species =======
  li_generalist <- unique((df_basic_info %>% filter(type == 'generalist'))$name)
  df_all_info <- df_all_info %>% filter(!species %in% li_generalist)
  
  ## ======= 2.2 balance # of records specialist vs non-specialist ========
  if (balance_specialist_generalist) {
    set.seed(bal_seed)
    records <- merge(df_all_info, df_basic_info, by.x='species', by.y='name')
    new_rec <- data.frame(matrix(nrow=0, ncol=ncol(records)))
    colnames(new_rec) <- colnames(records)
    for (t in unique(records$taxa)){
      print(t)
      sub_records <- records %>% filter(taxa == t)
      print(table(sub_records$type))
      set.seed(bal_seed)
      sub_records_other <- sub_records %>% filter(type == 'other') %>% sample_n(sum(sub_records$type == 'specialist'))
      sub_records <- rbind(sub_records %>% filter(type != 'other'), sub_records_other)
      new_rec <- rbind(new_rec, sub_records)
    }
  } else {
    new_rec <- merge(df_all_info, df_basic_info, by.x='species', by.y='name')
  }
  
  # ======= 3. NO: Remove desert and savanna =======
  if (remove_desert_savanna) {
    cat("Removing desert and savanna habitats\n")
    new_rec <- new_rec %>% select(-all_of(c('hab_2', 'hab_8')))
    colname_dataset <- c("lc_code", "n_samples", habitat_info1$habitat_code[!habitat_info1$habitat_code %in% c('hab_2', 'hab_8')], "auc")
  } else {
    colname_dataset <- c("lc_code", "n_samples", habitat_info1$habitat_code, "auc")
  }
  
  # ======= 4. Compute the transition matrix =======
  generate_matrix <- function(dat, landcover_colname, landcover_dataset, seed, cols){
    info_list <- lapply(unlist(unique(dat[landcover_colname])),
                        function(x) get_a_single_row(dat, x,
                                                     seed = seed,
                                                     lc_dataset = landcover_dataset))
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
  }
  
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
        # Get habitats with high odds ratios
        pos_odds <- habitat_cols[odds_ratios > 1]
        paste(pos_odds, collapse=", ")
      }),
      high_odds_habitats = sapply(1:nrow(df), function(i) {
        # Get habitat columns (exclude irrelevant columns)
        habitat_cols <- colnames(df)[!colnames(df) %in% c("Land Cover", "n_samples", "auc")]
        # Get odds ratios for this land cover type
        odds_ratios <- as.numeric(df[i, habitat_cols])
        # Get habitats with high odds ratios
        high_odds <- habitat_cols[odds_ratios > 1.71]
        paste(high_odds, collapse=", ")
      }),
      medium_odds_habitats = sapply(1:nrow(df), function(i) {
        habitat_cols <- colnames(df)[!colnames(df) %in% c("Land Cover", "n_samples", "auc")] 
        odds_ratios <- as.numeric(df[i, habitat_cols])
        # Get habitats with medium odds ratios
        medium_odds <- habitat_cols[odds_ratios > 1.35 & odds_ratios <= 1.71]
        paste(medium_odds, collapse=", ")
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
  
  ## ========= 4.1 bootstrap =======
  btst_ideam <- as.data.frame(matrix(ncol=8, nrow=0))
  colnames(btst_ideam) <- c('land_cover', 'n_samples', 'positive_habitats', 'high_odds_habitats', 'medium_odds_habitats','auc', 'seed', 'dataset')
  
  tert_ideam <- as.data.frame(matrix(ncol=5, nrow=0))
  colnames(tert_ideam) <- c('min', 'one-third', 'two-third', 'max', 'seed')
  
  set.seed(df_seed)
  li_seeds <- sample(1:1000000, n_bootstrap)
  
  cat("computing matrix with ", n_bootstrap, "iterations...\n")
  for(i in 1:length(li_seeds)){
    if(i %% 100 == 0) cat("Completed", i, "iterations\n")
    
    this_seed = li_seeds[i]
    
    dat <- new_rec
    
    this_ideam <- generate_matrix(dat, 'nvl_2_n', 'ideam', this_seed, colname_dataset)
    
    all_values <- as.vector(unlist(this_ideam %>% select(-`Land Cover`, auc, n_samples)))
    all_values <- all_values[all_values > 1] # keep positive values
    
    this_tert <- data.frame(
      min = min(all_values),
      `one-third` = quantile(all_values, 1/3),
      `two-third` = quantile(all_values, 2/3), 
      max = max(all_values),
      seed = this_seed
    )
    tert_ideam <- rbind(tert_ideam, this_tert)
    
    tmp <- getHabAssoc(this_ideam, this_seed, 'ideam')
    btst_ideam <- rbind(btst_ideam, tmp)
  }
  
  # ======= 5. Process and save results =======
  if (save_results) {
    # Create filename with parameters
    filename_base <- paste0(dft_folder, "/btst_ideam_", 
                            bal_seed, "_",
                            "gen", num_generalist, "_",
                            ifelse(d_near > 0, paste0("near", d_near, "_"), ""),
                            "boot", n_bootstrap,
                            ifelse(remove_desert_savanna, "_noDS", ""),
                            ifelse(!balance_specialist_generalist, "_nobal", ""),
                            output_suffix)
    
    write.csv(btst_ideam, paste0(filename_base, ".csv"))
    write.csv(tert_ideam, paste0(filename_base, "_tertiles.csv"))
    
    # process bootstrap results for visualization
    # count positive odds ratios pairs
    btst_ideam_pos <- selectPairs(btst_ideam, 'positive') 
    
    # count high odds ratios pairs
    btst_ideam_high <- selectPairs(btst_ideam, 'high')
    
    # save formatted tables
    #pos_graph <- format_habitat_table(btst_ideam_pos, n_bootstrap)
    #save_kable(pos_graph, file = paste0("aoh_results/viz/btst_ideam_positive_", bal_seed, "_", n_bootstrap, output_suffix, ".html"))
    
    #high_graph <- format_habitat_table(btst_ideam_high, n_bootstrap)
    #save_kable(high_graph, file = paste0("aoh_results/viz/btst_ideam_high_", bal_seed, "_", n_bootstrap, output_suffix, ".html"))
  }
  
  # return results
  return(list(
    pos = btst_ideam_pos,
    tert = tert_ideam
  ))
}

# ================== Actual run =====================
dft_folder <- 'aoh_testNumGen'
parameter_combinations <- list()
for (i in seq(3, 6, 1)){
  parameter_combinations[[i-2]] <- list(num_generalist = i, bal_seed = 2025, d_near = 0, n_bootstrap = 1000,
                                      dft_folder = dft_folder)
}

all_results <- list()
for(i in seq_along(parameter_combinations)) {
  all_results[[i]] <- do.call(run_analysis, parameter_combinations[[i]])
}

# save all_results
save(all_results, file = paste0(dft_folder,"/all_results.RData"))

graph3 <- format_habitat_table(all_results[[1]]$pos, 1000)
save_kable(graph3, file = paste0(dft_folder, "/viz/btst_ideam_positive_gen_3_1000.html"))

graph4 <- format_habitat_table(all_results[[2]]$pos, 1000)
save_kable(graph4, file = paste0(dft_folder, "/viz/btst_ideam_positive_gen_4_1000.html"))

graph5 <- format_habitat_table(all_results[[3]]$pos, 1000)
save_kable(graph5, file = paste0(dft_folder, "/viz/btst_ideam_positive_gen_5_1000.html"))

graph6 <- format_habitat_table(all_results[[4]]$pos, 1000)
save_kable(graph6, file = paste0(dft_folder, "/viz/btst_ideam_positive_gen_6_1000.html"))