addGeneralistInfo <- function(df_pref, num_generalist){
  drop_codes <- c(7, 9, 10, 12, 13, 16, 17, 18)
  drop_codes <- paste0('hab_', drop_codes)
  df_pref <- df_pref %>% select(-all_of(drop_codes))
  # merge several artificial habitat types
  df_pref <- df_pref %>% mutate(
    hab_14.12 = ifelse(hab_14.1+hab_14.2 >0, 1, 0),
    hab_14.36 = ifelse(hab_14.3+hab_14.6 >0, 1, 0),
    hab_14.45 = ifelse(hab_14.4+hab_14.5 >0, 1, 0)
  )

  df_pref <- df_pref %>% select(-all_of(c('hab_14.1', 'hab_14.2', 'hab_14.3', 'hab_14.4', 'hab_14.5', 'hab_14.6', 'X')))
  # row sum of df_pref for all columns except for the first column
  df_pref$row_sum <- rowSums(df_pref[, -1])
  df_pref <- df_pref %>% mutate(
    type = ifelse(row_sum > num_generalist, 'generalist', ifelse(row_sum==1, 'specialist', 'other') )
  )
  df_basic_info <- df_pref[c('name', 'type')]

  return(df_basic_info)
}



create_unbalanced_splits <- function(lc_dataset, df_all_info, lc_code, seed=226){
  set.seed(seed)
  # create presence/absence dataset for the specified land cover
  if(lc_dataset == 'ideam') {
    df_lc <- df_all_info %>% mutate(
      presence = ifelse(nvl_2_n == lc_code, 1, 0)
    )
  } else if(lc_dataset == 'cgls') {
    df_lc <- df_all_info %>% mutate(
      presence = ifelse(cgls_value == lc_code, 1, 0)
    )
  }

  df_pres <- df_lc %>%
    filter(presence == 1) %>%
    select(presence, starts_with('hab_'))

  df_abs <- df_lc %>%
    filter(presence == 0) %>%
    select(presence, starts_with('hab_'))

  # split presence cases
  n_pres <- nrow(df_pres)
  train_size_pres <- round(n_pres * 0.7)
  test_size_pres <- round((n_pres - train_size_pres)) # use all rest for test set

  set.seed(seed + 2)  # different seed for presence splits
  pres_indices <- 1:n_pres
  train_idx_pres <- sample(pres_indices, size = train_size_pres)
  remaining_idx_pres <- setdiff(pres_indices, train_idx_pres)
  test_idx_pres <- sample(remaining_idx_pres, size = test_size_pres)
  val_idx_pres <- setdiff(remaining_idx_pres, test_idx_pres)

  # split absence cases (using same sizes as presence for balance)
  n_abs <- nrow(df_abs)
  train_size_abs <- train_size_pres
  test_size_abs <- test_size_pres

  set.seed(seed + 3)  # different seed for absence splits
  abs_indices <- 1:n_abs
  train_idx_abs <- sample(abs_indices, size = train_size_abs)
  remaining_idx_abs <- setdiff(abs_indices, train_idx_abs)
  test_idx_abs <- sample(remaining_idx_abs, size = test_size_abs)
  val_idx_abs <- setdiff(remaining_idx_abs, test_idx_abs)

  # combine the splits
  train_set <- rbind(
    df_pres[train_idx_pres, ],
    df_abs[train_idx_abs, ]
  )

  test_set <- rbind(
    df_pres[test_idx_pres, ],
    df_abs[test_idx_abs, ]
  )

  val_set <- rbind(
    df_pres[val_idx_pres, ],
    df_abs[val_idx_abs, ]
  )

  return(list(
    train = train_set,
    test = test_set,
    validation = val_set,
    n_samples = list(
      total = nrow(df_pres) * 2,
      train = nrow(train_set),
      test = nrow(test_set),
      validation = nrow(val_set)
    )
  ))

}


create_splits_replace <- function(lc_dataset, df_all_info, lc_code,
                                  seed=226){
  set.seed(seed)  # set seed at the start

  # 1. Create presence/absence dataset for the specified land cover
  if(lc_dataset == 'ideam') {
    if('nvl_2_n' %in% colnames(df_all_info)){
      df_lc <- df_all_info %>% mutate(
        presence = ifelse(nvl_2_n == lc_code, 1, 0))
    }
    if('nvl_1_n' %in% colnames(df_all_info)){
      df_lc <- df_all_info %>% mutate(
        presence = ifelse(nvl_1_n == lc_code, 1, 0))
    }
  } else if(lc_dataset == 'cgls') {
    df_lc <- df_all_info %>% mutate(
      presence = ifelse(cgls_value == lc_code, 1, 0)
    )
  }

  # 2. Filter presence cases
  df_pres <- df_lc %>%
    filter(presence == 1) %>%
    select(presence, starts_with('hab_'))

  n_pres <- nrow(df_pres)

  # 3. Sample absences to match presence count (Balanced Background Selection)
  set.seed(seed + 1)  # different seed for absence sampling
  df_abs <- df_lc %>%
    filter(presence == 0) %>%
    sample_n(n_pres) %>%
    select(presence, starts_with('hab_'))

  n_abs <- nrow(df_abs)

  # 4. Bootstrap Presences (WITH REPLACEMENT)
  set.seed(seed + 2)
  train_idx_pres <- sample(1:n_pres, size = n_pres, replace = TRUE)
  # Test set is purely the Out-of-Bag (leftover) unique rows
  test_idx_pres  <- setdiff(1:n_pres, unique(train_idx_pres))

  # 5. Bootstrap Absences (WITH REPLACEMENT)
  set.seed(seed + 3)
  train_idx_abs <- sample(1:n_abs, size = n_abs, replace = TRUE)
  # Test set is purely the Out-of-Bag (leftover) unique rows
  test_idx_abs  <- setdiff(1:n_abs, unique(train_idx_abs))

  # 6. Combine the Bootstrap (In-Bag) and Out-of-Bag (OOB) splits
  train_set <- rbind(
    df_pres[train_idx_pres, ],
    df_abs[train_idx_abs, ]
  )

  test_set <- rbind(
    df_pres[test_idx_pres, ],
    df_abs[test_idx_abs, ]
  )


  # TEST FOR PROBLEM
  n_test_pres <- sum(test_set$presence == 1)
  n_test_abs  <- sum(test_set$presence == 0)

  if (n_test_pres == 0 || n_test_abs == 0) {
    # If a parallel node hits this, return NULL instead of letting it crash
    return(NULL)
  }



  return(list(
    train = train_set,
    test = test_set,  # This is your official Out-of-Bag validation set
    n_samples = list(
      total_pool = n_pres * 2,
      train      = nrow(train_set), # Will equal total_pool, but with duplicated rows
      test_oob   = nrow(test_set)   # Independent unique testing rows (~36.8% of pool)
    )
  ))
}


create_balanced_splits <- function(lc_dataset, df_all_info, lc_code,
                                   seed = 226) {
  set.seed(seed)  # set seed at the start

  # create presence/absence dataset for the specified land cover
  if(lc_dataset == 'ideam') {
    if('nvl_2_n' %in% colnames(df_all_info)){
      df_lc <- df_all_info %>% mutate(
        presence = ifelse(nvl_2_n == lc_code, 1, 0))
    }
    if('nvl_1_n' %in% colnames(df_all_info)){
      df_lc <- df_all_info %>% mutate(
        presence = ifelse(nvl_1_n == lc_code, 1, 0))
    }
  } else if(lc_dataset == 'cgls') {
    df_lc <- df_all_info %>% mutate(
      presence = ifelse(cgls_value == lc_code, 1, 0)
    )
  }

  # create balanced presence/absence datasets
  df_pres <- df_lc %>%
    filter(presence == 1) %>%
    select(presence, starts_with('hab_'))

  set.seed(seed + 1)  # different seed for absence sampling
  df_abs <- df_lc %>%
    filter(presence == 0) %>%
    sample_n(nrow(df_pres)) %>%
    select(presence, starts_with('hab_'))

  # split presence cases
  n_pres <- nrow(df_pres)
  train_size_pres <- round(n_pres * 0.7)
  test_size_pres <- round((n_pres - train_size_pres)) # use all rest for test set

  set.seed(seed + 2)  # different seed for presence splits
  pres_indices <- 1:n_pres
  train_idx_pres <- sample(pres_indices, size = train_size_pres)
  remaining_idx_pres <- setdiff(pres_indices, train_idx_pres)
  test_idx_pres <- sample(remaining_idx_pres, size = test_size_pres)
  val_idx_pres <- setdiff(remaining_idx_pres, test_idx_pres)

  # split absence cases (using same sizes as presence for balance)
  n_abs <- nrow(df_abs)
  train_size_abs <- train_size_pres
  test_size_abs <- test_size_pres

  set.seed(seed + 3)  # different seed for absence splits
  abs_indices <- 1:n_abs
  train_idx_abs <- sample(abs_indices, size = train_size_abs)
  remaining_idx_abs <- setdiff(abs_indices, train_idx_abs)
  test_idx_abs <- sample(remaining_idx_abs, size = test_size_abs)
  val_idx_abs <- setdiff(remaining_idx_abs, test_idx_abs)

  # combine the splits
  train_set <- rbind(
    df_pres[train_idx_pres, ],
    df_abs[train_idx_abs, ]
  )

  test_set <- rbind(
    df_pres[test_idx_pres, ],
    df_abs[test_idx_abs, ]
  )

  val_set <- rbind(
    df_pres[val_idx_pres, ],
    df_abs[val_idx_abs, ]
  )

  return(list(
    train = train_set,
    test = test_set,
    validation = val_set,
    n_samples = list(
      total = nrow(df_pres) * 2,
      train = nrow(train_set),
      test = nrow(test_set),
      validation = nrow(val_set)
    )
  ))
}

build_evaluate_model <- function(splits, lc_code, modtype) {
  # access the splits
  train_data <- splits$train
  test_data <- splits$test

  if(modtype == 'glm') {  # build logistic regression model
    model <- glm(presence ~ ., data = train_data, family = 'binomial')
    test_data$predicted_prob <- predict(model, test_data, type = 'response')
  }

  if(modtype == 'firth'){ # firth's penalty
    library(logistf)
    model <- logistf(presence ~ ., data = train_data)
    test_data$predicted_prob <- as.numeric(predict(model, newdata = test_data, type = "response"))
  }

  # The classification step remains type-safe
  test_data$predicted_class <- ifelse(test_data$predicted_prob > 0.5, 1, 0)

  # compute and print AUC
  roc_obj <- roc(test_data$presence, test_data$predicted_prob)
  auc_value <- auc(roc_obj)
  # print(paste("AUC for land cover", lc_code, ":", round(auc_value, 3)))

  return(list(
    model = model,
    auc = auc_value,
    test_predictions = test_data,
    n_samples = nrow(train_data)
  ))
}

get_odds_ratios_row <- function(model, lc_code, auc_value, n_samples,
                                modtype, filterpval=0) {
  # get variable names from the model
  var_names <- names(coef(model))[-1]  # exclude intercept

  # initialize vectors for odds ratios and p-values
  odds_ratios <- rep(0, length(var_names))
  p_values <- rep(1, length(var_names))  # default to 1 (not significant)
  names(odds_ratios) <- var_names
  names(p_values) <- var_names

  # fill in values for non-NA coefficients
  valid_coefs <- which(!is.na(coef(model)[-1]))
  if(length(valid_coefs) > 0) {
    odds_ratios[valid_coefs] <- exp(coef(model)[-1][valid_coefs])
    if(modtype=='glm'){
      model_summary <- summary(model)$coefficients
      p_values[valid_coefs] <- model_summary[-1, "Pr(>|z|)"][valid_coefs]
    }
    if(modtype=='firth'){
      model_summary <- model$prob
      p_values[valid_coefs] <- model_summary[-1][valid_coefs]
    }

  }

  # round values
  odds_ratios <- round(odds_ratios, 2)
  p_values <- round(p_values, 3)

  # create table with odds ratios and p-values
  odds_ratios_table <- data.frame(
    odds_ratios = odds_ratios,
    p_values = p_values,
    row.names = var_names
  )

  if(filterpval==1){
    # set odds ratios to 0 where p-values > 0.05
    odds_ratios_table$odds_ratios[odds_ratios_table$p_values > 0.05] <- 0
  }

  # set column names and get land cover names
  colnames(odds_ratios_table) <- c(paste0("lc_", lc_code), "p_values")
  lc_names <- rownames(odds_ratios_table)

  # create final transposed table
  odds_ratios_table <- odds_ratios_table[, -2]  # remove p_values column
  odds_ratios_row <- as.data.frame(t(odds_ratios_table))
  colnames(odds_ratios_row) <- lc_names
  odds_ratios_row$lc_code <- paste0("lc_", lc_code)

  # attach AUC and n_samples to the final output
  odds_ratios_row$auc <- auc_value
  odds_ratios_row$n_samples <- n_samples

  return(odds_ratios_row)
}

get_unbalanced_single_row <- function(df, lc_code, seed = 226, lc_dataset, modtype, filterpval) {
  # print(lc_code)
  splits_lc <- create_unbalanced_splits(lc_dataset, df, lc_code, seed)
  model_results_lc <- build_evaluate_model(splits_lc, lc_code)
  odds_ratios_row_lc <- get_odds_ratios_row(model_results_lc$model,
                                            lc_code,
                                            model_results_lc$auc,
                                            model_results_lc$n_samples,
                                            modtype,
                                            filterpval)
  return(odds_ratios_row_lc)
}

get_a_single_row <- function(df, lc_code, seed = 226, lc_dataset, modtype,
                             filterpval, resampling_approach) {
  # print(lc_code)
  if(resampling_approach == 'bootstrap'){
    splits_lc <- create_splits_replace(lc_dataset, df, lc_code, seed)
    if(is.null(splits_lc)){
      tmpnames <- c(colnames(df %>% select(starts_with("hab_"))), "lc_code", "auc", "n_samples")
      tmp <- data.frame(matrix(NA, nrow=1, ncol=length(tmpnames)))
      colnames(tmp) <- tmpnames
      tmp$lc_code <- lc_code

      return(tmp)
      }
  }

  if(resampling_approach == 'subsampling'){
    splits_lc <- create_balanced_splits(lc_dataset, df, lc_code, seed)
  }

  model_results_lc <- build_evaluate_model(splits_lc, lc_code, modtype)
  odds_ratios_row_lc <- get_odds_ratios_row(model_results_lc$model,
                                            lc_code,
                                            model_results_lc$auc,
                                            model_results_lc$n_samples,
                                            modtype,
                                            filterpval)
  return(odds_ratios_row_lc)
}

# ======= create tables for visualization ======
format_habitat_table <- function(df, ifbtst) {
  # round all numeric columns to 2 decimal places
  df_rounded <- df %>%
    mutate(across(where(is.numeric), ~round(., 2)))

  # create color formatting for habitat columns
  hab_cols <- grep("^hab_", names(df_rounded), value = TRUE)

  if(ifbtst==0){
    val1 = 1.71
    val2 = 1.35
  } else {
    val1 = 0.8*ifbtst
    val2 = 0.5*ifbtst
  }
  # create style function for conditional formatting
  style_df <- df_rounded %>%
    mutate(across(3:(ncol(df_rounded)-1),
                  ~cell_spec(.,
                             background = case_when(
                               . > val1 ~ "#006400", # Dark green
                               . > val2 ~ "#90EE90", # Light green
                               TRUE ~ "white"
                             ),
                             color = ifelse(. > val2, "white", "black")
                             #color = "black"
                             )))

  # convert to kable format with styling
  formatted_table <- kable(style_df, format = "html", escape = FALSE) %>%
    kable_styling(bootstrap_options = c("striped", "hover"))

  return(formatted_table)
}

#' calculate distances from points to nearest polygon boundary
#' @param points sf object containing points
#' @param polygons sf object containing polygons
#' @return numeric vector of distances
calc_dist_to_boundary <- function(points, polygons) {
  # Convert polygons to boundary lines
  polygon_boundaries <- st_boundary(polygons)

  # Calculate distance from each point to nearest boundary
  distances <- st_distance(points, polygon_boundaries)

  # Get minimum distance for each point
  min_distances <- apply(distances, 1, min)

  return(min_distances)
}


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
  } else if(posOrHigh=='high'){
    df_sel <- df %>%
      select(land_cover, high_odds_habitats) %>%
      separate_rows(high_odds_habitats, sep=", ") %>%
      group_by(land_cover, high_odds_habitats) %>%
      summarise(count=n(), .groups='drop') %>%
      right_join(
        expand.grid(
          land_cover=unique(df$land_cover),
          high_odds_habitats=unique(unlist(strsplit(df$high_odds_habitats, ", ")))
        ),
        by=c("land_cover", "high_odds_habitats")
      ) %>%
      mutate(count = replace_na(count, 0)) %>%
      pivot_wider(names_from=high_odds_habitats, values_from=count, values_fill=0)
  }

  # add missing habitats
  all_habitats <- unique(c(
    unlist(strsplit(df$pos_odds_habitats, ", ")),
    unlist(strsplit(df$high_odds_habitats, ", ")),
    unlist(strsplit(df$medium_odds_habitats, ", "))
  ))

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
}
