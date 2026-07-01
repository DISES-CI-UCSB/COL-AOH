packages <- c("sf", "terra", "here", "rgrass", "dplyr", "data.table", 
              "rredlist", "tidyr", "Matrix", "tidyverse", "pROC",
              "kableExtra", "ggplot2", "reshape2", "parallel")
lapply(packages, library, character.only = TRUE)

setwd('/')
setwd('Users/wenxinyang/Desktop/GitHub/colander')
source('scripts/refineBiomodelos/funcs.R')
source('scripts/refineBiomodelos/info.R')

nums <- c(5, 6, 7, 8)
folders <- paste0('gen', nums, '_glm_btst_2012_keephab/')
cvmats <- c()
for(i in 1:length(folders)){
  f <- folders[i]
  folder <- file.path('results/glm_btst_2012_keephab', f)
  cv_f <- read.csv(file.path(folder, 'cv_table.csv')) %>% filter(X!='Mine and Waste Dump')
  cv_f_mat <- as.matrix(cv_f %>% select(-any_of("X")))
  cvmats[[i]] <- cv_f_mat
}

# use sum, mean, se

# ----------- get sum of cv for all ----------
cv_sums <- sapply(cvmats, sum, na.rm = TRUE)

# ----------- get median of cv for all ---------
cv_meds <- sapply(cvmats, median, na.rm=TRUE)

# ----------- get mean of cv for all ---------
cv_mean <- sapply(cvmats, mean, na.rm=TRUE)

# ----------- get sd of cv for all ---------
cv_sd <- sapply(cvmats, sd, na.rm=TRUE)

# ----------- get sd of cv for all ---------
se <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
cv_se <- sapply(cvmats, se)
# ----------- get min of cv for all ---------
cv_min <- sapply(cvmats, min, na.rm=TRUE)

# ----------- get max of cv for all ---------
cv_max <- sapply(cvmats, max, na.rm=TRUE)

# ----------- get ranks ---------------
nrow <- nrow(cvmats[[1]])
ncol <- ncol(cvmats[[1]])
nmat <- length(cvmats)
cvarray <- array(unlist(cvmats), dim = c(nrow, ncol, nmat))
rank_array <- array(NA, dim = dims)
dims <- dim(cvarray)
# For each cell (i,j), rank across matrices
for (i in 1:dims[1]) {
  for (j in 1:dims[2]) {
    # Get values at position (i,j) from all matrices
    values <- cvarray[i, j, ]
    # Rank them and store back
    rank_array[i, j, ] <- rank(values, ties.method = "average")
  }
}

rank_mats <- lapply(1:dims[3], function(m) rank_array[, , m])

cv_ranks <- sapply(rank_mats, sum, na.rm=TRUE)

# ------------------ combine them into a single dataframe -------------------
cv_info <- data.frame(
  N_gen = nums,
  cv_sums = cv_sums,
  cv_min = cv_min,
  #cv_meds = cv_meds,
  cv_mean = cv_mean,
  cv_se = cv_se,
  cv_max = cv_max,
  cv_sd = cv_sd,
  cv_ranks = cv_ranks
)


