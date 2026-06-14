packages <- c("sf", "terra", "here", "rgrass", "dplyr", "data.table", 
              "rredlist", "tidyr", "Matrix", "tidyverse", "pROC",
              "kableExtra", "ggplot2", "reshape2", "parallel")
lapply(packages, library, character.only = TRUE)

setwd('/')
setwd('Users/wenxinyang/Desktop/GitHub/colander')
source('scripts/refineBiomodelos/funcs.R')
source('scripts/refineBiomodelos/info.R')

folder_gen7 <- 'results/7gen_glm_pval_oob_2012_2022_bal_keepaa_spthi/'
count_gen7 <- read.csv(file.path(folder_gen7, 'count_above_1_matrix.csv'))
ci_gen7 <- read.csv(file.path(folder_gen7, 'ci_bounds_table.csv'))
cv_gen7 <- read.csv(file.path(folder_gen7, 'cv_table.csv'))
se_gen7 <- read.csv(file.path(folder_gen7, 'se_table.csv'))



folder_gen5 <- 'results/5gen_glm_pval_oob_2012_2022_bal_keepaa_spthi/'
count_gen5 <- read.csv(file.path(folder_gen5, 'count_above_1_matrix.csv'))
ci_gen5 <- read.csv(file.path(folder_gen5, 'ci_bounds_table.csv'))
cv_gen5 <- read.csv(file.path(folder_gen5, 'cv_table.csv'))
se_gen5 <- read.csv(file.path(folder_gen5, 'se_table.csv'))

rowSums(cv_gen5 %>% select(-X))

## ------------ compare cv -------------------
cv_gen7[is.na(cv_gen7)] <- 0
cv_gen5[is.na(cv_gen5)] <- 0

cv_gen7_1 <- as.matrix(cv_gen7 %>% select(-X))
cv_gen5_1 <- as.matrix(cv_gen5 %>% select(-X))

delta_cv <- cv_gen7_1 - cv_gen5_1
sum(delta_cv)

## ------------ compare se -------------------
se_gen7[is.na(se_gen7)] <- 0
se_gen5[is.na(se_gen5)] <- 0

se_gen7_1 <- as.matrix(se_gen7 %>% select(-X))
se_gen5_1 <- as.matrix(se_gen5 %>% select(-X))

delta_se <- se_gen7_1 - se_gen5_1
sum(delta_se)

sum(delta_se[1:14,])
