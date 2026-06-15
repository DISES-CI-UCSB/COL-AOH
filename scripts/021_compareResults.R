packages <- c("sf", "terra", "here", "rgrass", "dplyr", "data.table", 
              "rredlist", "tidyr", "Matrix", "tidyverse", "pROC",
              "kableExtra", "ggplot2", "reshape2", "parallel")
lapply(packages, library, character.only = TRUE)

setwd('/')
setwd('Users/wenxinyang/Desktop/GitHub/colander')
source('scripts/refineBiomodelos/funcs.R')
source('scripts/refineBiomodelos/info.R')

folder_gen7 <- 'results/glm_btst_2012_keephab/gen7_glm_btst_2012_keephab/'
cv_gen7 <- read.csv(file.path(folder_gen7, 'cv_table.csv'))
se_gen7 <- read.csv(file.path(folder_gen7, 'se_table.csv'))

folder_gen5 <- 'results/glm_btst_2012_keephab/gen8_glm_btst_2012_keephab/'
cv_gen5 <- read.csv(file.path(folder_gen5, 'cv_table.csv'))
se_gen5 <- read.csv(file.path(folder_gen5, 'se_table.csv'))

## ------------ compare cv -------------------
cv_gen7_1 <- as.matrix(cv_gen7 %>% select(-any_of("X")))
cv_gen5_1 <- as.matrix(cv_gen5 %>% select(-any_of("X")))

both_finite <- is.finite(cv_gen7_1) & is.finite(cv_gen5_1)
delta_cv <- matrix(0, nrow = nrow(cv_gen7_1), ncol = ncol(cv_gen7_1))
delta_cv[both_finite] <- cv_gen7_1[both_finite] - cv_gen5_1[both_finite]
cv_gen7_1[!both_finite] <- 0
cv_gen5_1[!both_finite] <- 0
sum(delta_cv)

## ------------ compare se -------------------
se_gen7_1 <- as.matrix(se_gen7 %>% select(-any_of("X")))
se_gen5_1 <- as.matrix(se_gen5 %>% select(-any_of("X")))

both_finite <- is.finite(se_gen7_1) & is.finite(se_gen5_1)
delta_se <- matrix(0, nrow = nrow(cv_gen7_1), ncol = ncol(se_gen7_1))
delta_se[both_finite] <- se_gen7_1[both_finite] - se_gen5_1[both_finite]
se_gen7_1[!both_finite] <- 0
se_gen5_1[!both_finite] <- 0
sum(delta_se)
