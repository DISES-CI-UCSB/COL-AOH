packages <- c("sf", "terra", "here", "rgrass", "dplyr", "data.table", 
              "rredlist", "tidyr", "Matrix", "tidyverse", "pROC",
              "kableExtra", "ggplot2", "reshape2", "parallel")
lapply(packages, library, character.only = TRUE)

setwd('/')
setwd('Users/wenxinyang/Desktop/GitHub/colander')
source('scripts/refineBiomodelos/funcs.R')
source('scripts/refineBiomodelos/info.R')

count_gen7 <- read.csv('results/aoh_results_randomCI_7gen_2012_2022_bal_rmaa_spthi/count_above_1_matrix.csv')
ci_gen7 <- read.csv('results/aoh_results_randomCI_7gen_2012_2022_bal_keepaa_spthi/ci_bounds_table.csv')

count_gen6 <- read.csv('results/aoh_results_randomCI_6gen_2012_2022_bal_rmaa_spthi/count_above_1_matrix.csv') 
ci_gen6 <- read.csv('results/aoh_results_randomCI_6gen_2012_2022_bal_keepaa_spthi/ci_bounds_table.csv')

count_gen5 <- read.csv('results/aoh_results_randomCI_5gen_2012_2022_bal_rmaa_spthi/count_above_1_matrix.csv')
ci_gen5 <- read.csv('results/aoh_results_randomCI_5gen_2012_2022_bal_keepaa_spthi/ci_bounds_table.csv')


sum7 <- sum(as.matrix(count_gen7 %>% select(-X)))
sum6 <- sum(as.matrix(count_gen6 %>% select(-X)))
sum5 <- sum(as.matrix(count_gen5 %>% select(-X)))
cat(sum7, ' ', sum6, ' ', sum5)

habsum_gen7 <- colSums(count_gen7 %>% select(-X))
habsum_gen6 <- colSums(count_gen6 %>% select(-X))
habsum_gen5 <- colSums(count_gen5 %>% select(-X))

habsum_info <- data.frame()
habsum_info <- rbind(habsum_info, habsum_gen7)
habsum_info <- rbind(habsum_info, habsum_gen6)
habsum_info <- rbind(habsum_info, habsum_gen5)
colnames(habsum_info) <- colnames(count_gen7 %>% select(-X))
rownames(habsum_info) <- 8-(1:nrow(habsum_info))
