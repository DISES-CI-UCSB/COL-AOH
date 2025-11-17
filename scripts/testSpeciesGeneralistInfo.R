# Add specialist vs generalist info
# Author: Wenxin Yang
# Date: April, 2025

# Load libraries
packages <- c("sf", "terra", "here", "rgrass", "dplyr", "data.table", 
              "rredlist", "tidyr", "Matrix", "tidyverse", "pROC",
              "kableExtra")
lapply(packages, library, character.only = TRUE)

df_pref <- read.csv('data/occ_pts/col_animal_pref_cleaned.csv')

colnames(df_pref)

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



n_generalist <- 3
df_test <- df_pref %>% mutate(
  type = ifelse(row_sum>n_generalist, 'generalist', ifelse(row_sum==1, 'specialist', 'other') )
)

table(df_test$type)

gen_spp <- df_test %>% filter(type=='generalist')
num_gen_spp <- gen_spp %>% select(where(is.numeric))
colSums(num_gen_spp)/nrow(gen_spp)


# write.csv(df_pref, 'data/occ_pts/col_animal_pref_cleaned_type.csv')

