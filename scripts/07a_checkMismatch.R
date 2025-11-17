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
df_basic_info <- addGeneralistInfo(df_pref, num_generalist=5)
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
remove_desert_rocky_aa = FALSE
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

# ======== 3 get validated aoh data ==========
validated <- read.csv('results/validation/aoh_validation_300m.csv')
validated$error_message <- NULL
validated <- validated[!is.na(validated$aoh_area_km2),] %>% 
  mutate(ifOkay = point_prevalence >= model_prevalence)
table(validated$taxa)
na_validated <- validated[is.na(validated$ifOkay), ]
validated <- validated[!is.na(validated$ifOkay),]
table(validated$taxa)

okay <- validated[validated$ifOkay==1,]
okay <- okay[!is.na(okay$species),]
table(okay$taxa)

li_spp_validaoh <- unique(validated$species)
li_spp_okay <- unique(okay$species)

# ======== 4 compare ==========
li_spp_all <- unique(df_all_info$species) # 1735

getPercAOH <- function(hab_class, validOrOkay){
  df_tmp_info <- df_all_info[df_all_info[hab_class] == 1, ]
  li_tmp <- unique(df_tmp_info$species)
  
  if(validOrOkay == 'valid'){
    li_valid <- intersect(li_tmp, li_spp_validaoh)
    print(length(li_valid)/length(li_tmp))
  }
  if(validOrOkay == 'okay'){
    li_okay <- intersect(li_tmp, li_spp_okay)
    print(length(li_okay)/length(li_tmp))
  }
  
}

for(i in 1:nrow(habitat_info1)){
  code <- habitat_info1[i, ]$habitat_code
  name <- habitat_info1[i, ]$habitat_name
  if(!name %in% c('Rocky Areas', 'Desert', 'Artificial-Aquatic')){
    print(name)
    getPercAOH(code, 'valid')
    getPercAOH(code, 'okay')
  }
}
