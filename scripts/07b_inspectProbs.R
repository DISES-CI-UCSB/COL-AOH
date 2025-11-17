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
df_all <- read.csv('data/occ_pts/allinfo_ideam_cgls_coords_2022.csv') %>% select(-X)
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
remove_desert_rocky_aa = TRUE
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


# ======= 3 check Savanna (habitat) =======
df_all_sa <- df_all_info %>% filter(hab_2 == 1)
colnames(df_all_sa)
df_info_sa <- df_all_sa %>% select(c(species, hab_1, hab_3, hab_2, hab_4, hab_5, hab_14.12, hab_14.36, hab_14.45)) %>% unique()
colSums(df_info_sa[, sapply(df_info_sa, is.numeric)])

df_info_sa$sum <- rowSums(df_info_sa[,grep("^hab_", colnames(df_info_sa))])
table(df_info_sa$sum)

# check overall specialist species composition
df_all$sum <- rowSums(df_all[, grep("^hab_", colnames(df_all))])
df_sp <- df_all %>% filter(sum==1) %>% select(-c(X, occ.ID, nivel_1_num, longitude, nivel_2_num, nivel_3_num, longitude, latitude))
df_sp_unique <- df_sp %>% unique()
colSums(df_sp_unique[, sapply(df_sp_unique, is.numeric)])

# ======= 4 check Continental Waters (lc) =======
colnames(df_all_info)
df_ctwt <- df_all_info %>% filter(nvl_2_n == '51')
colnames(df_ctwt)

# Ensure longitude and latitude columns exist and are numeric
df_ctwt$longitude <- as.numeric(df_ctwt$longitude)
df_ctwt$latitude <- as.numeric(df_ctwt$latitude)

# Remove rows with missing coordinates
df_ctwt_sf <- df_ctwt %>% filter(!is.na(longitude) & !is.na(latitude))

# Convert to sf object (GeoDataFrame)
df_ctwt_sf <- sf::st_as_sf(df_ctwt_sf, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)

df_ctwt_sa_sf <- df_ctwt_sf %>% filter(hab_2==1)
df_ctwt_sa <- df_ctwt_sa_sf[c('species', 'taxa')]
df_ctwt_sa$geometry <- NULL
df_ctwt_sa <- unique(df_ctwt_sa)
table(df_ctwt_sa$taxa)

# Export as shapefile
dir.create("data/tmp")
sf::st_write(df_ctwt_sf, "data/tmp/continental_waters_points.shp", delete_layer = TRUE)
sf::st_write(df_ctwt_sa_sf, "data/tmp/continental_waters_sav_points.shp", delete_layer = TRUE)



ctwt_sp <- df_all_info %>% select(c(species, taxa, hab_1, hab_3, hab_2, hab_4, hab_5, hab_14.12, hab_14.36, hab_14.45)) %>% unique()
table(ctwt_sp$taxa)
colSums(ctwt_sp[, sapply(ctwt_sp, is.numeric)])/nrow(ctwt_sp)

prob <- ctwt_sp %>% filter(!hab_5==1)
table(prob$taxa)
prob_aves <- prob %>% filter(taxa == 'aves')
prob_aves_name <- unique(prob_aves$species)
head(prob_aves_name)

prob_2 <- prob %>% filter(taxa == 'mamiferos')
prob_2_name <- unique(prob_2$species)
head(prob_2_name)


