# Prepare species occurrence data & match /w land cover data
# Author: Wenxin Yang
# Date: April, 2025

# ================== Prep =====================
# Load libraries
packages <- c("sf", "terra", "here", "rgrass", "dplyr", "data.table", 
              "rredlist", "tidyr", "Matrix", "tidyverse", "pROC")
lapply(packages, library, character.only = TRUE)

setwd('/')
setwd('Users/wenxinyang/Desktop/GitHub/colander')
source('scripts/refineBiomodelos/funcs.R')
source('scripts/refineBiomodelos/info.R')

dfanf <- read.csv('data/occ_pts/anfibios_occ_pts.csv')
length(unique(dfanf$scientificName))
dfanf$taxa <- 'anfibios'
dfavs <- read.csv('data/occ_pts/aves_occ_pts.csv')
length(unique(dfavs$species))
dfavs$taxa <- 'aves'
dfmam <- read.csv('data/occ_pts/mamiferos_occ_pts.csv')
length(unique(dfmam$species))
dfmam$taxa <- 'mamiferos'
dfrep <- read.csv('data/occ_pts/squamata_occ_pts.csv')
length(unique(dfrep$species))
dfrep$taxa <- 'squamata'
colnames(dfanf) <- colnames(dfavs)

# ================= 1. Merge occurrence point files ==================
full_list <- rbind(dfanf, dfavs, dfmam, dfrep)
nrow(full_list) == nrow(dfanf) + nrow(dfavs) + nrow(dfmam) + nrow(dfrep)

names_li <- unique(full_list$species)


# check for spatial coverage criteria
# count species with < 10 pts
occ_count <- as.data.frame(table(full_list$species))
nrow(occ_count[occ_count$Freq < 10,])/nrow(occ_count)
# ~ 10% sepcies have < 10 points, we are not removing them
perc_occ <- 0.02*nrow(full_list)
nrow(occ_count[occ_count$Freq > perc_occ,])
# no species have more than 2% occ points

# ================= 2. Join points with species preference file ================ 

pref_path <- file.path('data/occ_pts/', 'animals_preference.csv')
df_pref <- read.csv(pref_path, sep=';')

names_li_2 <- unique(df_pref$name)

getlv1 <- function(x) {return(strsplit(as.character(x), '[.]')[[1]][1])}

# set NA elevations as 0 - 9000 m
pref_file <- df_pref %>%
  mutate(elevation_lower = ifelse(is.na(elevation_lower), 0, elevation_lower)) %>%
  mutate(elevation_upper = ifelse(is.na(elevation_upper), 9000, elevation_upper))

# extract level 1 iucn habitat
pref_file$iucn_code_lv1 <- lapply(pref_file$habitats.code, getlv1)

# extract IUCN habitat classes to match based on supp material (Table 1) from Lumbierres et al. 2021
pref_file <- pref_file %>% 
  mutate(iucn_code = ifelse(iucn_code_lv1 == '14', habitats.code, iucn_code_lv1)) %>%
  mutate(iucn_code = unlist(iucn_code))

# get iucn_codes from species preference information
iucn_hab_info <- unique(pref_file[c('habitats.code', 'iucn_code')]) %>% 
  mutate(iucn_code = unlist(iucn_code))

# convert to a different format
transform_preferences <- function(raw_pref_file) {
  
  # extract unique species names and habitat codes
  species_names <- unique(raw_pref_file$name)
  habitat_codes <- unique(raw_pref_file$iucn_code)
  habitat_codes <- habitat_codes[!is.na(habitat_codes)]
  
  # create an empty matrix with species as rows and habitat codes as columns
  result_matrix <- matrix(0, 
                          nrow = length(species_names), 
                          ncol = length(habitat_codes))
  
  # convert to data frame and set row/column names
  result_df <- as.data.frame(result_matrix)
  rownames(result_df) <- species_names
  colnames(result_df) <- habitat_codes
  
  # fill the matrix with binary values
  for (i in 1:nrow(raw_pref_file)) {
    species <- raw_pref_file$name[i]
    habitat <- raw_pref_file$iucn_code[i]
    if(!is.na(habitat)){
      result_df[species, habitat] <- 1
    }
  }
  
  result_df$name <- rownames(result_df)
  
  # reorder columns to put name first
  result_df <- result_df %>%
    select(name, everything())
  
  rownames(result_df) <- seq_len(nrow(result_df))
  
  return(result_df)
}

df_pref_sp <- transform_preferences(pref_file)
# write.csv(df_pref_sp, 'data/occ_pts/col_animal_pref_cleaned.csv')
# rename column names with prefix 'hab_' except for name
colnames(df_pref_sp) <- c('name', paste0('hab_', colnames(df_pref_sp)[-1]))

# ================= 3. Extract land cover values at pts ==================
# ============= 3.1 IDEAM land cover data =============
# read in land cover data
# geometry was fixed using QGIS
# lc <- read_sf('data/IDEAM_landcover_2018/fixed.shp') # year 2018
lc <- read_sf('data/Corine_hab_COL/Cobertura_tierra_100K_periodo_2022_limite_administrativo/ECOSISTEMAS_18062025/ECOSISTEMAS_18062025.gpkg') # year 2022
# plot(lc['nivel_2'])

colnames(lc)
head(lc)

lc <- lc %>% mutate(leyenda_num = sapply(strsplit(leyenda, ' '), `[`, 1),
                    nivel_1_num = sapply(strsplit(nivel_1, ' '), `[`, 1),
                    nivel_3_num = sapply(strsplit(nivel_3, ' '), `[`, 1),
                    nivel_2_num = sapply(strsplit(nivel_2, ' '), `[`, 1))

# ========== side task 1: get low confiabili polygons ========
# first get area info for all level 2 classes
#lc_nivel2 <- lc %>% select(nivel_2_num, Shape_Area)
#lc_nivel2$geometry <- NULL
#lc_nivel2_area <- lc_nivel2 %>% group_by(nivel_2_num) %>% summarise(sum_area = sum(Shape_Area))
#rm(lc_nivel2)

# a <- unique(lc$nivel_2_num)

# learn about polygons with low confiability
#lc_no <- lc %>% filter(confiabili == "NO")
#lc_no1 <- lc_no %>% select(nivel_2_num, Shape_Area)
#lc_no1$geometry <- NULL
#lc_no1_area<- lc_no1 %>% group_by(nivel_2_num) %>% summarise(low_area = sum(Shape_Area))
#rm(lc_no)

#low_confi <- merge(lc_nivel2_area, lc_no1_area, all.x=TRUE)
# fill na with 0
#low_confi$low_area <- ifelse(is.na(low_confi$low_area), 0, low_confi$low_area)
# get the percentage of low confiability polygons
#low_confi$perc_low <- 100*low_confi$low_area / low_confi$sum_area

#100*sum(low_confi$low_area)/sum(low_confi$sum_area)

# ========== side task completed ======
lc_1 <- lc %>% select(leyenda_num, nivel_1_num, nivel_3_num, nivel_2_num)
# , geometry
# extract land cover types for pt values

# 1. convert full_list into a point sf object
# first clean and validate coordinates
full_list <- full_list %>%
  # remove rows with NA coordinates
  filter(!is.na(decimalLongitude) & !is.na(decimalLatitude)) %>%
  # filter to only valid coordinate ranges
  filter(decimalLongitude >= -180 & decimalLongitude <= 180) %>%
  filter(decimalLatitude >= -90 & decimalLatitude <= 90)

# convert to sf with error checking
all_pts <- full_list %>% 
    st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), 
             crs = 4326) %>%
  st_transform(crs = st_crs(lc_1))
# export all_pts
# st_write(all_pts, 'data/occ_pts/all_pts.shp')
  
# ---------------- side task 2: compute distance to nearest boundary --------------
# dissolve lc_1 by nivel_1_num
#lc_1_dissolve <- lc_1 %>% group_by(nivel_1_num) %>% summarise()
lc_1_dissolve <- st_read('data/IDEAM_landcover_2018/level1.shp')
#st_write(lc_1_dissolve, 'data/IDEAM_landcover_2018/level1.shp')

unique(lc_1$nivel_1_num)
#lc_2_dissolve <- lc %>% group_by(nivel_2_num) %>% summarise()
#st_write(lc_2_dissolve, 'data/IDEAM_landcover_2018/level2.shp')


# check crs before running the analysis
crs(all_pts)
crs(lc_1_dissolve)
crs(all_pts) == crs(lc_1_dissolve)

# compute distance from each point to lc_1_dissolve
all_pts_lc_lv1 <- calc_dist_to_boundary(all_pts, lc_1_dissolve)
all_pts$dist_to_boundary <- all_pts_lc_lv1
# st_write(all_pts, 'data/occ_pts/all_pts_dist_to_boundary.shp')
# remove pts outside of boundary 180 -180 and -90 90
coords <- st_coordinates(all_pts)
all_pts <- all_pts[coords[,1] >= -180 & coords[,1] <= 180 & 
                   coords[,2] >= -90 & coords[,2] <= 90,]

head(all_pts)
write.csv(all_pts, 'data/occ_pts/all_pts_dist_to_boundary_lv1.csv')
# ---------------- side task completed ------------------
# 2. spatial join all_pts with lc_1
# First make sure geometries are valid
lc_1 <- st_make_valid(lc_1)
all_pts <- st_make_valid(all_pts)

# Perform the spatial join with error handling
all_pts_lc <- try({
  st_join(all_pts, lc_1, join = st_intersects)
})

if(inherits(all_pts_lc, "try-error")) {
  stop("Error in spatial join. Please check geometry validity.")
}
# export all_pts_lc
# st_write(all_pts_lc, 'data/occ_pts/all_pts_lc.shp')
# keep only rows with non na leyenda_num
all_pts_lc_col <- all_pts_lc %>% filter(!is.na(leyenda_num))
# export all_pts_lc_col
summary(all_pts_lc_col$dist_to_boundary)
st_write(all_pts_lc_col, 'data/occ_pts/all_pts_lc_col.shp', append = FALSE)

# ============= 3.2 CGLS dataset =============
# read in all_pts_lc_col
all_pts_lc_col <- vect('data/occ_pts/all_pts_lc_col.shp')

# read in CGLS dataset
cgls <- rast('data/PROBAV_LC100_global_v3.0.1_2015-base_Discrete-Classification-map_EPSG-4326.tif')

# convert the crs of all_pts_lc_col to the crs of cgls
all_pts_lc_col <- project(all_pts_lc_col, cgls)

# get the values of cgls at the points of all_pts_lc_col
# and rename column names to "cgls_value"
cgls_pts <- terra::extract(cgls, all_pts_lc_col)
colnames(cgls_pts) <- c('ID', 'value')

# add the values to all_pts_lc_col
all_pts_lc_col$cgls_value <- cgls_pts$value

# export all_pts_lc_col with terra package into shapefile
terra::writeVector(all_pts_lc_col, 'data/occ_pts/all_pts_lc_col_cgls.shp', overwrite=TRUE)

# ================= 4. Merge and format data ==================
all_pts_lc_col <- st_read('data/occ_pts/all_pts_lc_col_cgls.shp')

# merge all_pts_lc_col with df_pref_sp
all_pts_lc_col_pref <- merge(all_pts_lc_col, df_pref_sp, by.x = 'species', by.y = 'name', all.x = TRUE)
# remove geogmetry of all_pts_lc_col_pref
df_all_info <- all_pts_lc_col_pref
# Add longitude and latitude columns from geometry
coords <- st_coordinates(df_all_info$geometry)
df_all_info$longitude <- coords[,1]
df_all_info$latitude <- coords[,2]

df_all_info$geometry <- NULL
write.csv(df_all_info, 'data/occ_pts/allinfo_ideam_cgls_coords_2022.csv')

gc()


summary(df_all_info$dst_t_b)

# df_not_bound <- df_all_info %>% filter(dst_t_b > 0.001)

# ================= 5. Check problematic habitat types ==================
df <- read.csv('data/occ_pts/allinfo_ideam_cgls.csv')
df[is.na(df)] = 0

getInfo <- function(habitat_code){
  thisdf <- df[df[habitat_code]==1,]
  # drop the columns that are not needed
  thisdf <- thisdf[,-which(names(thisdf) %in% c('X', 'occ_ID', habitat_code, 'sorc_fl', 'taxa', 'lynd_nm','nvl_3_n', 'nvl_2_n', 'cgls_value'))]
  
  thisinfo <- as.data.frame(colSums(thisdf %>% select(-species)))
  colnames(thisinfo) <- c('freq')
  rownames(thisinfo) <- seq(1, nrow(thisinfo), 1)
  thisinfo$hab_type <- colnames(thisdf)[-1]
  thisinfo <- merge(thisinfo, habitat_info, by.x = 'hab_type', by.y = 'habitat_code')
  thisinfo$perc <- thisinfo$freq/nrow(thisdf)
  
  return(thisinfo)
}

checkSpecialist <- function(habitat_code){
  thisdf <- df[df[habitat_code]==1,]
  # drop the columns that are not needed
  thisdf <- thisdf[,-which(names(thisdf) %in% c('X', 'occ_ID', habitat_code, 'sorc_fl', 'taxa', 'lynd_nm','nvl_3_n', 'nvl_2_n', 'cgls_value'))]
  
  thisinfo <- getInfo(habitat_code)
  thisdf_0 <- thisdf %>% filter(rowSums(thisdf %>% select(-species))==0)
  
  print(nrow(thisdf_0))
  return(nrow(thisdf_0))
}


habitat_info_sp <- habitat_info
habitat_info_sp$sp_num <- 0
for(i in 1:nrow(habitat_info_sp)){
  code_i <- habitat_info_sp[i,]$habitat_code
  num_i <- checkSpecialist(code_i)
  habitat_info_sp[i,]$sp_num <- num_i
}

# rocky areas
info_6 <- getInfo('hab_6')
ggplot(info_6, aes(x=habitat_name, y=freq)) +
  geom_bar(stat='identity') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# artificial aquatic
info_15 <- getInfo('hab_15')
ggplot(info_15, aes(x=habitat_name, y=freq)) +
  geom_bar(stat='identity') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# desert
info_8 <- getInfo('hab_8')
ggplot(info_8, aes(x=habitat_name, y=freq)) +
  geom_bar(stat='identity') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
