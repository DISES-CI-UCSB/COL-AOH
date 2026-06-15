# Prepare species occurrence data & match /w land cover data
# Author: Wenxin Yang
# Revision Date: June, 2026

# ================== Prep =====================
# Load libraries
packages <- c("sf", "terra", "here", "rgrass", "dplyr", "data.table", 
              "rredlist", "tidyr", "Matrix", "tidyverse", "pROC", "stringr")
lapply(packages, library, character.only = TRUE)

setwd('/')
setwd('Users/wenxinyang/Desktop/GitHub/colander')
source('scripts/refineBiomodelos/funcs.R')
source('scripts/refineBiomodelos/info.R')

col_to_keep <- c('scientificName', 'decimalLatitude', 'decimalLongitude', 
                 'year', 'eventDate', 
                 #'createdDate', 'reportedDate', 'yyyy','eventTime',
                 'coordinateUncertaintyInMeters', 'taxonomicStatus', 'duplicated',
                 'alt', 'extremo', 'spatialDuplicated',
                 'source_file')

analysis_cols <- c('scientificName', 'decimalLatitude', 'decimalLongitude',
                   'year', 'taxa', 'eventYear')

# ================= 0. Read in and clean occurrence points ==================
## -------------- amphibians ---------
dfanf <- read.csv('data/raw_occ_pts/anfibios_raw_pts.csv') %>% select(all_of(col_to_keep))
length(unique(dfanf$scientificName))
dfanf$taxa <- 'anfibios'
summary(dfanf$year) # no NA values
dfanf <- dfanf %>% mutate(
  eventYear = year
) 
nrow(dfanf %>% filter(scientificName==""))
dfanfAll <- dfanf %>% select(all_of(analysis_cols))

## ------------- birds ---------------
dfavs <- read.csv('data/raw_occ_pts/aves_raw_pts.csv') %>% 
  # if a row does not have scientific name extract it from source file name
  select(all_of(col_to_keep)) %>% mutate( 
    scientificName = ifelse(
      scientificName=="", 
      source_file %>% str_remove("\\.csv$") %>% str_replace_all("_", " "), 
      scientificName)
)
nrow(dfavs %>% filter(scientificName==""))
length(unique(dfavs$scientificName))
dfavs$taxa <- 'aves'
summary(dfavs$year) # has 0 and NA values

dfavs <- dfavs %>%
  mutate(
    year = case_when(
      year > 0 ~ year,
      TRUE ~ NA
    ),
    eventDate = case_when(!eventDate %in% c("NA--NA", "") ~ eventDate, 
                          TRUE ~ '0'),
    eventYear = case_when(
      str_detect(eventDate, "^\\d{4}-") ~ as.character(str_extract(eventDate, "\\d{4}")),
      str_detect(eventDate, "^[^/]*/[^/]*/") ~ as.character(str_match(eventDate, "^[^/]*/[^/]*/(\\d{4})")[,2]),
      str_detect(eventDate, "^\\d{4}/*") ~ as.character(str_extract(eventDate, "^\\d{4}[^/]*")),
      TRUE ~ as.character(eventDate)
    ),
    eventYear = as.numeric(eventYear),
    eventYear = case_when(eventYear>1700 & eventYear<2025 ~ eventYear, 
                          TRUE ~ NA)
    )

dfavsAll <- dfavs %>% select(all_of(analysis_cols))

## ------------- mammals ---------------
dfmam <- fread('data/raw_occ_pts/mamiferos_raw_pts.csv') %>% 
  select(all_of(col_to_keep),'createdDate','reportedDate') %>%
  mutate(reportedDate = as.character(reportedDate),
         createdDate = as.character(createdDate),
         scientificName = ifelse(
           scientificName=="", 
           source_file %>% str_remove("\\.csv$") %>% str_replace_all("_", " "), 
           scientificName))
length(unique(dfmam$scientificName))
dfmam$taxa <- 'mamiferos'
summary(dfmam$year) # has 0 and NA values
nrow(dfmam %>% filter(scientificName==""))

dfmam <- dfmam %>%
  mutate(
    eventDate = case_when(!eventDate %in% c("NA--NA", "") | (reportedDate!="") & (createdDate!="") ~ eventDate, 
                          TRUE ~ '0'),
    eventYear = case_when(
      str_detect(eventDate, "^\\d{4}-") ~ as.character(str_extract(eventDate, "\\d{4}")),
      str_detect(eventDate, "^[^/]*/[^/]*/") ~ as.character(str_match(eventDate, "^[^/]*/[^/]*/(\\d{4})")[,2]),
      str_detect(eventDate, "^\\d{4}/*") ~ as.character(str_extract(eventDate, "^\\d{4}[^/]*")),
      eventDate=='0' & reportedDate !="" ~ as.character(str_extract(reportedDate, "\\d{4}")),
      eventDate=='0' & createdDate !="" ~ as.character(str_extract(createdDate, "\\d{4}")),
      TRUE ~ as.character(eventDate)
  ),
  eventYear = as.numeric(eventYear),
  eventYear = case_when(eventYear>1700 & eventYear<2025 ~ eventYear, 
                        TRUE ~ NA)
) %>% select(
  -c(reportedDate, createdDate)
)
summary(dfmam$eventYear)

dfmamAll <- dfmam %>% select(all_of(analysis_cols))

## ------------------- reptiles -------------------
dfrep <- read.csv('data/raw_occ_pts/squamata_raw_pts.csv') %>% 
  select(all_of(col_to_keep)) %>% mutate(
    scientificName = ifelse(
      scientificName=="", 
      source_file %>% str_remove("\\.csv$") %>% str_replace_all("_", " "), 
      scientificName)
  )
length(unique(dfrep$scientificName))
dfrep$taxa <- 'squamata'
summary(dfrep$year)
nrow(dfrep %>% filter(scientificName==""))


dfrep <- dfrep %>%
  mutate(
    eventDate = case_when(!eventDate %in% c("NA--NA", "") ~ eventDate, 
                          TRUE ~ '0'),
    eventYear = case_when(
      str_detect(eventDate, "^\\d{4}-") ~ as.character(str_extract(eventDate, "\\d{4}")),
      str_detect(eventDate, "^[^/]*/[^/]*/") ~ as.character(str_match(eventDate, "^[^/]*/[^/]*/(\\d{4})")[,2]),
      str_detect(eventDate, "^\\d{4}/*") ~ as.character(str_extract(eventDate, "^\\d{4}[^/]*")),
      TRUE ~ as.character(eventDate)
    ),
    eventYear = as.numeric(eventYear),
    eventYear = case_when(eventYear>1700 & eventYear<2025 ~ eventYear, 
                          TRUE ~ NA)
  )

dfrepAll <- dfrep %>% select(all_of(analysis_cols))

# ================= 1. Merge occurrence point files ==================
full_list <- rbind(dfanfAll, dfavsAll, dfmamAll, dfrepAll) %>% mutate(
  finalYear = ifelse(!is.na(year), year, eventYear)
)
nrow(full_list) == nrow(dfanfAll) + nrow(dfavsAll) + nrow(dfmamAll) + nrow(dfrepAll)

full_list_time <- full_list %>% filter(!is.na(finalYear) & finalYear >0)

nrow(full_list) - nrow(full_list_time)

names_li <- unique(full_list_time$species)

## ---------------- temporal filtering ---------------------
year_rec = 2012
recent_list = full_list_time %>% filter(finalYear >= year_rec)
#recent_list <- full_list
recent_list$id <- 1:nrow(recent_list)
## ---------------- spatial thinning --------------------
# need a fishnet of Colombia at 100m resolution to match with the AOH maps
rast_tmplt <- rast('data/Corine_hab_COL/ideam_2022_level2_100m.tif')
rast_tmplt <- project(rast_tmplt, 'EPSG:4326')
cell_ids <- cellFromXY(rast_tmplt, recent_list[, c('decimalLongitude', 'decimalLatitude')])
length(cell_ids) == nrow(recent_list)
sum(is.na(cell_ids))

recent_list_thinned <- recent_list %>%
  mutate(cell_id = cell_ids) %>%
  filter(!is.na(cell_id)) %>%
  group_by(scientificName, cell_id) %>%
  slice(1) %>%
  ungroup() %>%
  select(-cell_id)

names_before <- unique(recent_list$scientificName)
names_thinned <- unique(recent_list_thinned$scientificName)
setdiff(names_before, names_thinned)

  
nrow(recent_list_thinned)/nrow(recent_list)
length(unique(recent_list_thinned$scientificName))
length(unique(recent_list$scientificName))

## --------------- other check points -----------------
# check for spatial coverage criteria
# count species with < 10 pts
occ_count <- as.data.frame(table(recent_list_thinned$scientificName))
nrow(occ_count[occ_count$Freq < 10,])/nrow(occ_count)
# ~ 40% sepcies have < 10 points, we are not removing them
perc_occ <- 0.02*nrow(full_list)
nrow(occ_count[occ_count$Freq > perc_occ,])
# 0 species has more than 2% occ points

# ================= 2. Join points with species preference file ================ 
pref_path <- file.path('data/occ_pts', 'animals_preference.csv')
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


lc_1 <- lc %>% select(leyenda_num, nivel_1_num, nivel_3_num, nivel_2_num)
# , geometry
# extract land cover types for pt values

## ------------ (1) convert full_list into a point sf object -------------
# first clean and validate coordinates
final_list <- recent_list_thinned %>%
  # remove rows with NA coordinates
  filter(!is.na(decimalLongitude) & !is.na(decimalLatitude)) %>%
  # filter to only valid coordinate ranges
  filter(decimalLongitude >= -180 & decimalLongitude <= 180) %>%
  filter(decimalLatitude >= -90 & decimalLatitude <= 90)

# convert to sf with error checking
all_pts <- final_list %>% 
    st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), 
             crs = 4326) %>%
  st_transform(crs = st_crs(lc_1))
# export all_pts
# st_write(all_pts, 'data/occ_pts/all_pts.shp')

## --------------- (2) spatial join all_pts with lc_1 ---------------
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
st_write(all_pts_lc_col, 'data/occ_pts/all_pts_lc_col_0605_2012.shp', append = FALSE)
# ================= 4. Merge and format data ==================
#all_pts_lc_col <- st_read('data/occ_pts/all_pts_lc_col_0605_2012.shp')

# merge all_pts_lc_col with df_pref_sp
spp_field <- ifelse('scntfcN' %in% colnames(all_pts_lc_col), 'scntfcN', 'scientificName')
all_pts_lc_col_pref <- merge(all_pts_lc_col, df_pref_sp, by.x = spp_field, by.y = 'name', all.x = TRUE)
# remove geogmetry of all_pts_lc_col_pref
df_all_info <- all_pts_lc_col_pref
# Add longitude and latitude columns from geometry
coords <- st_coordinates(df_all_info$geometry)
df_all_info$longitude <- coords[,1]
df_all_info$latitude <- coords[,2]

df_all_info$geometry <- NULL
write.csv(df_all_info, 'data/occ_pts/allinfo_ideam_coords_2012_0605.csv')

gc()

# ================= 5. Check problematic habitat types ==================
df <- read.csv('data/occ_pts/allinfo_ideam_coords_2012_0605.csv') %>% drop_na()
df <- df %>% mutate(
  nvl_2_n = gsub('.', '', nvl_2_n, fixed = TRUE)
)

unique(df$nvl_2_n)

cols_remove <- c('X', 'occ_ID', 'sorc_fl', 'taxa', 'leyenda_num','nivel_3_num', 'nivel_2_num', 'nivel_1_num', 'finalYear', 'longitude', 'latitude', 'lynd_nm','nvl_1_n', 'nvl_3_n', 'nvl_2_n', 'finalYr', 'eventYr', 'year', 'id', 'eventYear')
## --------------- get how many species per hab pref ---------------
head(df)
df_spp_hab <- df %>% select(-any_of(cols_remove)) %>% unique()
info_spp_hab <- as.data.frame(colSums(df_spp_hab %>% select(-any_of(c('scntfcN', 'scientificName')))))
colnames(info_spp_hab) <- 'N_speices'
info_spp_hab$hab_code <- rownames(info_spp_hab)
rownames(info_spp_hab) <- 1:nrow(info_spp_hab)
info_spp_hab <- merge(info_spp_hab, habitat_info, by.x='hab_code', by.y='habitat_code')
## --------------- get how many data per hab pref ---------------
df_N_hab <- as.data.frame(colSums(df %>% select(starts_with('hab_')) %>% drop_na()))
colnames(df_N_hab) <- 'N_occ'
df_N_hab$hab_code <- rownames(df_N_hab)
rownames(df_N_hab) <- 1:nrow(df_N_hab)
df_N_hab <- merge(df_N_hab, habitat_info, by.x='hab_code', by.y='habitat_code')

## --------------- get pt breakdown by habitat and landcover ----------------
# 1. Aggregate and sum the habitat columns by your land cover groups
heatmap_data <- df %>%
  group_by(nvl_2_n) %>%
  summarise(across(starts_with("hab_"), ~ sum(., na.rm = TRUE))) %>%
  
  # 2. Pivot the data to a long format ready for ggplot
  pivot_longer(
    cols = starts_with("hab_"), 
    names_to = "Habitat_Variable", 
    values_to = "Total_Sum"
  ) %>% mutate(
    lc_code = paste0('lc_', as.character(nvl_2_n))
  ) %>% group_by(nvl_2_n) %>%
  mutate(
    Row_Total = sum(Total_Sum),
    Percentage = ifelse(Row_Total==0, 0, (Total_Sum/Row_Total)*100),
    Percentage = round(Percentage, 2)
  ) %>% ungroup()

heatmap_data <- merge(heatmap_data, habitat_info1, by.x='Habitat_Variable', by.y='habitat_code')
heatmap_data <- merge(heatmap_data, ideam_lc_info, by.x='lc_code', by.y='ideam_lc_code')


# 3. Create the ggplot heatmap
ggplot(heatmap_data, aes(x = habitat_name, y = as.factor(ideam_lc_name), fill = Total_Sum)) +
  geom_tile(color = "white", lwd = 0.5, linetype = 1) + 
  
  # Add the text labels to each cell
  geom_text(
    aes(label = Total_Sum), 
    color = "white",       # Sets text color (change to "black" if using a light palette)
    fontface = "bold", 
    size = 3.5             # Adjust size to fit your grid perfectly
  ) +
  
  scale_fill_viridis_c(name = "Total Count", option = "mako") + 
  theme_minimal() +
  labs(
    title = "Habitat Column Density by Land Cover Class (nvl_2_n)",
    x = "Habitat Variables",
    y = "Land Cover Class (nvl_2_n)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14)
  )


# 3. Create the ggplot heatmap
ggplot(heatmap_data, aes(x = habitat_name, y = as.factor(ideam_lc_name), fill = round(Percentage))) +
  geom_tile(color = "white", lwd = 0.5, linetype = 1) + 
  
  # Add the text labels to each cell
  geom_text(
    aes(label = Percentage), 
    color = "white",       # Sets text color (change to "black" if using a light palette)
    fontface = "bold", 
    size = 3.5             # Adjust size to fit your grid perfectly
  ) +
  
  scale_fill_viridis_c(name = "Percentage", option = "mako") + 
  theme_minimal() +
  labs(
    title = "Habitat Column Density by Land Cover Class (nvl_2_n)",
    x = "Habitat",
    y = "Land Cover Class (nvl_2_n)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14)
  )


habitats_with_all_low_df <- heatmap_data %>% 
  group_by(Habitat_Variable) %>% 
  summarise(maxperc = max(Percentage)) %>% filter(maxperc<1)

habitats_with_all_low_df <- unique(habitats_with_all_low_df$Habitat_Variable)

## ------------------ other ---------------------

getInfo <- function(habitat_code){
  thisdf <- df[df[habitat_code]==1,]
  # drop the columns that are not needed
  thisdf <- thisdf[,-which(names(thisdf) %in% c(cols_remove, habitat_code))]
  
  thisinfo <- as.data.frame(colSums(thisdf %>% select(-any_of(c('scntfcN','scientificName')))))
  colnames(thisinfo) <- c('freq')
  rownames(thisinfo) <- 1:nrow(thisinfo)
  thisinfo$hab_type <- colnames(thisdf)[-1]
  thisinfo <- merge(thisinfo, habitat_info, by.x = 'hab_type', by.y = 'habitat_code')
  thisinfo$perc <- thisinfo$freq/nrow(thisdf)
  
  return(thisinfo)
}

checkSpecialist <- function(habitat_code){
  thisdf <- df[df[habitat_code]==1,]
  # drop the columns that are not needed
  thisdf <- thisdf[,-which(names(thisdf) %in% c(cols_remove, habitat_code))]
  
  thisinfo <- getInfo(habitat_code)
  thisdf_0 <- thisdf %>% filter(rowSums(thisdf %>% select(-any_of(c('scntfcN','scientificName'))))==0)
  
  print(nrow(thisdf_0))
  return(nrow(thisdf_0))
}


habitat_info_sp <- habitat_info
habitat_info_sp$sp_num <- 0
for(i in 1:nrow(habitat_info_sp)){
  print(habitat_info_sp[i,]$habitat_name)
  code_i <- habitat_info_sp[i,]$habitat_code
  num_i <- checkSpecialist(code_i)
  habitat_info_sp[i,]$sp_num <- num_i
}



write.csv(habitat_info_sp, 'results/habitat_info_specialist_2012onwards.csv')
