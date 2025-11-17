# This script takes a translation matrix and creates subsequent habitat layers
# from a land cover layer
# Author: Wenxin Yang
# Date: July, 2025

# ====== load libraries =======
packages <- c("sf", "terra", "here", "rgrass", "dplyr", "data.table", 
              "rredlist", "tidyr", "Matrix", "tidyverse", "pROC",
              "kableExtra", "ggplot2", "reshape2")
lapply(packages, library, character.only = TRUE)

setwd('/')
setwd('Users/wenxinyang/Desktop/GitHub/colander')
source('scripts/refineBiomodelos/funcs.R')
source('scripts/refineBiomodelos/info.R')

# ======= read in and preprocess translation matrix =======
df.mat <- read.csv('results/aoh_results_randomCI_final_7gen_2022/btst_ideam_randomCI_gen7__noDesRockAA_1104_2025_pos.csv')
df.mat$X <- NULL

habitat_cols <- c("Artificial-arable and pasture", "Artificial-degraded forest and plantation", "Artificial-urban areas and rural gardens", "Grassland", "Savanna", "Shrubland", "Wetlands (inland)", "Forest")

mat_cols <- c("land_cover", habitat_cols, "auc", "n_samples")
colnames(df.mat) <- mat_cols

## ====== convert df.mat to long =======
df.mat.long <- df.mat %>%
  tidyr::pivot_longer(cols = all_of(habitat_cols), names_to = "habitat",
                      values_to = "count") %>% select(-c("auc", "n_samples"))

# keep pairs
df.mat.sel <- df.mat.long %>%
  filter(count > 300) %>%
  mutate(uncertainty_level = ifelse(count > 900, 1, 
                              ifelse(count > 600, 2, 
                                     ifelse(count > 300, 3, 0)))) %>%
  select(-c('count'))

## ====== convert land cover name to code =======
df.mat.sel <- merge(df.mat.sel, ideam_lc_info, by.x='land_cover', 
                    by.y='ideam_lc_name', all.x=TRUE)
df.mat.sel <- df.mat.sel %>% mutate(
  lc_code = sapply(strsplit(ideam_lc_code, "_"), function(x) x[2])
) %>% select(!ideam_lc_code)


# ============ pre-proecess land cover layer ========
lc <- vect('data/Corine_hab_COL/Cobertura_tierra_100K_periodo_2022_limite_administrativo/ECOSISTEMAS_18062025/ECOSISTEMAS_18062025.gpkg')

lc_prj <- project(lc, "epsg:3116")

lc_prj$nivel_2 <- as.integer(gsub("[^0-9]", "", lc_prj$nivel_2))

# unique(lc_prj$nivel_2) %in% df.mat.sel$lc_code

# rasterize to EPSG 3116 and 100 m resolution using terra package
# create a template raster with the desired CRS and resolution
lc_bbox <- st_bbox(lc_prj)
template_rast <- rast(crs = "EPSG:3116", 
                      xmin = lc_bbox[["xmin"]], 
                      xmax = lc_bbox[["xmax"]], 
                      ymin = lc_bbox[["ymin"]], 
                      ymax = lc_bbox[["ymax"]], 
                      resolution = 100)

# rasterize the land cover vector data to the template
lc_rast <- terra::rasterize(lc_prj, template_rast, field = 'nivel_2')

writeRaster(lc_rast, 'data/Corine_hab_COL/ideam_2022_level2_100m.tif', overwrite = TRUE)
# plot(lc_rast)

# lc_rast <- rast('data/Corine_hab_COL/ideam_2022_level2_100m.tif')
# ====== convert land cover to habitat layers with uncertainty =======
for(hab in habitat_cols){
  print(hab)
  sel_hab <- df.mat.sel %>% filter(habitat %in% hab)
  ## ====== if base layer is vector data =======
  #lc_hab_shp <- merge(lc, sel_hab, by.x='nivel_2', by.y = 'lc_code')
  #colnames(lc_hab_shp)
  
  # dissolve lc_hab_shp by uncertainty level
  #lc_hab_shp_diss <- lc_hab_shp %>% group_by(uncertainty_level) %>%
  #  summarise(geometry = st_union(geometry))
  
  ## ====== if base layer is raster data ======
  # convert habitat land cover info to a matrix
  sel_hab1 <- sel_hab %>% select(lc_code, uncertainty_level) %>% mutate(
    from = as.integer(lc_code)-0.5,
    to = as.integer(lc_code)+0.5,
    values = as.integer(uncertainty_level)
  ) %>% select(!c(lc_code, uncertainty_level))
  sel_mt <- as.matrix(sel_hab1)
  sel_lc_codes <- unique(as.integer(sel_hab$lc_code))
  
  hab_msk <- ifel(lc_rast %in% sel_lc_codes, lc_rast, NA)
  # plot(hab_msk)
  hab_rast <- mask(lc_rast, hab_msk)
  
  recl <- classify(hab_rast, sel_mt, right = TRUE)
  # plot(recl)
  
  if(!dir.exists('data/hab_layers')){
    dir.create('data/hab_layers')
  }
  
  if(!dir.exists('data/hab_layers/2022')){
    dir.create('data/hab_layers/2022')
  }
  
  writeRaster(recl, file.path('data/hab_layers/2022/', paste0(hab, '_layer.tif')), overwrite = TRUE)
  
  gc()
}


test <- rast('data/hab_layers/2022/Savanna_layer.tif')
plot(test)
