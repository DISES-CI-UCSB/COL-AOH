# This script takes a translation matrix and creates subsequent habitat layers
# with continuous stability from a land cover layer
# Author: Wenxin Yang
# Date: July, 2025
# Modified: June, 2026

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
df.mat <- read.csv('results/glm_btst_2012_keephab/gen7_glm_btst_2012_keephab/btst_ideam_randomCI_gen7__0615_2026_pos.csv')
df.mat$X <- NULL

if("Artificial.Aquatic" %in% colnames(df.mat)){
  df.mat <- df.mat %>% mutate(Artificial.aquatic = Artificial.Aquatic) %>% select(-Artificial.Aquatic)
}
if("Rocky.Areas" %in% colnames(df.mat)){
  df.mat <- df.mat %>% mutate(Rocky.areas = Rocky.Areas) %>% select(-Rocky.Areas)
}

colnames(df.mat) <- gsub("\\.", " ", colnames(df.mat))
df.mat <- df.mat %>% rename(
  `Wetlands inland` = `Wetlands  inland `
)

habitat_cols <- c("Artificial aquatic","Artificial arable and pasture", "Artificial degraded forest and plantation", "Artificial urban areas and rural gardens", "Desert", "Forest", "Grassland", "Rocky areas", "Savanna", "Shrubland", "Wetlands inland")

mat_cols <- c("land_cover", habitat_cols, "auc", "n_samples")
df.mat <- df.mat %>% select(all_of(mat_cols))

## ====== convert df.mat to long =======
df.mat.long <- df.mat %>%
  tidyr::pivot_longer(cols = all_of(habitat_cols), names_to = "habitat",
                      values_to = "count") %>% select(-c("auc", "n_samples"))

## ====== convert land cover name to code =======
df.mat.long <- merge(df.mat.long, ideam_lc_info, by.x='land_cover', 
                    by.y='ideam_lc_name', all.x=TRUE)
df.mat.long <- df.mat.long %>% mutate(
  lc_code = sapply(strsplit(ideam_lc_code, "_"), function(x) x[2])) %>% 
  select(-ideam_lc_code)


# ============ pre-proecess land cover layer ========
# SKIP THIS IF ALREADY DONE
if(!file.exists('data/Corine_hab_COL/ideam_2022_level2_100m.tif')){
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
} else{
  lc_rast <- rast('data/Corine_hab_COL/ideam_2022_level2_100m.tif')
}

# ====== convert land cover to habitat layers with uncertainty =======
for(hab in habitat_cols){
  print(hab)
  sel_hab <- df.mat.long %>% filter(habitat %in% hab)

  # convert habitat land cover info to a matrix
  sel_hab1 <- sel_hab %>% select(lc_code, count) %>% mutate(
    from = as.integer(lc_code)-0.5,
    to = as.integer(lc_code)+0.5,
    values = as.integer(count)
  ) %>% select(!c(lc_code, count))
  sel_mt <- as.matrix(sel_hab1)
  sel_lc_codes <- unique(as.integer(sel_hab$lc_code))
  
  
  count_mt <- as.matrix(sel_hab %>% 
                          select(lc_code, count) %>% 
                          mutate(
                            lc_code = as.numeric(lc_code),
                            count = as.numeric(count)
                            )
  )
  
  recl <- classify(lc_rast, count_mt)
  # plot(recl)
  
  if(!dir.exists('data/hab_layers')){
    dir.create('data/hab_layers')
  }
  
  if(!dir.exists('data/hab_layers/2012_btst_7gen_keephab')){
    dir.create('data/hab_layers/2012_btst_7gen_keephab')
  }
  
  writeRaster(recl, file.path('data/hab_layers/2012_btst_7gen_keephab', paste0(hab, '_layer.tif')), overwrite = TRUE)
  
  gc()
}


#test <- rast('data/hab_layers/2012_btst_7gen_keephab/Savanna_layer.tif')
#plot(test)
