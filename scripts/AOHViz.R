# not used in the end
library(sf)
library(dplyr)
library(here)
library(terra)
library(ggplot2)
library(patchwork)
library(tidyterra)
library(ggspatial)
library(prettymapr)

setwd('/')
setwd('Users/wenxinyang/Desktop/GitHub/colander')
source('scripts/refineBiomodelos/funcs.R')
source('scripts/refineBiomodelos/info.R')
#col_bound <- vect('C:/Users/wenxinyang/Desktop/GitHub/predFHD/data/Colombia/boundary/Colombia_bound.shp') # desktop
col_bound <- vect('/Users/wenxinyang/Desktop/Dissertation/DATA/Colombia_bound/Colombia_bound.shp') # my laptop

# Function to get range maps for different taxa
getRanges <- function(folder_name){
  taxa <- toupper(folder_name)
  if(folder_name == 'Mammals'){
    ranges <- st_read('data/IUCN_range_maps/cleaned_ranges/Cleaned_MAMMALS.shp')
  }
  
  if(folder_name == 'Birds'){
    ranges <- st_read('data/IUCN_range_maps/cleaned_ranges/Cleaned_kindof_birds.gpkg')
  }
  
  if(folder_name %in% c('Reptiles', 'Amphibians')){
    ranges_1 <- st_read(paste0('data/IUCN_range_maps/cleaned_ranges/Cleaned_', taxa, '_PART1.shp'))
    ranges_2 <- st_read(paste0('data/IUCN_range_maps/cleaned_ranges/Cleaned_', taxa, '_PART2.shp'))
    ranges <- list(ranges_1, ranges_2)
  }
  
  return(ranges)
}

# Function to get species range map
getSpeciesRange <- function(spp_name, ranges, folder_name){
  if(folder_name %in% c("Birds", "Mammals")){
    sp_range <- ranges %>% filter(.data$sci_name == spp_name) %>% sf::st_make_valid()
  } else{
    ranges_1 <- ranges[[1]]
    ranges_2 <- ranges[[2]]
    sp_range1 <- ranges_1 %>% filter(.data$sci_name == spp_name) %>% sf::st_make_valid()
    sp_range2 <- ranges_2 %>% filter(.data$sci_name == spp_name) %>% sf::st_make_valid()
    if(nrow(sp_range1)+nrow(sp_range2)>0){
      sp_range <- rbind(sp_range1, sp_range2)
    } else{
      sp_range <- NULL
    }
  }
  
  if(!is.null(sp_range) && nrow(sp_range) > 1){
    sp_range <- sp_range %>% group_by(.data$sci_name) %>% summarise()
  }
  
  # Clip species range by Colombia boundary
  if(!is.null(sp_range) && nrow(sp_range) > 0){
    # Convert col_bound from SpatVector to sf object for clipping
    col_bound_sf <- st_as_sf(col_bound)
    sp_range <- st_intersection(sp_range, col_bound_sf)
  }
  
  return(sp_range)
}

# Define the four species to visualize
species_list <- list(
  list(taxa = 'Birds', name = 'Zenaida auriculata'),
  list(taxa = 'Mammals', name = 'Tapirus pinchaque'),
  list(taxa = 'Reptiles', name = 'Echinosaura orcesi'),
  list(taxa = 'Amphibians', name = 'Andinobates bombetes')
)

# Create plots for each species
plots_list <- list()

for(i in seq_along(species_list)){
  taxa <- species_list[[i]]$taxa
  name <- species_list[[i]]$name
  
  cat("Creating maps for:", taxa, "-", name, "\n")
  
  # Load range maps
  ranges <- getRanges(taxa)
  sp_range <- getSpeciesRange(name, ranges, taxa)

  final_mask <- rast(paste0('data/dsc_aoh/', gsub(" ","_", name ), '.tif'))
  orig_aoh_col <- rast(paste0('data/IUCN_AOH_100m/', taxa, '/', gsub(" ","_", name ), '.tif'))
  orig_aoh_col <- crop(orig_aoh_col, col_bound, mask = TRUE)
  orig_aoh_col <- resample(orig_aoh_col, final_mask, method='near') # make them matching extent
  orig_aoh_col <- ifel(orig_aoh_col == 1, 1, NA)
  
  # Create p1 (National approach)
  if(i != 4 ){
    p1 <- ggplot() +
      annotation_map_tile(type = "cartolight", zoom = 8) +
      geom_spatraster(data = final_mask) +
      scale_fill_stepsn(colors = c("#669c6a", "#669c6a","#dcb202", "#02b3ff", "#02b3ff"),
                        breaks = c(0.5, 1.5, 2.5),
                        #labels = c("1", "2", "3"),
                        na.value = "transparent") +
      # labs(title ='National approach', fill = 'Uncertainty level') +
      theme_minimal() +
      theme(legend.position = "none")
  } else{
    p1 <- ggplot() +
      annotation_map_tile(type = "cartolight", zoom = 8) +
      geom_spatraster(data = final_mask) +
      scale_fill_gradient(low = "#669c6a", high = "#669c6a", na.value = "transparent") +
      # labs(title ='National approach', fill = 'Uncertainty level') +
      theme_minimal() +
      theme(legend.position = "none")
  }

  
  # Add range map to p1 if available
  if(!is.null(sp_range) && nrow(sp_range) > 0){
    p1 <- p1 + geom_sf(data = sp_range, fill = "transparent", color = "black", linewidth = 1)
  }
  
  # Create p2 (Global approach)
  p2 <- ggplot() +
    annotation_map_tile(type = "cartolight", zoom = 8) +
    geom_spatraster(data = orig_aoh_col) +
    scale_fill_gradient(low = "#669c6a", high = "#669c6a", na.value = "transparent") +
    # labs(title = 'Global approach', fill = 'AOH') +
    theme_minimal() +
    theme(legend.position = "none") 
    
  # Add range map to p2 if available
  if(!is.null(sp_range) && nrow(sp_range) > 0){
    p2 <- p2 + geom_sf(data = sp_range, fill = "transparent", color = "black", linewidth = 1)
  }
  
  # Combine p1 and p2
  combined_plot <- p1 + p2 + plot_layout(ncol = 2)
  
  # Add the main legend only to the right side
  #combined_plot <- combined_plot + 
  #  plot_annotation(
  #    theme = theme(plot.title = element_text(hjust = 0.5, size = 14))
  #  ) &
  #  theme(legend.position = "right") &
  #  guides(fill = guide_colorbar(title = "AOH", 
  #                              title.position = "top",
  #                              title.hjust = 0.5,
  #                              barwidth = 1,
  #                              barheight = 10))
  
  plots_list[[i]] <- combined_plot
  
  # Save files to data/viz folder
  speciesname <- gsub(" ", "_", name)
  
  # Create data/viz folder if it doesn't exist
  if(!dir.exists('data/viz')){
    dir.create('data/viz', recursive = TRUE)
  }
  
  # Save range map (1)
  if(!is.null(sp_range) && nrow(sp_range) > 0){
    st_write(sp_range, paste0('data/viz/', speciesname, '_range.shp'), delete_dsn = TRUE)
    cat("Saved range map:", paste0('data/viz/', speciesname, '_range.shp'), "\n")
  }
  
  # Save final_mask (2)
  writeRaster(final_mask, paste0('data/viz/', speciesname, '_national.tif'), overwrite = TRUE)
  cat("Saved final mask:", paste0('data/viz/', speciesname, '_national.tif'), "\n")
  
  # Save orig_aoh_col (3)
  writeRaster(orig_aoh_col, paste0('data/viz/', speciesname, '_lumbierres.tif'), overwrite = TRUE)
  cat("Saved original AOH:", paste0('data/viz/', speciesname, '_lumbierres.tif'), "\n")
  
  # Print the plot
  print(combined_plot)
  rm(ranges)
  gc()
}

# Create a final combined plot showing all four species
if(length(plots_list) == 4){
  final_plot <- plots_list[[1]] / plots_list[[2]] / plots_list[[3]] / plots_list[[4]]
  print(final_plot)
}