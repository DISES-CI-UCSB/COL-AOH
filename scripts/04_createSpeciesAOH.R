# This script creates AOH maps (with uncertainty) for each species
# Author: Wenxin Yang
# Date: July, 2025

# inputs needed
# - species name
# - species range map/IUCN range map or biomodelos in raster format
# - species habitat preference
# - species elevation preference
# - elevation file
# - habitat layers

library(sf)
library(dplyr)
library(here)
library(terra)
library(ggplot2)
library(patchwork)
library(tidyterra)
library(ggspatial)
library(prettymapr)
library(parallel)

setwd('/')
setwd('Users/wenxinyang/Desktop/GitHub/colander')
source('scripts/refineBiomodelos/funcs.R')
source('scripts/refineBiomodelos/info.R')

cleanRangeMaps <- function(folder_name, folder_ranges, if_introduced){
  # -------------------------------------------------------
  # Purpose: filter range maps according to the paper
  # Inputs:
  # folder_name (str): name of the folder of the specific type of range maps, e.g., REPTILES
  # folder_ranges (str): directory where all range maps are stored
  # if_introduced (boolean): if the user wants ranges of introduced species, set the parameter
  # to be 1, else 0
  
  # Output:
  # names_li (list): a list of species scientific names to gather info in the next step
  # -------------------------------------------------------
  tgt_folder <- file.path(folder_ranges, folder_name)
  out_folder <- file.path(folder_ranges, 'cleaned_ranges')
  dir.create(out_folder, showWarnings = FALSE)
  
  range_maps <- list.files(path=tgt_folder, pattern=".shp$")
  n_maps <- length(range_maps)
  
  # add a for loop here
  names_li <- list()
  for(i in 1:length(range_maps)){
    range_map_name <- range_maps[i]
    # read in 
    range_map <- st_read(file.path(tgt_folder, range_map_name))
    # different for birds
    # info <- st_layers('data/IUCN_range_maps/BOTW.gdb')
    # range_map <- st_read(dsn='data/IUCN_range_maps/BOTW.gdb', layer='All_Species')
    #colnames(range_map)
    print(length(unique(range_map$sci_name)))
    
    # filter for extinct or extinct in the wild
    range_map <- range_map %>% dplyr::filter((category!='EX') & (category!='EW'))
    print(length(unique(range_map$sci_name)))
    
    # filter for seasonality
    range_map <- range_map %>% dplyr::filter(seasonal == 1)
    print(length(unique(range_map$sci_name)))
    
    # check origin to get introduced species or not-introduced species
    if(if_introduced==1){
      range_map <- range_map %>% dplyr::filter(origin %in% c(3))
    }else{
      range_map <- range_map %>% dplyr::filter(origin %in% c(1,2,6))
    }
    print(length(unique(range_map$sci_name)))
    
    extant_ranges <- range_map[grep("Extant", range_map$legend), ]
    extant_ranges <- range_map %>% dplyr::filter(presence %in% c(1,2,3))
    print(length(unique(extant_ranges$sci_name)))
    
    # include possibly extinct polygons if a species is critically endangered and does not have any extant or possibly extant polygons
    cr_ranges <- range_map[range_map$category=='CR', ]
    li_sp_extant <- unique(extant_ranges$sci_name)
    li_sp_cr <- unique(cr_ranges$sci_name)
    li_missed_cr <- setdiff(li_sp_cr, li_sp_extant)
    cr_psbext_ranges <- cr_ranges[cr_ranges$sci_name %in% li_missed_cr, ] %>% dplyr::filter(legend =="Possibly Extinct")
    print(length(unique(cr_psbext_ranges$sci_name)))
    
    # combine them
    extant_ranges <- rbind(extant_ranges, cr_psbext_ranges)
    print(length(unique(extant_ranges$sci_name)))
    
    # fix geometries
    cleaned_ranges <- extant_ranges %>% sf::st_make_valid()
    
    # write out cleaned range maps
    if(if_introduced==1){
      write_sf(cleaned_ranges, file.path(out_folder, paste('Cleaned_intr_', range_map_name, sep='')))
    } else{
      write_sf(cleaned_ranges, file.path(out_folder, paste('Cleaned_', range_map_name, sep='')))
    }
    #st_write(obj = extant_ranges, dsn = 'data/IUCN_range_maps/Cleaned_ranges/Cleaned_kindof_birds.gpkg', layer='birds',
    #         delete_dsn = TRUE)
    # got some warnings here, probably due to previous no-data values  
    
    # update species name list
    names_li <- append(names_li, list(unique(range_map$sci_name))[[1]])
  }
  return(names_li)
}

createAOH <- function(spp_name, ranges, folder_name, ifComp, col_srtm, col_bound){
  taxa <- toupper(folder_name)
  
  # print(spp_name)
  # get elev
  tmp_elev <- elev_info %>% filter(name == spp_name) %>% unique()
  elev <- c(tmp_elev$elevation_lower, tmp_elev$elevation_upper)
  
  # get pref
  tmp_pref <- pref_info %>% filter(species == spp_name) %>% unique()
  pref_code <- colnames(tmp_pref)[which(tmp_pref[1, ] == 1)]
  pref <- habitat_info1$habitat_name[habitat_info1$habitat_code %in% pref_code]
  pref <- pref[!pref %in% c('Artificial-Aquatic','Desert','Rocky Areas')]
  
  
  if(length(pref)== 0 | is.null(pref)){
    cat(spp_name, ' does not have matching information, requires manual confirmation')
    cat('\n')
    df <- data.frame(taxa = t, species = spp_name, reason = 'no pref', row.names = FALSE)
    if(!file.exists('data/pref_needs_info.csv')){
      write.table(df, 'data/pref_needs_info.csv', sep = ';', row.names = FALSE)
    }else{
      write.table(df, 'data/pref_needs_info.csv', append=TRUE, col.names = FALSE,
                  row.names = FALSE, sep = ';')
    }
  } else{
    ### ========= 1. create a range subset ========
    tryCatch(
      expr = {
        if(folder_name %in% c("Birds", "Mammals")){
          sp_range <- ranges %>% filter(sci_name == spp_name) %>% sf::st_make_valid()
          rm(ranges)
        } else{
          ranges_1 <- ranges[[1]]
          ranges_2 <- ranges[[2]]
          sp_range1 <- ranges_1 %>% filter(sci_name == spp_name) %>% sf::st_make_valid()
          sp_range2 <- ranges_2 %>% filter(sci_name == spp_name) %>% sf::st_make_valid()
          if(nrow(sp_range1)+nrow(sp_range2)>0){
            sp_range <- rbind(sp_range1, sp_range2)
          } else{
            sp_range <- NULL
          }
          rm(sp_range1, sp_range2, ranges_1, ranges_2)
        }
      },
      error = function(e){
        print(e)
      }
    )
    
    if(is.null(sp_range) | NA %in% st_is_valid(sp_range)){
      if(is.null(sp_range)){
        cat(spp, ' does not have range map \n')
        df <- data.frame(taxa = t, species = spp_name, reason = 'no range map')
      } else{
        cat(spp, ' does not have valid geom')
        df <- data.frame(taxa = t, species = spp_name, reason = 'no valid geom')
      }
      if(!file.exists('data/pref_needs_info.csv')){
        write.table(df, 'data/pref_needs_info.csv', sep = ';', row.names = FALSE)
      }else{
        write.table(df, 'data/pref_needs_info.csv', append=TRUE, col.names = FALSE,
                    row.names = FALSE, sep = ';')
      }
    }
    else{
      if(FALSE %in% st_is_valid(sp_range)){
        sf_use_s2(FALSE)
        sp_range <- st_make_valid(sp_range)
      } 
      if(nrow(sp_range) > 1){
        sp_range <- sp_range %>% group_by(sci_name) %>% summarise()
      }
      
      # plot(sp_range['sci_name'])
      
      ### ======== 2. rasterize ========
      # first read in a template raster
      # then rasterize
      problem <- 0
      sp_range_ras <- rasterize(sp_range, col_srtm)
      tryCatch(expr = {sp_range_ras <- crop(sp_range_ras, sp_range, extend=FALSE)},
               error = function(e){
                 print(e)
                 problem <<- 1
                 df <- data.frame(taxa = t, species = spp_name, reason = as.character(e))
                 write.table(df, 'data/pref_needs_info.csv', append=TRUE, col.names = FALSE,
                             row.names = FALSE, sep=';')
                 })
      # plot(sp_range_ras)
      if(problem==0){
        ### ======== 3. add elev mask ========
        elev_mask <- ifel(col_srtm>=elev[1] & col_srtm <=elev[2], 1, NA)
        
        ### ======= 4. create habitat layer =======
        for(i in 1:length(pref)){
          gc()
          print(pref[i])
          prj_path <- paste0('data/hab_layers/2022/', pref[i], '_layer_prj.tif')
          if(!file.exists(prj_path)){
            prefi_dft <- rast(paste0('data/hab_layers/2022/', pref[i], '_layer.tif'))
            prefi_prj <- project(prefi_dft, "EPSG:4326", method = "near")
            writeRaster(prefi_prj, prj_path, datatype = "INT1U")
            rm(prefi_dft)
            rm(prefi_prj)
            gc()
          }
          prefi <- rast(prj_path)
          
          rsmpl_pref <- resample(prefi, sp_range_ras, method='near')
          
          if(i==1){
            pref_1 = mask(rsmpl_pref, sp_range_ras)
            pref_final = pref_1
            rm(pref_1)
          } else {
            pref_2 <- mask(rsmpl_pref, sp_range_ras)
            pref_final <- ifel(is.na(pref_2) | pref_final<=pref_2, pref_final, pref_2)
            rm(pref_2)
          }
        }
        
        elev_mask <- resample(elev_mask, pref_final)
        final_mask <- mask(pref_final, elev_mask)
        # plot(final_mask)
        # gc()
        
        
        writeRaster(final_mask, paste0('data/dsc_aoh/', gsub(" ","_", spp_name ), '.tif'), overwrite=TRUE)
        
        ### ========= 5. get global aoh map ==========
        if(ifComp == 1){
          orig_aoh <- rast(paste0('data/IUCN_AOH_100m/',taxa, '/', gsub(" ","_", name ), '.tif'))
          orig_aoh_col <- crop(orig_aoh, col_bound, mask = TRUE)
          orig_aoh_col <- resample(orig_aoh_col, final_mask, method='near') # make them matching extent
          orig_aoh_col <- ifel(orig_aoh_col == 1, 1, NA)
          # plot(orig_aoh_col)
          gc()
          
          p1 <- ggplot() +
            annotation_map_tile(type = "osm", zoom = 8) +  # OpenStreetMap topographic base
            geom_spatraster(data = final_mask) +
            scale_fill_gradient(low = "blue", high = "gray", na.value = "transparent") +
            labs(title = 'National approach', fill = 'AOH') +
            theme_minimal()
          
          p2 <- ggplot() +
            annotation_map_tile(type = "osm", zoom = 8) +  # OpenStreetMap topographic base
            geom_spatraster(data = orig_aoh_col) +
            scale_fill_gradient(low = "blue", high = "blue", na.value = "transparent") +
            labs(title = 'Lumbierres et al.', fill = 'AOH') +
            theme_minimal()
          
          
          p1 + p2 + plot_layout(ncol=2)
        }
      }
      
    }
    
  }

}

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


# ============= Actual running ==============

# crop DEM by Colombia boundary
col_bound <- vect('C:/Users/wenxinyang/Desktop/GitHub/predFHD/data/Colombia/boundary/Colombia_bound.shp')

#srtm <- rast('data/srtm.tif')
#col_srtm <- crop(srtm, col_bound, mask = TRUE)
#writeRaster(col_srtm, 'data/col_srtm.tif')

col_srtm <- rast('data/col_srtm.tif')
folder_ranges <- here('data/IUCN_range_maps')
# if need to clean anything, use the above function

# get species list
star_species <- read.csv('data/CI_STAR_T_Species_List.csv')

# get species elevation info
elev_info <- read.csv('data/animals_preference.csv', sep=';') %>% select(name, elevation_upper, elevation_lower)
# set NA elevations as 0 - 9000 m
elev_info <- elev_info %>%
  mutate(elevation_lower = ifelse(is.na(elevation_lower), 0, elevation_lower)) %>%
  mutate(elevation_upper = ifelse(is.na(elevation_upper), 9000, elevation_upper)) %>% unique()

# get species preference info
pts <- read.csv('data/occ_pts/allinfo_ideam_cgls_coords_2022.csv')
cols <- colnames(pts)
cols_to_keep <- cols[grepl("^hab", cols)]
cols_to_keep <- c('species', cols_to_keep)
pref_info <- pts %>% select(all_of(cols_to_keep)) %>% unique()
# merge several artificial habitat types
pref_info <- pref_info %>% mutate(
  hab_14.12 = ifelse(hab_14.1+hab_14.2 >0, 1, 0),
  hab_14.36 = ifelse(hab_14.3+hab_14.6 >0, 1, 0),
  hab_14.45 = ifelse(hab_14.4+hab_14.5 >0, 1, 0)
)

# remove the original ones
pref_info <- pref_info %>% select(-all_of(c('hab_14.1', 'hab_14.2', 'hab_14.3', 'hab_14.4', 'hab_14.5', 'hab_14.6')))

all_taxa <- c('Birds', 'Mammals', 'Amphibians', 'Reptiles')
for(t in all_taxa){
  print(t)
  if(file.exists('data/pref_needs_info.csv')){
    prob_sp <- read.csv('data/pref_needs_info.csv', sep = ';')
    sp_li <- unique(prob_sp$species)
  } else {
    sp_li <- character(0)  # empty vector if file doesn't exist
  }
  # prob_sp <- read.csv('data/pref_needs_info.csv')
  sp_t <- star_species %>% filter(taxonomic_group == t) %>% 
    select('species') %>% unique()
  sp_t <- sp_t$species
  
  if(t == 'Amphibians'){t_es = 'anfibios'}
  if(t == 'Reptiles'){t_es = 'squamata'}
  if(t == 'Birds'){t_es = 'aves'}
  if(t == 'Mammals'){t_es = 'mamiferos'}
  sp_t_occ <- read.csv(paste0('data/occ_pts/', t_es, '_occ_pts.csv'))
  if('species' %in% colnames(sp_t_occ)){
    sp_t_occ <- unique(sp_t_occ$species)
  } else{
    sp_t_occ <- unique(sp_t_occ$scientificName)
  }
  
  sp_both <- intersect(sp_t, sp_t_occ)
  
  # print(length(sp_both))
  
  #ranges <- getRanges(t)
  
  # progress file to track what's happening
  progress_file <- paste0('data/tmp/progress_', t, '.txt')
  dir.create('data/tmp', showWarnings = FALSE, recursive = TRUE)
  
  # Function to process a single species (for parallel processing)
  process_species <- function(spp) {
    tryCatch({
      gc()
      # Write to progress file (this will definitely work)
      write(paste(Sys.time(), spp, "STARTED", sep = " | "), 
            file = progress_file, append = TRUE)
      
      if(!file.exists(paste0('data/dsc_aoh/', gsub(" ","_", spp ), '.tif'))){
        if(!spp %in% sp_li){
          # Write progress: loading ranges
          write(paste(Sys.time(), spp, "LOADING_RANGES", sep = " | "), 
                file = progress_file, append = TRUE)
          
          # Load ranges on this worker instead of exporting the large object
          worker_ranges <- getRanges(t)
          
          # Load spatial objects for this worker
          col_bound <- vect('C:/Users/wenxinyang/Desktop/GitHub/predFHD/data/Colombia/boundary/Colombia_bound.shp')
          col_srtm <- rast('data/col_srtm.tif')
          
          # Write progress: ranges loaded, starting createAOH
          write(paste(Sys.time(), spp, "RANGES_LOADED_STARTING_AOH", sep = " | "), 
                file = progress_file, append = TRUE)
          
          createAOH(spp, worker_ranges, t, 0, col_srtm, col_bound)
          
          # Write progress: createAOH completed
          write(paste(Sys.time(), spp, "AOH_COMPLETED", sep = " | "), 
                file = progress_file, append = TRUE)
          
          rm(worker_ranges)
          gc()
        } else {
          write(paste(Sys.time(), spp, "SKIPPED_PROBLEM_LIST", sep = " | "), 
                file = progress_file, append = TRUE)
        }
      } else {
        write(paste(Sys.time(), spp, "SKIPPED_FILE_EXISTS", sep = " | "), 
              file = progress_file, append = TRUE)
      }
      
      write(paste(Sys.time(), spp, "FINISHED", sep = " | "), 
            file = progress_file, append = TRUE)
      return(list(species = spp, status = "success"))
    }, error = function(e) {
      write(paste(Sys.time(), spp, "ERROR", e$message, sep = " | "), 
            file = progress_file, append = TRUE)
      return(list(species = spp, status = "error", error = e$message))
    })
  }
  
  # Set up parallel processing
  # Use fewer cores to avoid memory issues (each worker loads full datasets)
  max_cores <- min(3, detectCores() - 2)  # Use max 6 cores or leave 2 cores free, whichever is smaller
  n_cores <- min(max_cores, length(sp_t_occ))  # Don't use more cores than species
  cat("Processing", length(sp_t_occ), "species using", n_cores, "cores\n")
  
  if (n_cores > 1 && length(sp_t_occ) > 1) {
    # Create cluster with output enabled
    cl <- makeCluster(n_cores, outfile = "")
    
    # Export necessary variables and functions to cluster
    # Note: ranges, col_srtm, and col_bound are NOT exported - they contain external pointers
    # Each worker will load them to avoid serialization issues
    clusterExport(cl, c("createAOH", "getRanges", "t", "sp_li", 
                        "elev_info", "pref_info", "habitat_info1", 
                        "progress_file"),
                  envir = environment())
    
    # Load required packages and source files on cluster
    # Also load spatial objects that can't be serialized
    clusterEvalQ(cl, {
      library(sf)
      library(dplyr)
      library(here)
      library(terra)
      library(ggplot2)
      library(patchwork)
      library(tidyterra)
      library(ggspatial)
      library(prettymapr)
      
      # Source required files
      setwd('/')
      setwd('Users/wenxinyang/Desktop/GitHub/colander')
      source('scripts/refineBiomodelos/funcs.R')
      source('scripts/refineBiomodelos/info.R')
      
      # Note: col_bound and col_srtm are loaded fresh for each species in process_species
      # to avoid external pointer issues with terra objects in parallel workers
    })
    
    cat("Starting parallel processing...\n")
    cat("Progress is being written to:", progress_file, "\n")
    cat("You can check this file to see what workers are doing:\n")
    cat("  tail -f", progress_file, "\n")
    cat("========================================\n")
    flush.console()
    
    # Process species in parallel with load balancing for better progress visibility
    results <- parLapplyLB(cl, sp_t_occ, process_species)
    
    # Print final progress summary
    if(file.exists(progress_file)) {
      cat("\n=== Progress Summary ===\n")
      progress_data <- readLines(progress_file)
      started <- sum(grepl("STARTED", progress_data))
      loading <- sum(grepl("LOADING_RANGES", progress_data))
      ranges_loaded <- sum(grepl("RANGES_LOADED", progress_data))
      aoh_completed <- sum(grepl("AOH_COMPLETED", progress_data))
      errors <- sum(grepl("ERROR", progress_data))
      cat("Started:", started, "\n")
      cat("Loading ranges:", loading, "\n")
      cat("Ranges loaded:", ranges_loaded, "\n")
      cat("AOH completed:", aoh_completed, "\n")
      cat("Errors:", errors, "\n")
      cat("Check", progress_file, "for details\n")
    }
    
    # Stop cluster
    stopCluster(cl)
    
    # Print summary
    success_count <- sum(sapply(results, function(x) x$status == "success"))
    error_count <- sum(sapply(results, function(x) x$status == "error"))
    cat("\n=== Processing Summary ===\n")
    cat("Successfully processed:", success_count, "species\n")
    cat("Errors:", error_count, "species\n")
    
  } else {
    # Fall back to sequential processing if only 1 core or 1 species
    cat("Using sequential processing (1 core or 1 species)\n")
    results <- lapply(sp_t_occ, process_species)
  }
}

gc()

