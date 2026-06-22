# Test and select threshold to convert continuous maps into binary maps
# based on validation
# Author: Wenxin Yang
# Date: September, 2025
# Modified: June, 2026 (Added parallel processing)

# Model prevalence: AOH size / range size (km^2)
# Point prevalence: (# points within AOH) / (# points within actual range)

packages <- c("sf", "terra", "dplyr", "readr", "tictoc", "foreach", "doParallel", "parallel")
lapply(packages, library, character.only = TRUE)

# keep consistent working directory style with other scripts
setwd('/')
setwd('C:/Users/wenxinyang/Desktop/GitHub/colander')
# setwd('Users/wenxinyang/Desktop/GitHub/colander')
source('scripts/refineBiomodelos/funcs.R')
source('scripts/refineBiomodelos/info.R')

dir_aoh <- 'data/col_aoh/2012_7gen_btst_keephab/'
dir_ranges <- 'data/IUCN_range_maps/cleaned_ranges'
pts_path_all <- 'data/occ_pts/allinfo_ideam_coords_all_0605.csv'
pts_path_temp <- 'data/occ_pts/allinfo_ideam_coords_2012_0605.csv'
out_dir <- 'results/validation'
out_csv <- file.path(out_dir, 'aoh_validation_0622.csv')

if(!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# initialize CSV with headers if it doesn't exist
if (!file.exists(out_csv)) {
  df_headers <- data.frame(
    threshold = character(0),
    species = character(0),
    taxa = character(0),
    range_area_km2 = numeric(0),
    aoh_area_km2 = numeric(0),
    model_prevalence = numeric(0),
    points_in_range = integer(0),
    points_in_aoh = integer(0),
    point_prevalence = numeric(0),
    error_message = character(0),
    stringsAsFactors = FALSE
  )
  readr::write_csv(df_headers, out_csv)
}

# check if species already exists in CSV
species_exists_in_csv <- function(species_name, csv_path) {
  if (!file.exists(csv_path)) return(FALSE)
  existing_df <- readr::read_csv(csv_path, show_col_types = FALSE)
  return(species_name %in% existing_df$species)
}

# append result to CSV
append_to_csv <- function(result_df, csv_path) {
  if (file.exists(csv_path)) {
    existing_df <- readr::read_csv(csv_path, show_col_types = FALSE)
    # Remove any existing entry for this species
    existing_df <- existing_df[existing_df$species != result_df$species, ]
    # Combine and write
    combined_df <- rbind(existing_df, result_df)
    readr::write_csv(combined_df, csv_path)
  } else {
    readr::write_csv(result_df, csv_path)
  }
}

# helper to safely read a vector layer if it exists
read_vector_if_exists <- function(path) {
  if (file.exists(path)) {
    return(sf::st_read(path, quiet = TRUE))
  }
  return(NULL)
}

# map Spanish taxa labels in points to range groups
taxa_map <- c(
  aves = 'Birds',
  mamiferos = 'Mammals',
  squamata = 'Reptiles',
  anfibios = 'Amphibians'
)

# return species range (sf) merged to a single geometry, or NULL if not found
get_species_range <- function(species_name, taxa_label, range_data) {
  if (is.null(range_data)) return(NULL)
  
  hit <- trySuppressWarnings(range_data[range_data$sci_name == species_name, ])
  if (nrow(hit) > 0) {
    sp_range <- sf::st_make_valid(hit)
    # Check if geometry is valid and fix if needed
    if(FALSE %in% sf::st_is_valid(sp_range)){
      sf::sf_use_s2(FALSE)
      sp_range <- sf::st_make_valid(sp_range)
    }
    return(sp_range)
  }
  return(NULL)
}

trySuppressWarnings <- function(expr) {
  suppressWarnings(suppressMessages(expr))
}

# read occurrence points from CSV
pts_df <- if (file.exists(pts_path_temp)) read.csv(pts_path_temp) else NULL
if (is.null(pts_df)) stop('Occurrence points file not found: ', pts_path_temp)
if('scntfcN' %in% colnames(pts_df)) pts_df <- pts_df %>% rename(species = scntfcN)
if('scientificName' %in% colnames(pts_df)) pts_df <- pts_df %>% rename(species = scientificName)
# ensure standard column names
if(!('species' %in% names(pts_df))) stop('Points file missing "species" attribute')
if(!('taxa' %in% names(pts_df))) pts_df$taxa <- NA_character_
if(!('longitude' %in% names(pts_df))) stop('Points file missing "longitude" attribute')
if(!('latitude' %in% names(pts_df))) stop('Points file missing "latitude" attribute')

# convert to sf object
pts_sf <- sf::st_as_sf(pts_df, coords = c("longitude", "latitude"), crs = 4326)

# species from available AOH rasters
aoh_files <- list.files(dir_aoh, pattern = '\\.(tif|tiff)$', full.names = TRUE)
if(length(aoh_files) == 0) stop('No AOH rasters found in ', dir_aoh)

sp_from_files <- basename(aoh_files)
sp_from_files <- gsub('\\.(tif|tiff)$', '', sp_from_files)
species_list <- gsub('_', ' ', sp_from_files)

# group species by taxa to process efficiently
species_by_taxa <- list()
for(sp_name in species_list) {
  sp_taxa <- unique(pts_sf[pts_sf$species == sp_name, ]$taxa)
  sp_taxa <- if(length(sp_taxa) > 0) sp_taxa[1] else NA_character_
  
  if(!is.na(sp_taxa) && sp_taxa %in% names(taxa_map)) {
    taxa_group <- taxa_map[[sp_taxa]]
    if(!(taxa_group %in% names(species_by_taxa))) {
      species_by_taxa[[taxa_group]] <- list()
    }
    species_by_taxa[[taxa_group]] <- append(species_by_taxa[[taxa_group]], sp_name)
  } else {
    # if taxa unknown, add to a catch-all group
    if(!("Unknown" %in% names(species_by_taxa))) {
      species_by_taxa[["Unknown"]] <- list()
    }
    species_by_taxa[["Unknown"]] <- append(species_by_taxa[["Unknown"]], sp_name)
  }
}

# Function to load range data for a taxa group
load_range_data <- function(taxa_group) {
  if(taxa_group == 'Mammals') {
    return(read_vector_if_exists(file.path(dir_ranges, 'Cleaned_MAMMALS.shp')))
  } else if(taxa_group == 'Birds') {
    return(read_vector_if_exists(file.path(dir_ranges, 'Cleaned_kindof_birds.gpkg')))
  } else if(taxa_group == 'Reptiles') {
    range_data1 <- read_vector_if_exists(file.path(dir_ranges, 'Cleaned_REPTILES_PART1.shp'))
    range_data2 <- read_vector_if_exists(file.path(dir_ranges, 'Cleaned_REPTILES_PART2.shp'))
    if(!is.null(range_data1) && !is.null(range_data2)) {
      return(rbind(range_data1, range_data2))
    } else if(!is.null(range_data1)) {
      return(range_data1)
    } else if(!is.null(range_data2)) {
      return(range_data2)
    }
  } else if(taxa_group == 'Amphibians') {
    range_data1 <- read_vector_if_exists(file.path(dir_ranges, 'Cleaned_AMPHIBIANS_PART1.shp'))
    range_data2 <- read_vector_if_exists(file.path(dir_ranges, 'Cleaned_AMPHIBIANS_PART2.shp'))
    if(!is.null(range_data1) && !is.null(range_data2)) {
      return(rbind(range_data1, range_data2))
    } else if(!is.null(range_data1)) {
      return(range_data1)
    } else if(!is.null(range_data2)) {
      return(range_data2)
    }
  }
  return(NULL)
}

# Function to process a single species
process_species <- function(sp_name, thresholds, range_data, out_csv, pts_sf, aoh_files, taxa_map) {
  tryCatch({
    # Check if species already exists in CSV
    if (species_exists_in_csv(sp_name, out_csv)) {
      message('Species ', sp_name, ' already processed; skipping')
      return(NULL)
    }
    
    message('Processing: ', sp_name)
    
    # Get species taxa
    sp_taxa <- unique(pts_sf[pts_sf$species == sp_name, ]$taxa)
    sp_taxa <- if(length(sp_taxa) > 0) sp_taxa[1] else NA_character_
    
    # Find AOH file
    aoh_path <- aoh_files[gsub(' ', '_', sp_name) == gsub('\\.(tif|tiff)$', '', basename(aoh_files))]
    if(length(aoh_path) == 0) {
      stop('AOH file not found for species: ', sp_name)
    }
    aoh_path <- aoh_path[1]
    
    # Get range polygons
    sp_range <- get_species_range(sp_name, sp_taxa, range_data)
    if (is.null(sp_range)) {
      stop('Range not found for ', sp_name)
    }
    
    # Dissolve to single geometry
    sp_range_u <- trySuppressWarnings(sf::st_union(sf::st_make_valid(sp_range)))
    
    # Compute range area (km^2)
    range_area_km2 <- tryCatch({
      rng3116 <- sf::st_transform(sp_range_u, 3116)
      as.numeric(sum(sf::st_area(rng3116))) / 1e6
    }, error = function(e) NA_real_)
    
    # Read AOH raster
    aoh <- tryCatch({ terra::rast(aoh_path) }, error = function(e) NULL)
    if (is.null(aoh)) {
      stop('Failed to read AOH for ', sp_name)
    }
    
    # Initialize results list
    results_list <- list()
    
    # Process each threshold
    for(threshold in thresholds) {
      # Initialize result
      result_df <- data.frame(
        threshold = threshold,
        species = sp_name,
        taxa = sp_taxa,
        range_area_km2 = range_area_km2,
        aoh_area_km2 = NA_real_,
        model_prevalence = NA_real_,
        points_in_range = NA_integer_,
        points_in_aoh = NA_integer_,
        point_prevalence = NA_real_,
        error_message = "",
        stringsAsFactors = FALSE
      )
      
      # Binarize AOH
      aoh_mask <- terra::app(aoh, fun = function(x) 
        ifelse(is.na(x), NA, ifelse(x >= threshold, 1, NA)))
      
      # Calculate AOH area
      cell_area_km2 <- terra::cellSize(aoh_mask, unit = 'km')
      aoh_area_km2 <- as.numeric(terra::global(terra::mask(cell_area_km2, aoh_mask), 
                                               'sum', na.rm = TRUE)[1,1])
      result_df$aoh_area_km2 <- aoh_area_km2
      
      # Model prevalence
      model_prev <- if(!is.na(range_area_km2) && range_area_km2 > 0) aoh_area_km2 / range_area_km2 else NA_real_
      result_df$model_prevalence <- model_prev
      
      # Point prevalence
      pts_sp <- pts_sf[pts_sf$species == sp_name, ]
      if (nrow(pts_sp) == 0) {
        num_pts_range <- 0
        num_pts_in_aoh <- 0
        point_prev <- NA_real_
      } else {
        # Points within range polygon
        pts_in_range <- trySuppressWarnings(pts_sp[sf::st_within(pts_sp, sp_range_u, sparse = FALSE)[,1], ])
        num_pts_range <- nrow(pts_in_range)
        if (num_pts_range == 0) {
          num_pts_in_aoh <- 0
          point_prev <- NA_real_
        } else {
          # Create 300m buffer around points
          pts_projected <- sf::st_transform(pts_in_range, 32618)
          pts_buffered <- sf::st_buffer(pts_projected, dist = 300)
          
          # Project to AOH CRS
          pts_buffered_aoh_crs <- sf::st_transform(pts_buffered, terra::crs(aoh))
          pts_buffered_v <- terra::vect(pts_buffered_aoh_crs)
          
          # Extract AOH values within buffered areas
          vals <- tryCatch({ 
            terra::extract(aoh, pts_buffered_v, fun = max, na.rm = TRUE)[,2] 
          }, error = function(e) rep(NA, num_pts_range))
          
          # Count points where buffered area intersects with AOH
          num_pts_in_aoh <- sum(!is.na(vals) & vals > 0)
          point_prev <- num_pts_in_aoh / num_pts_range
        }
      }
      
      result_df$points_in_range <- num_pts_range
      result_df$points_in_aoh <- num_pts_in_aoh
      result_df$point_prevalence <- point_prev
      
      results_list[[as.character(threshold)]] <- result_df
    }
    
    # Combine all thresholds and write to CSV
    combined_results <- do.call(rbind, results_list)
    
    # Write with lock to prevent conflicts
    if (file.exists(out_csv)) {
      # Use file lock with retry
      max_retries <- 5
      for (i in 1:max_retries) {
        tryCatch({
          existing_df <- readr::read_csv(out_csv, show_col_types = FALSE)
          existing_df <- existing_df[existing_df$species != sp_name, ]
          combined_df <- rbind(existing_df, combined_results)
          readr::write_csv(combined_df, out_csv)
          break
        }, error = function(e) {
          if (i == max_retries) stop(e)
          Sys.sleep(1)  # Wait before retry
        })
      }
    } else {
      readr::write_csv(combined_results, out_csv)
    }
    
    # Clean up
    rm(aoh, aoh_mask, cell_area_km2)
    gc()
    
    return(paste0("Success: ", sp_name))
    
  }, error = function(e) {
    # Record error
    error_msg <- paste0("Error processing ", sp_name, ": ", e$message)
    message(error_msg)
    
    # Create error result
    error_result <- data.frame(
      threshold = NA_character_,
      species = sp_name,
      taxa = NA_character_,
      range_area_km2 = NA_real_,
      aoh_area_km2 = NA_real_,
      model_prevalence = NA_real_,
      points_in_range = NA_integer_,
      points_in_aoh = NA_integer_,
      point_prevalence = NA_real_,
      error_message = e$message,
      stringsAsFactors = FALSE
    )
    
    tryCatch({
      append_to_csv(error_result, out_csv)
    }, error = function(e2) {
      message("Could not write error result: ", e2$message)
    })
    
    return(paste0("Failed: ", sp_name))
  })
}

# Main function with parallel processing
AOHValidationParallel <- function(thresholds, out_csv, n_cores = NULL) {
  # Determine number of cores
  if (is.null(n_cores)) {
    n_cores <- max(1, detectCores() - 1)  # Leave one core for system
  }
  
  message("Using ", n_cores, " cores for parallel processing")
  
  # Setup parallel backend
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  # Define thresholds as a variable in the global environment of each worker
  # This is the key fix for the 'thresholds' not found error
  clusterExport(cl, c("thresholds"), envir = environment())
  
  # Export all necessary variables to cluster
  clusterExport(cl, c(
    "pts_sf", "aoh_files", "out_csv", 
    "dir_ranges", "taxa_map", "species_by_taxa",
    "species_exists_in_csv", "get_species_range", 
    "load_range_data", "append_to_csv", "process_species",
    "read_vector_if_exists", "trySuppressWarnings"
  ))
  
  # Load packages on each worker
  clusterEvalQ(cl, {
    library(sf)
    library(terra)
    library(dplyr)
    library(readr)
  })
  
  # Process each taxa group in parallel
  results <- tryCatch({
    foreach(taxa_group = names(species_by_taxa), 
            .packages = c("sf", "terra", "dplyr", "readr"),
            .errorhandling = "pass") %dopar% {
              
              message("Processing taxa group: ", taxa_group)
              
              # Load range data for this taxon
              range_data <- load_range_data(taxa_group)
              if (is.null(range_data)) {
                warning('No range data found for taxa group: ', taxa_group)
                return(NULL)
              }
              
              # Process species in this taxa group
              species_list <- species_by_taxa[[taxa_group]]
              results_summary <- list()
              
              for(sp_name in species_list) {
                # Pass thresholds explicitly to the function
                result <- process_species(sp_name, thresholds, range_data, out_csv, 
                                          pts_sf, aoh_files, taxa_map)
                results_summary[[sp_name]] <- result
              }
              
              return(results_summary)
            }
  }, error = function(e) {
    message("Error in parallel processing: ", e$message)
    return(NULL)
  })
  
  # Stop cluster
  stopCluster(cl)
  
  return(results)
}

# Define thresholds
li_thresholds <- seq(800, 1000, by = 50)

# Run with parallel processing
tic()
message("Starting parallel processing...")
results <- AOHValidationParallel(li_thresholds, out_csv, n_cores = 4)
toc()

message('Validation complete. Results saved to: ', out_csv)

# Print summary
if (!is.null(results)) {
  message("Processing summary:")
  for (i in seq_along(results)) {
    taxa_group <- names(species_by_taxa)[i]
    if (!is.null(results[[i]])) {
      summary_stats <- table(unlist(results[[i]]))
      message("  ", taxa_group, ": ", paste(names(summary_stats), summary_stats, collapse = ", "))
    }
  }
}