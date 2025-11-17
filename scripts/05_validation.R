# Compute model and point prevalence for each species
# Author: Wenxin Yang
# Date: September, 2025

# Model prevalence: AOH size / range size (km^2)
# Point prevalence: (# points within AOH) / (# points within actual range)

packages <- c("sf", "terra", "dplyr", "readr")
lapply(packages, library, character.only = TRUE)

# keep consistent working directory style with other scripts
setwd('/')
setwd('Users/wenxinyang/Desktop/GitHub/colander')

dir_aoh <- 'data/dsc_aoh'
dir_ranges <- 'data/IUCN_range_maps/cleaned_ranges'
sp_list_path <- 'data/CI_STAR_T_Species_List.csv'
pts_path <- 'data/occ_pts/allinfo_ideam_cgls_coords_2022.csv'
out_dir <- 'results/validation'
out_csv <- file.path(out_dir, 'aoh_validation_300m.csv')

if(!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# initialize CSV with headers if it doesn't exist
if (!file.exists(out_csv)) {
  df_headers <- data.frame(
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
pts_df <- if (file.exists(pts_path)) read.csv(pts_path) else NULL
if (is.null(pts_df)) stop('Occurrence points file not found: ', pts_path)

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

# process each taxon group separately to save memory
for(taxa_group in names(species_by_taxa)) {
  message('Processing taxa group: ', taxa_group)
  
  # load range data for this taxon only
  range_data <- NULL
  if(taxa_group == 'Mammals') {
    range_data <- read_vector_if_exists(file.path(dir_ranges, 'Cleaned_MAMMALS.shp'))
  } else if(taxa_group == 'Birds') {
    range_data <- read_vector_if_exists(file.path(dir_ranges, 'Cleaned_kindof_birds.gpkg'))
  } else if(taxa_group == 'Reptiles') {
    range_data1 <- read_vector_if_exists(file.path(dir_ranges, 'Cleaned_REPTILES_PART1.shp'))
    range_data2 <- read_vector_if_exists(file.path(dir_ranges, 'Cleaned_REPTILES_PART2.shp'))
    if(!is.null(range_data1) && !is.null(range_data2)) {
      range_data <- rbind(range_data1, range_data2)
    } else if(!is.null(range_data1)) {
      range_data <- range_data1
    } else if(!is.null(range_data2)) {
      range_data <- range_data2
    }
  } else if(taxa_group == 'Amphibians') {
    range_data1 <- read_vector_if_exists(file.path(dir_ranges, 'Cleaned_AMPHIBIANS_PART1.shp'))
    range_data2 <- read_vector_if_exists(file.path(dir_ranges, 'Cleaned_AMPHIBIANS_PART2.shp'))
    if(!is.null(range_data1) && !is.null(range_data2)) {
      range_data <- rbind(range_data1, range_data2)
    } else if(!is.null(range_data1)) {
      range_data <- range_data1
    } else if(!is.null(range_data2)) {
      range_data <- range_data2
    }
  }
  
  if(is.null(range_data)) {
    warning('No range data found for taxa group: ', taxa_group)
    next
  }
  
  # process each species in this taxon group
  for(sp_name in species_by_taxa[[taxa_group]]) {
    # Check if species already exists in CSV
    if (species_exists_in_csv(sp_name, out_csv)) {
      message('Species ', sp_name, ' already processed; skipping')
      next
    }
    
    message('Processing: ', sp_name)
    
    # initialize result with NA values and empty error message
    result_df <- data.frame(
      species = sp_name,
      taxa = NA_character_,
      range_area_km2 = NA_real_,
      aoh_area_km2 = NA_real_,
      model_prevalence = NA_real_,
      points_in_range = NA_integer_,
      points_in_aoh = NA_integer_,
      point_prevalence = NA_real_,
      error_message = "",
      stringsAsFactors = FALSE
    )
    
    # Wrap entire processing in tryCatch for error handling
    tryCatch({
      aoh_path <- aoh_files[gsub(' ', '_', sp_name) == gsub('\\.(tif|tiff)$', '', basename(aoh_files))]
      if(length(aoh_path) == 0) {
        stop('AOH file not found for species: ', sp_name)
      }
      aoh_path <- aoh_path[1]  # take first match
      
      sp_taxa <- unique(pts_sf[pts_sf$species == sp_name, ]$taxa)
      sp_taxa <- if(length(sp_taxa) > 0) sp_taxa[1] else NA_character_
      result_df$taxa <- sp_taxa

      # get range polygons
      sp_range <- get_species_range(sp_name, sp_taxa, range_data)
      if (is.null(sp_range)) {
        stop('Range not found for ', sp_name)
      }
      # dissolve to single geometry for robust containment tests
      sp_range_u <- trySuppressWarnings(sf::st_union(sf::st_make_valid(sp_range)))

      # compute range area (km^2) using EPSG:3116 (meters)
      range_area_km2 <- tryCatch({
        rng3116 <- sf::st_transform(sp_range_u, 3116)
        as.numeric(sum(sf::st_area(rng3116))) / 1e6
      }, error = function(e) NA_real_)
      result_df$range_area_km2 <- range_area_km2

      # read AOH raster and compute area (km^2)
      aoh <- tryCatch({ terra::rast(aoh_path) }, error = function(e) NULL)
      if (is.null(aoh)) {
        stop('Failed to read AOH for ', sp_name)
      }
      # non-NA cells are considered AOH presence
      aoh_mask <- terra::app(aoh, fun = function(x) ifelse(is.na(x), NA, 1))
      cell_area_km2 <- terra::cellSize(aoh_mask, unit = 'km')
      aoh_area_km2 <- as.numeric(terra::global(terra::mask(cell_area_km2, aoh_mask), 'sum', na.rm = TRUE)[1,1])
      result_df$aoh_area_km2 <- aoh_area_km2

      model_prev <- if(!is.na(range_area_km2) && range_area_km2 > 0) aoh_area_km2 / range_area_km2 else NA_real_
      result_df$model_prevalence <- model_prev

      # point prevalence
      pts_sp <- pts_sf[pts_sf$species == sp_name, ]
      if (nrow(pts_sp) == 0) {
        num_pts_range <- 0
        num_pts_in_aoh <- 0
        point_prev <- NA_real_
      } else {
        # points within range polygon
        pts_in_range <- trySuppressWarnings(pts_sp[sf::st_within(pts_sp, sp_range_u, sparse = FALSE)[,1], ])
        num_pts_range <- nrow(pts_in_range)
        if (num_pts_range == 0) {
          num_pts_in_aoh <- 0
          point_prev <- NA_real_
        } else {
          # create 300m buffer around points and check for AOH intersection
          # project points to a projected CRS for accurate buffering (using UTM zone 18N for Colombia)
          pts_projected <- sf::st_transform(pts_in_range, 32618)  # UTM 18N
          pts_buffered <- sf::st_buffer(pts_projected, dist = 300)  # 300m buffer
          
          # project buffered points back to AOH raster CRS
          pts_buffered_aoh_crs <- sf::st_transform(pts_buffered, terra::crs(aoh))
          pts_buffered_v <- terra::vect(pts_buffered_aoh_crs)
          
          # extract AOH values within buffered areas
          vals <- tryCatch({ 
            terra::extract(aoh, pts_buffered_v, fun = max, na.rm = TRUE)[,2] 
          }, error = function(e) rep(NA, num_pts_range))
          
          # count points where buffered area intersects with AOH (any non-NA, non-zero value)
          num_pts_in_aoh <- sum(!is.na(vals) & vals > 0)
          point_prev <- num_pts_in_aoh / num_pts_range
        }
      }
      
      result_df$points_in_range <- num_pts_range
      result_df$points_in_aoh <- num_pts_in_aoh
      result_df$point_prevalence <- point_prev
      
    }, error = function(e) {
      # Record error message
      result_df$error_message <<- as.character(e$message)
      message('Error processing ', sp_name, ': ', e$message)
    })
    
    # Write result to CSV after each species
    append_to_csv(result_df, out_csv)
    message('Saved result for ', sp_name, ' to CSV')
  }
  
  # clear range data from memory after processing each taxon
  rm(range_data)
  gc()
}

message('Validation complete. Results saved to: ', out_csv)


