# ======= prep habitat info ========
# create a dataframe to store the habitat codes and habitat names
habitat_codes <- c(1,2,3,4,5,6,7,8,9,10,12,13,14.1, 14.2,14.3, 14.4,14.5,14.6,15,16,17,18)
# concatenate "hab_" to the habitat codes
habitat_codes <- paste0("hab_", habitat_codes)
habitat_names <- c("Forest", "Savanna", "Shrubland", "Grassland", "Wetlands (inland)", "Rocky areas", 
                   "drop", #"Caves and Subterranean",
                   "Desert",
                   "drop","drop","drop","drop", #"Marine Neritic", "Marine Oceanic", "Marine Intertidal", "Marine Coastal/Supratidal", 
                   "Artificial-arable and pasture", "Artificial-arable and pasture", #"Artificial-Arable Land", "Artificial-Pastureland",
                   "Artificial-degraded forest and plantation", # "Artificial-Plantations", merged with degraded forest
                   "Artificial-urban areas and rural gardens", "Artificial-urban areas and rural gardens", # "Artificial-Rural Gardens", "Artificial-Urban Areas",
                   "Artificial-degraded forest and plantation", #"Artificial-Degraded Former Forest", merged with plantations
                   "Artificial-aquatic", 
                   "drop","drop","drop"# "Introduced Vegetation", "Other", "Unknown"
)
# create a dataframe to store the habitat codes and habitat names
habitat_info <- data.frame(habitat_code = habitat_codes,
                           habitat_name = habitat_names)
drop_cols <- habitat_info[habitat_info$habitat_name=='drop',]$habitat_code

# cleaned up habitat info
habitat_info1 <- habitat_info[habitat_info$habitat_name!='drop',]
habitat_info1[nrow(habitat_info1)+1,] <- c('hab_14.12', 'Artificial-arable and pasture')
habitat_info1[nrow(habitat_info1)+1,] <- c('hab_14.36', 'Artificial-degraded forest and plantation')
habitat_info1[nrow(habitat_info1)+1,] <- c('hab_14.45', 'Artificial-urban areas and rural gardens')
habitat_info1 <- habitat_info1 %>% filter(!habitat_code %in% c('hab_14.1', 'hab_14.2', 'hab_14.3', 'hab_14.4', 'hab_14.5', 'hab_14.6'))



# ======== prep ideam land cover info =======
ideam_lc_code <- unique(read.csv('data/occ_pts/allinfo_ideam_cgls.csv')$nvl_2_n)
ideam_lc_code <- paste0("lc_", ideam_lc_code)
ideam_lc_names <- c("Forest", "Heterogeneous Ag", "Grass/Shrub", "Pasture", "Wetlands", "Artificial Green Space", "Permanent Crops", "Industrial/Commercial Zones", "Urban", "Open Areas", "Continental Waters", "Transition Crops", "Coastal Marshes", "Aquaculture", "Mine and Waste Dump")
ideam_lc_info <- data.frame(ideam_lc_code = ideam_lc_code,
                            ideam_lc_name = ideam_lc_names)

ideam_lc_code_lv1 <- unique(read.csv('data/occ_pts/allinfo_ideam_cgls.csv')$nvl_1_n)
ideam_lc_code_lv1 <- paste0("lc_", ideam_lc_code_lv1)
ideam_lc_names_lv1 <- c("Forest and semi-natural areas", "Agricultural land", "Wetlands", "Artificial land", "Water surfaces")
ideam_lc_info_lv1 <- data.frame(ideam_lc_code = ideam_lc_code_lv1,
                            ideam_lc_name = ideam_lc_names_lv1)


# ======== prep cgls land cover info =======
cgls_lc_code <- unique(read.csv('data/occ_pts/allinfo_ideam_cgls.csv')$cgls_value)
# order cgls_lc_code ascending
cgls_lc_code <- sort(cgls_lc_code)
cgls_lc_code <- paste0('lc_', cgls_lc_code)

cgls_lc_names <- c("Unknown", "Shrubs", "Herbaceous vegetation","Cultivated and managed vegetation/ag","Urban/built-up", "Bare/sparse vegetation", "Snow and ice", "Permanent water bodies", "Herbaceous wetland", "Closed forest, evergreen broadleaf", "Closed forest, deciduous broadleaf", "Closed forest, other", "Open forest, evergreen broadleaf", "Open forest, deciduous broadleaf", "Open forest, other", "Oceans and seas")
cgls_lc_info <- data.frame(cgls_lc_code = cgls_lc_code,
                           cgls_lc_name = cgls_lc_names)
