packages <- c("sf", "terra", "here", "rgrass", "dplyr", "data.table", 
              "rredlist", "tidyr", "Matrix", "tidyverse", "pROC",
              "kableExtra", "ggplot2", "reshape2", "parallel")
lapply(packages, library, character.only = TRUE)

setwd('/')
setwd('Users/wenxinyang/Desktop/GitHub/colander')
source('scripts/refineBiomodelos/funcs.R')
source('scripts/refineBiomodelos/info.R')


dfpref <- read.csv('data/animals_preference.csv', sep=';') %>%
  mutate(habitats.code = as.character(habitats.code))
dfpref1 <- dfpref %>% mutate(
  habcode = case_when(
    habitats.code=='14.1' | habitats.code=='14.2'   ~ '14.12',
    habitats.code=='14.3' | habitats.code=='14.6'   ~ '14.36',
    habitats.code=='14.4' | habitats.code=='14.5'   ~ '14.45',
    str_detect(habitats.code, '15.') ~ '15',
    str_detect(habitats.code, '1.') ~ '1',
    str_detect(habitats.code, '2.') ~ '2',
    str_detect(habitats.code, '3.') ~ '3',
    str_detect(habitats.code, '4.') ~ '4',
    str_detect(habitats.code, '5.') ~ '5',
    str_detect(habitats.code, '6.') ~ '6',
    str_detect(habitats.code, '8.') ~ '8',
    TRUE ~ NA
  )
) %>% filter(!is.na(habcode))
unique(dfpref1$habcode)
dfpref1 <- dfpref1 %>% select(name, habcode) %>% unique()
rm(dfpref)
#df_pref <- read.csv('data/occ_pts/col_animal_pref_cleaned.csv')
#df_pref1 <- df
df_basic_info <- addGeneralistInfo(df_pref, 7)
# rm(df_pref)

df <- read.csv('results/validation/aoh_validation_0626.csv') 
as.datam.frame(table(df$threshold))
df <- df %>% filter(!is.na(point_prevalence) & !is.na(model_prevalence))
df_spp <- df %>% select(species, taxa) %>% unique()

table(df_spp$taxa) 
#aves: 1086; mamiferos: 145; anfibios: 67; squamata: 55
unique(df$threshold)

## ---------------- 1. threshold vs. validation accuracy/count -----------------
dfthrshval <- as.data.frame(matrix(data=NA, nrow=0, ncol=3))
colnames(dfthrshval) <- c('taxa', 'threshold0', 'count_passing')
for (t in c('aves', 'mamiferos', 'anfibios', 'squamata', 'all')){
  cat(t, "\n")
  for(i in seq(500, 1000, 50)){
    if(t != 'all'){
      tmp <- nrow(df %>% filter(taxa == t & threshold == i & point_prevalence>model_prevalence))
      } else {
        tmp <- nrow(df %>% filter(threshold == i & point_prevalence>model_prevalence))
      }
    cat(i, ":", tmp, "\n")
    dfthrshval[nrow(dfthrshval)+1, ] <- c(t, as.integer(i), as.integer(tmp))
    }
}

dfthrshval <- dfthrshval %>% 
  mutate(
    threshold0 = factor(threshold0, levels = sort(unique(as.numeric(threshold0)), decreasing = FALSE))
  )

### --------------- 1.1. threshold vs. validation accuracy/count ---------------
ggplot(dfthrshval, aes(x=threshold0, y=count_passing, group=taxa))+
  geom_point()+
  geom_line()+
  facet_wrap(~taxa, scales="free")

all750 <- df %>% filter(threshold==750 & model_prevalence<point_prevalence)
all800 <- df %>% filter(threshold==800 & model_prevalence<point_prevalence)
#all850 <- df %>% filter(threshold==850 & model_prevalence<point_prevalence)
#all900 <- df %>% filter(threshold==900 & model_prevalence<point_prevalence)
li1 <- setdiff(unique(all750$species), unique(all800$species))
li2 <- setdiff(unique(all800$species), unique(all750$species))

# ruled out generalist spp. problem
genli1 <- df_basic_info %>% filter(name %in% li1)
table(genli1$type)
genli2 <- df_basic_info %>% filter(name %in% li2)
table(genli2$type)

### --------------- 1.2. threshold vs. points -------------------
ggplot(df, aes(x=point_prevalence, y=model_prevalence, group=threshold))+
  geom_point()+
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed")+
  facet_wrap(~threshold)



## ---------------- 2. threshold vs. sum of aoh area -----------------
hist(df$range_area_km2)
summary(df$range_area_km2)
nrow(df %>% filter(range_area_km2>3*10^7))
ggplot(df %>% filter(range_area_km2 <3*10^7), aes(x=range_area_km2, y=aoh_area_km2, color=as.factor(threshold)))+
  geom_point() +
  facet_wrap(~threshold)


dfsuminfo <- df %>% group_by(threshold, taxa) %>% summarize(
  sumptprv = sum(point_prevalence),
  summdprv = sum(model_prevalence),
  sumaoharea = sum(aoh_area_km2)
) %>% mutate(diffprv = sumptprv-summdprv)


allsuminfo <- df %>% group_by(threshold) %>% summarize(
  sumptprv = sum(point_prevalence),
  summdprv = sum(model_prevalence),
  sumaoharea = sum(aoh_area_km2)
) %>% mutate(taxa = 'all',
             diffprv = sumptprv-summdprv)

dfsuminfo_all <- rbind(dfsuminfo, allsuminfo)

### ----------------- 2.1. plot threshold vs sum of aoh area ------------------
ggplot(dfsuminfo_all %>% filter(threshold!=1000), aes(x=threshold, y=sumaoharea))+
  geom_point()+
  geom_line()+
  facet_wrap(~taxa, scales="free")

### ---------------- 2.2. compute change rate of sum of aoh area ---------------
dfdaoh <- as.data.frame(matrix(data=NA, nrow=0, ncol=4))
colnames(dfdaoh) <- c('taxa', 'threshold', 'delta','rate')
for(t in c('aves', 'mamiferos', 'anfibios', 'squamata', 'all')){
  for(i in seq(800, 900, 50)){
    area_0 <- dfsuminfo_all %>% filter(taxa == t & threshold == i) %>% pull(sumaoharea)
    area_1 <- dfsuminfo_all %>% filter(taxa == t & threshold == i+50) %>% pull(sumaoharea)
    delta = area_0-area_1
    rate = delta/area_0
    dfdaoh[nrow(dfdaoh)+1, ] <- c(t, i, delta, rate)
  }
}

ggplot(dfdaoh %>% filter(threshold!=1000), aes(x=threshold, y=rate))+
  geom_point()+
  facet_wrap(~taxa, scales="free")

## --------------------- 3. number of matched pairs ---------------------
mat <- read.csv('results/glm_btst_2012_keephab/gen7_glm_btst_2012_keephab/count_above_1_matrix.csv')
dfcountall <- as.data.frame(matrix(data=NA, nrow=0, ncol=2))
colnames(dfcountall) <- c('threshold', 'count')
for(i in seq(500,1000,50)){
  cat('threshold is', i, '\n')
  a <- mat %>%
    select(where(is.numeric)) %>%
    summarise(across(everything(), ~ sum(. >= i, na.rm = TRUE)))
  print(a)
  print(sum(a))
  dfcountall[nrow(dfcountall)+1, ] <- c(i, sum(a))
}
ggplot(dfcountall, aes(x=threshold, y=count))+
  geom_point()+
  geom_line()+
  geom_label(aes(label=count))


## ----------------------- make final matrix figure ---------------------
pos_data <- read.csv('results/glm_btst_2012_keephab/gen7_glm_btst_2012_keephab/btst_ideam_randomCI_gen7__0615_2026_pos.csv')
good_col_names <- gsub('\\.', ' ', colnames(pos_data))
good_col_names <- gsub("Wetlands  inland", "Wetland (inland)", good_col_names)
colnames(pos_data) <- good_col_names
  
if("X" %in% colnames(pos_data)){
  pos_data$X <- NULL
  }
  
# Get habitat columns (exclude land_cover, auc, n_samples)
habitat_cols <- colnames(pos_data)[!colnames(pos_data) %in% c("", "land_cover", "auc", "n_samples")]
  
# Convert to long format for plotting
pos_long <- pos_data %>%
  select(land_cover, all_of(habitat_cols)) %>%
  tidyr::gather(habitat, count_value, -land_cover)
  
# Create color categories for habitat columns
pos_long$color_category <- cut(pos_long$count_value, 
                               breaks = c(-Inf, 850, Inf),
                               labels = c("Not a pair (≤850)", "AOH (>850)"),
                               include.lowest = TRUE)
  
# Add AUC and n_samples data
auc_data <- data.frame(
  land_cover = pos_data$land_cover,
  habitat = "AUC",
  count_value = pos_data$auc,
  color_category = "AUC"
)
  
n_samples_data <- data.frame(
  land_cover = pos_data$land_cover,
  habitat = "n_samples",
  count_value = round(pos_data$n_samples),
  color_category = "n_samples"
)
  
# Combine all data
plot_data <- rbind(pos_long, auc_data, n_samples_data)
  
# Reorder color categories for legend: Low at top, Not a pair at bottom
plot_data$color_category <- factor(plot_data$color_category,
                                     levels = c("AOH (>850)", "Not a pair (≤850)", "AUC", "n_samples"))
  
# Order land cover types alphabetically in descending order (Z to A)
plot_data$land_cover <- factor(plot_data$land_cover, 
                               levels = sort(unique(plot_data$land_cover), decreasing = TRUE))
  
# Set factor levels for habitat to include AUC and n_samples at the end
habitat_levels <- c(habitat_cols, "AUC", "n_samples")
plot_data$habitat <- factor(plot_data$habitat, levels = habitat_levels)
  
  # Create heatmap
p <- ggplot2::ggplot(plot_data, aes(x = habitat, y = land_cover, fill = color_category)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(color_category == "AUC", 
                               sprintf("%.2f", count_value),
                               ifelse(color_category == "n_samples",
                                      as.character(count_value),
                                      as.character(count_value)))), 
            size = 5, fontface = "bold", color = "black", family = "Arial") +
  scale_fill_manual(
    values = c("Not a pair (≤850)" = "#f7f7f7", "AOH (>850)" = "#659c6a",
               "AUC" = "white", "n_samples" = "white"),
    name = "Certainty level (counts)"
    ) +
  scale_x_discrete(position = "bottom") +
  theme_minimal() +
  theme(
    text = element_text(family = "Arial"),
    axis.text.x = element_text(size = 12, family = "Arial", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 15, family = "Arial"),
    axis.title.x = element_text(size = 16, face = 'bold', family = "Arial"),
    axis.title.y = element_text(size = 16, face = 'bold', family = "Arial"),
    legend.title = element_text(size = 16, family = "Arial"),
    legend.text = element_text(size = 14, family = "Arial"),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold", family = "Arial"),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA)
    ) +
  labs(
    # title = "Habitat-land cover translation matrix",
    x = "Habitat Classes",
    y = "Land Cover Classes"
    ) +
  coord_fixed(ratio = 0.6)
p
