# RGB Kew MSc Spatial Analysis Course 2024
# Joseph White, Carolina Tovar, Tarciso Leao, Felix Lim
# 2024/11/27

install.packages('rgbif')
install.packages('CoordinateCleaner')
install.packages('rnaturalearthdata')
install.packages('tidysdm')
install.packages('tidyterra')
install.packages('overlapping')
install.packages('ranger')
install.packages('xgboost')
install.packages('DALEX')

# Load in libraries
library(rgbif)
library(CoordinateCleaner)
library(rnaturalearth)
library(rnaturalearthdata)
library(tidyverse)
library(tidysdm)
library(tidyterra)
library(terra)
library(MASS)
library(patchwork)
library(sf)
library(DALEX)

### 1. Load in data ----
#### 1a. Country boundary ----
mad <- vect(ne_countries(country = 'Madagascar', scale = 'medium'))
crs(mad) <- "lonlat"

#### 1a. occurrence data ----
# Name your species
species_name <- c("Aerangis citrata")

# download GBIF occurrence data for this species; this takes time if there are many data points!
gbif_download <- occ_data(scientificName = species_name, hasCoordinate = TRUE, country = 'MG', limit = 10000)

# select only the data (ignore the other metadata for now)
gbif_data <- gbif_download$data
head(gbif_data)

#### 1b. environmental data ----
rast_files <- list.files(paste0(getwd(),'/data/environmental_data/'), full.names = T, recursive = T)
env_vars <- app(sds(rast_files),c)
names(env_vars) <- c("mean_ann_t",'mean_t_warm_q','mean_t_cold_q','ann_p', 'p_wet_m','p_dry_m','p_seas','p_wet_q','p_dry_q','p_warm_q','p_cold_q',"mean_diurnal_t_range","isothermality", "t_seas", 'max_t_warm_m','min_t_cold_m',"t_ann_range",'mean_t_wet_q','mean_t_dry_q', 'elev')
names(env_vars)

### 2. Clean occurrence data ----
# keep the necessary columns and rename them
spp_raw <- gbif_data %>% 
              dplyr::select(gbifID, scientificName, latitude = decimalLatitude, longitude = decimalLongitude)

# CoordinateCleaner
spp_flag <- spp_raw %>% 
                clean_coordinates(
                  species = 'scientificName',
                  lon = "longitude",
                  lat = "latitude",
                  value = "spatialvalid")
head(spp_flag)

# filter to only include valid records
spp <- spp_flag %>%
          filter(.summary == TRUE) %>%
          dplyr::select(gbifID:longitude)

# spatial thin
set.seed(1234567)
spp <- thin_by_cell(spp, raster = env_vars[[1]])
nrow(spp)

# thin by distance
set.seed(1234567)
spp_thin <- thin_by_dist(spp, dist_min = km2m(5))
nrow(spp_thin)

# convert to sf 
spp_thin <- st_as_sf(spp_thin, coords = c('longitude', 'latitude'), crs = 'EPSG:4326')

# plot output
ggplot() +
  geom_spatraster(data = env_vars, aes(fill = ann_p)) +
  scale_fill_cross_blended_c(palette = 'arid', direction = -1) +
  geom_sf(data = spp_thin) + 
  theme_void() +
  guides(fill="none")

### 3. Generate background/pseudo-absence points ----
# download GBIF occurrence data for the genus
genus_name <- word(species_name, 1)
gbif_download_genus <- occ_data(scientificName = genus_name, hasCoordinate = TRUE, country = 'MG', limit = 10000)
gbif_data_genus <- gbif_download_genus$data
# remove our model species
gbif_data_genus <- gbif_data_genus %>% filter(!species %in% species_name)
table(gbif_data_genus$scientificName)

# keep the necessary columns
genus_raw <- gbif_data_genus %>% 
                      dplyr::select(gbifID, 
                                    latitude = decimalLatitude, 
                                    longitude = decimalLongitude)

# convert to spatvect
genus_vect <- vect(genus_raw, geom = c('longitude', 'latitude'), crs = 'EPSG:4326')

# taxon targeted collection density
genus_density <- rasterize(genus_vect, env_vars[[1]], 
                                    fun = "count")

output_resolution <- res(env_vars)[1]*100
df <- as.data.frame(genus_density, xy = T)
kde <- kde2d(df$x, df$y,
             n = c(nrow(env_vars), ncol(env_vars)),
             h = rep(output_resolution,2))
df = expand.grid(kde$x, kde$y, KEEP.OUT.ATTRS = FALSE)
df$z = as.vector(kde$z)
genus_density_rast <- rast(df, type = 'xyz', crs = crs(env_vars))
genus_density_rast <- project(genus_density_rast, env_vars[[1]], method = 'bilinear')
genus_density_rast <- crop(genus_density_rast, env_vars[[1]], mask = T)
plot(genus_density_rast)

# sample background with taxon targeted collection density weights
set.seed(1234567)
spp_all <- sample_background(data = spp_thin, 
                              raster = genus_density_rast,
                              n = 3 * nrow(spp_thin),
                              method = "bias",
                              class_label = "background",
                              return_pres = TRUE)

set.seed(1234567)
spp_all <- sample_pseudoabs(data = spp_thin, 
                             raster = genus_density_rast,
                             n = 1 * nrow(spp_thin),
                             method = c('dist_min', 10000),
                             class_label = "pseudoabs",
                             return_pres = TRUE)

# plot output
ggplot() +
  geom_spatraster(data = env_vars, aes(fill = ann_p)) +
  scale_fill_cross_blended_c(palette = 'arid', direction = -1) +
  geom_sf(data = spp_all[spp_all$class == 'presence',]) + 
  theme_void() +
  guides(fill="none") +

ggplot() +
  geom_spatraster(data = env_vars, aes(fill = ann_p)) +
  scale_fill_cross_blended_c(palette = 'arid', direction = -1) +
  geom_sf(data = spp_all, aes(col = class)) + 
  theme_void() +
  guides(fill="none")

### 4. Process environmental variables
spp_all <- spp_all %>%
  bind_cols(terra::extract(env_vars, spp_all, ID = FALSE))

# plot differences between env_vars X class
spp_all %>%
  plot_pres_vs_bg(class)

# find distances between density distributions of env_vars X class
dist_env_vars <- spp_all %>%
    dist_pres_vs_bg(class)

# suggested variables based on distance discrimination
top8_vars <- names(dist_env_vars[1:8])

# view pairs
pairs(env_vars[[top8_vars]], maxcells = 10000)

# identify multi-collinearity in environmental variables
vars_uncor <- filter_collinear(env_vars[[top5_vars]], cutoff = 0.8, method = "cor_caret")
vars_uncor

# clean datasets based on filtered environmental variables
spp_all <- spp_all %>% dplyr::select(all_of(c(vars_uncor, "class")))
env_vars <- env_vars[[vars_uncor]]

# plot differences between env_vars X class
spp_all %>%
  plot_pres_vs_bg(class)

### 4. Spatial cross-validation design ----
n_blocks <- 5
set.seed(100)
# cluster points using kmeans into n blocks
spp_cv <- spatial_clustering_cv(data = spp_all, 
                                cluster_function = 'kmeans',
                                v = n_blocks)

# view folds
autoplot(spp_cv, aes(shape = class)) +
  geom_spatvector(data = mad, fill = NA) +
  theme_void()

### 5. Run distribution models ----
# create the model formual
spp_rec <- recipe(spp_all, formula = class ~ .)

# design the model workflows
spp_models <- workflow_set(
  preproc = list(default = spp_rec),
  models = list(maxent = sdm_spec_maxent(),
                rf = sdm_spec_rf(),
                gbm = sdm_spec_boost_tree()))  %>%
  # tweak controls to store information needed later to create the ensemble
  option_add(control = control_ensemble_grid())

# run the models
set.seed(1234567)
spp_models <- spp_models %>%
                workflow_map('tune_grid',
                              resamples = spp_cv,
                              grid = 5,
                              metrics = sdm_metric_set(), 
                              verbose = T)

# show mean metric values for each model configuration
spp_models %>% collect_metrics()
autoplot(spp_models)

# extract hyper-parameters and metrics for individual models
model_results <- spp_models$result
names(model_results) <- spp_models$wflow_id
model_parameters <- map_dfr(model_results, ~.x, .id = 'model') %>%
  dplyr::select(model, id, .metrics) %>%
  unnest(.metrics)
model_parameters

# create the model ensemble by selecting the best model using boyce_cont
spp_ensemble <- simple_ensemble() %>%
  add_member(spp_models, metric = "boyce_cont")
autoplot(spp_ensemble)

### 6. Predict species suitability ----
# use a metric threshold e.g. boyce_cont >= 0.4
prediction_present_boyce <- predict_raster(
                spp_ensemble, 
                env_vars,
                metric_thresh = c("boyce_cont", 0.4),
                fun = "median")

# view prediction
ggplot() +
  geom_spatraster(data = prediction_present_boyce, 
                  aes(fill = median)) +
  scale_fill_whitebox_c(palette = 'high_relief', direction = -1,
                        name = 'Habitat\nsuitability',
                        limits = c(0,1)) +
  geom_sf(data = spp_all[spp_all$class == 'presence',], 
          size = 0.5, alpha = 0.25) +
  theme_void()

# select threshold value
spp_ensemble <- calib_class_thresh(
                    spp_ensemble,
                    class_thresh = c('sensitivity', 0.75), 
                    metric_thresh = c("boyce_cont", 0.4))

# create a binary map
prediction_present_boyce_binary <- predict_raster(
  spp_ensemble, 
  env_vars,
  type = "class",
  class_thresh = c('sensitivity', 0.75), 
  metric_thresh = c("boyce_cont", 0.4),
  fun = "median")

# view prediction
ggplot() +
  geom_spatraster(data = prediction_present_boyce_binary, 
                  aes(fill = binary_median)) +
  scale_fill_manual(values = c('orange', 'gray90'), 
                    na.value = 'transparent',
                    name = 'Habitat\nsuitability',
                    labels = c('present', 'absent')) +
  geom_sf(data = spp_all[spp_all$class == 'presence',], 
          size = 0.5, alpha = 0.25, col = 'purple') +
  theme_void()

### 7. Variable importance ----
explain_models <- explain_tidysdm(spp_ensemble)
vi <- model_parts(explain_models)

# view variable importance
vi %>%
  filter(variable %in% names(env_vars)) %>%
  ggplot() +
  geom_boxplot(aes(x = reorder(variable, dropout_loss), y = dropout_loss),
               fill = 'lightblue') +
  coord_flip() +
  labs(x = '', y = '1 - AUC loss after permutations') +
  theme_bw()

### 8. Partial dependence plots ----
pdp <- model_profile(explain_models, variables = names(env_vars))
agg_data <- as.data.frame(pdp$agr_profiles)

agg_data %>%
  filter(`_vname_` == names(env_vars)[1]) %>%
  ggplot(aes(x = `_x_`, y = `_yhat_`)) + 
  geom_line(lwd = 1, col = 'lightblue') +
  labs(x = names(env_vars)[1], y = 'mean prediction') +
  scale_y_continuous(limits = c(0, max(agg_data$`_yhat_`))) +
  theme_bw() +

agg_data %>%
  filter(`_vname_` == names(env_vars)[2]) %>%
  ggplot(aes(x = `_x_`, y = `_yhat_`)) + 
  geom_line(lwd = 1, col = 'lightblue') +
  labs(x = names(env_vars)[2], y = 'mean prediction') +
  scale_y_continuous(limits = c(0, max(agg_data$`_yhat_`))) +
  theme_bw() +

agg_data %>%
  filter(`_vname_` == names(env_vars)[3]) %>%
  ggplot(aes(x = `_x_`, y = `_yhat_`)) + 
  geom_line(lwd = 1, col = 'lightblue') +
  labs(x = names(env_vars)[3], y = 'mean prediction') +
  scale_y_continuous(limits = c(0, max(agg_data$`_yhat_`))) +
  theme_bw() +

agg_data %>%
  filter(`_vname_` == names(env_vars)[4]) %>%
  ggplot(aes(x = `_x_`, y = `_yhat_`)) + 
  geom_line(lwd = 1, col = 'lightblue') +
  labs(x = names(env_vars)[4], y = 'mean prediction') +
  scale_y_continuous(limits = c(0, max(agg_data$`_yhat_`))) +
  theme_bw() +
  plot_annotation(tag_levels = 'a', tag_suffix = ')')
