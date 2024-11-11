# RGB Kew MSc Spatial Analysis Course 2024
# Joseph White, Carolina Tovar, Moabe Fernandes, Felix Lim
# 2024/11/27

#### Load in libraries ----
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
library(exactextractr)
source('scripts/helper_functions/mtp.R')

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
env_vars
plot(env_vars[[1]], main = names(env_vars)[[1]])

### 2. Clean occurrence data ----
# keep the necessary columns and rename them
spp_raw <- gbif_data %>% 
              dplyr::select(gbifID, scientificName, latitude = decimalLatitude, longitude = decimalLongitude)

# CoordinateCleaner
spp_raw <- spp_raw %>% 
                clean_coordinates(
                  species = 'scientificName',
                  lon = "longitude",
                  lat = "latitude",
                  value = "spatialvalid")
head(spp_raw)

# filter to only include valid records
spp_clean <- spp_raw %>%
          filter(.summary == TRUE) %>%
          dplyr::select(gbifID:longitude)

# spatial thin
set.seed(1234567)
spp_thin <- thin_by_cell(spp_clean, raster = env_vars[[1]])
nrow(spp)

# convert to sf 
spp_thin <- st_as_sf(spp_thin, coords = c('longitude', 'latitude'), crs = 'EPSG:4326')

# plot output
ggplot() +
  geom_spatraster(data = env_vars[[1]]) +
  scale_fill_cross_blended_c(palette = 'arid', direction = -1) +
  geom_sf(data = spp_thin) + 
  theme_void() +
  guides(fill="none")

### 3. Generate background/pseudo-absence points ----
set.seed(1)
spp_all <- sample_pseudoabs(data = spp_thin, 
                             raster = env_vars,
                             n = 1 * nrow(spp_thin),
                             method = c('dist_min', 10000),
                             class_label = "pseudoabs",
                             return_pres = TRUE)

# Alternative method
# set.seed(1)
# spp_all <- sample_background(data = spp_thin, 
#                             raster = env_vars,
#                             n = 10000,
#                             method = 'random',
#                             class_label = "background",
#                             return_pres = TRUE)

# plot output
ggplot() +
  geom_spatraster(data = env_vars[[1]]) +
  scale_fill_cross_blended_c(palette = 'arid', direction = -1) +
  geom_sf(data = spp_all[spp_all$class == 'presence',]) + 
  theme_void() +
  guides(fill="none") +

ggplot() +
  geom_spatraster(data = env_vars[[1]]) +
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
top15_vars <- names(dist_env_vars[1:15])

# view pairs
pairs(env_vars[[top15_vars]], maxcells = 10000)

# identify multi-collinearity in environmental variables
vars_uncor <- filter_collinear(env_vars[[top15_vars]], cutoff = 0.7, method = "cor_caret")
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
# create the model formula
spp_rec <- recipe(spp_all, formula = class ~ .)

# design the model workflows
spp_models <- workflow_set(
  preproc = list(default = spp_rec),
  models = list(maxent = sdm_spec_maxent(tune = 'sdm'),
                rf = sdm_spec_rf(tune = 'sdm'),
                gbm = sdm_spec_boost_tree(tune = 'sdm')),
  cross = TRUE)  %>%
  # tweak controls to store information needed later to create the ensemble
  option_add(control = control_ensemble_grid())

# run the models
set.seed(1)
spp_models <- spp_models %>%
                workflow_map('tune_grid',
                            resamples = spp_cv,
                            grid = 5,
                            metrics = sdm_metric_set(),
                            verbose = T)

# show mean metric values for each model configuration
# spp_models %>% collect_metrics()
autoplot(spp_models)

# extract hyper-parameters and metrics for individual models
model_results <- spp_models$result
names(model_results) <- spp_models$wflow_id
model_parameters <- map_dfr(model_results, ~.x, .id = 'model') %>%
  dplyr::select(model, id, .metrics) %>%
  unnest(.metrics)

#### what are all of the hyper parameters used in the models?
params <- names(model_parameters)[which(!names(model_parameters) %in% c('model', 'id', '.metric', '.estimator', '.estimate', '.config'))]

(hyper_parameters_used <- model_parameters %>%
  group_by(model) %>%
  summarise(across(params, ~ toString(sort(unique(.))))))

# create the model ensemble by selecting the best model using boyce_cont
spp_ensemble <- simple_ensemble() %>%
  add_member(spp_models, metric = "boyce_cont")
autoplot(spp_ensemble)

# what are the best hyper-parameters for these models?
final_parameters <- lapply(spp_ensemble$workflow, extract_spec_parsnip)
species_final_params <- data.frame(spp = species_name, 
           maxent_fc = final_parameters[[1]]$args$feature_classes[[2]],
           maxent_reg = final_parameters[[1]]$args$regularization_multiplier[[2]],
           rf_mtry = final_parameters[[2]]$args$mtry[[2]],
           gbm_mtry = final_parameters[[3]]$args$mtry[[2]],
           gbm_trees = final_parameters[[3]]$args$trees[[2]],
           gbm_tdepth = final_parameters[[3]]$args$tree_depth[[2]],
           gbm_lrate = final_parameters[[3]]$args$learn_rate[[2]],
           gbm_lossr = final_parameters[[3]]$args$loss_reduction[[2]],
           gbm_stop = final_parameters[[3]]$args$stop_iter[[2]])
species_final_params

### 6. Predict species suitability ----
# use all models or a metric threshold e.g. boyce_cont >= 0.25
pred_prob <- predict_raster(
                spp_ensemble, 
                env_vars,
                metric_thresh = c("boyce_cont", 0.25),
                fun = "median")

# view prediction
ggplot() +
  geom_spatraster(data = pred_prob, 
                  aes(fill = median)) +
  scale_fill_whitebox_c(palette = 'high_relief', direction = -1,
                        name = 'Habitat\nsuitability',
                        limits = c(0,1)) +
  geom_sf(data = spp_all[spp_all$class == 'presence',], 
          size = 0.5, alpha = 0.25) +
  theme_void()

# apply the P10 threshold
pred_binary <- sdm_threshold(pred_prob, st_coordinates(spp_thin), binary = TRUE)

# select threshold value
# spp_ensemble <- calib_class_thresh(
#                     spp_ensemble,
#                     class_thresh = c('sensitivity', 0.75), 
#                     metric_thresh = c("boyce_cont", 0.4))

# # create a binary map
# pred_binary <- predict_raster(
#   spp_ensemble, 
#   env_vars,
#   type = "class",
#   class_thresh = c('sensitivity', 0.75), 
#   metric_thresh = c("boyce_cont", 0.4),
#   fun = "median")

# view prediction
ggplot() +
  geom_spatraster(data = as.factor(pred_binary)) +
  scale_fill_manual(values = c('gray90', 'orange'), 
                    na.value = 'transparent',
                    name = 'Habitat\nsuitability',
                    labels = c('absent', 'present', '')) +
  geom_sf(data = spp_all[spp_all$class == 'presence',], 
          size = 0.5, alpha = 0.25, col = 'black') +
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
  labs(x = '', y = 'Mean dropout loss (1 - AUC)') +
  theme_bw()

# top 4 important variables
vi_top4 <- vi %>%
  filter(variable %in% names(env_vars)) %>%
  as.data.frame() %>%
  group_by(variable) %>%
  summarise(mean_dropout_loss = mean(dropout_loss)) %>%
  dplyr::arrange(-mean_dropout_loss) %>%
  top_n(4)

### 8. Partial dependence plots ----
pdp <- model_profile(explain_models, variables = vi_top4$variable)
agg_data <- as.data.frame(pdp$agr_profiles)

agg_data %>%
  filter(`_vname_` == vi_top4$variable[1]) %>%
  ggplot(aes(x = `_x_`, y = `_yhat_`)) + 
  geom_line(lwd = 1, col = 'lightblue') +
  labs(x = vi_top4$variable[1], y = 'mean prediction') +
  scale_y_continuous(limits = c(0, max(agg_data$`_yhat_`))) +
  theme_bw() +

agg_data %>%
  filter(`_vname_` == vi_top4$variable[2]) %>%
  ggplot(aes(x = `_x_`, y = `_yhat_`)) + 
  geom_line(lwd = 1, col = 'lightblue') +
  labs(x = vi_top4$variable[2], y = 'mean prediction') +
  scale_y_continuous(limits = c(0, max(agg_data$`_yhat_`))) +
  theme_bw() +

agg_data %>%
  filter(`_vname_` == vi_top4$variable[3]) %>%
  ggplot(aes(x = `_x_`, y = `_yhat_`)) + 
  geom_line(lwd = 1, col = 'lightblue') +
  labs(x = vi_top4$variable[3], y = 'mean prediction') +
  scale_y_continuous(limits = c(0, max(agg_data$`_yhat_`))) +
  theme_bw() +

agg_data %>%
  filter(`_vname_` == vi_top4$variable[4]) %>%
  ggplot(aes(x = `_x_`, y = `_yhat_`)) + 
  geom_line(lwd = 1, col = 'lightblue') +
  labs(x = vi_top4$variable[4], y = 'mean prediction') +
  scale_y_continuous(limits = c(0, max(agg_data$`_yhat_`))) +
  theme_bw() +
  plot_annotation(tag_levels = 'a', tag_suffix = ')')

spp_all %>%
  dplyr::select(class, vi_top4$variable) %>%
  plot_pres_vs_bg(class)

#### Identify overlap with Protected Areas ----
# load in protected area shapefiles. There are 3 files, so we want to load them all in together and then bind them into one file
prot_areas <- list.files('data/WDPA_Madagascar', pattern = '*.shp', full.names = TRUE)
prot_areas_list <- lapply(prot_areas, read_sf)
# bind the 3 files togther
prot_areas_all <- bind_rows(prot_areas_list) %>% filter(MARINE == 0)

#### convert to equal area projection
# convert the protected areas
prot_areas_all %>% 
  st_transform(crs = 'EPSG:29702') -> prot_areas_all_proj

# convert the presence/absence raster
pred_binary %>% 
  project(.,vect(prot_areas_all_proj), method = 'near') -> pred_binary_proj

# visualise the different projections
par(mfrow=c(1,2))
plot(pred_binary)
plot(vect(prot_areas_all), add = TRUE)
plot(pred_binary_proj)
plot(vect(prot_areas_all_proj), add = TRUE)

# What is the area of species presences?
# we select and sum only the cells with 1's, then multiply this by the size of the raster cells and lastly divide this by meters to get a result in km2.
pres_area <- (sum(pred_binary_proj[] == 1, na.rm = TRUE) * (res(pred_binary_proj)[1]*res(pred_binary_proj)[2]) / (1000^2))
paste('The area of species presences is',pres_area, 'km2')

# Calculate the area of all cells
all_area <- (sum(!is.na(pred_binary_proj[])) * (res(pred_binary_proj)[1]*res(pred_binary_proj)[2]) / (1000^2))
paste('The area of all cells is',all_area, 'km2')

# And lastly calculate the percentage of coverage of our species across all of Madagascar
paste('The species presences cover',round(pres_area/all_area*100, 2), '% of Madagascar')

#### We now want to work out what % of our species is found within Protected Areas

# create custom function to calculate the proportion of area covered by each Protected Area
sum_cover <- function(x){
  list(x %>%
         group_by(value) %>%
         summarize(total_area = sum(coverage_area)))
}

# extract the amount of area covered 
extract_all <- exact_extract(pred_binary_proj, prot_areas_all_proj, coverage_area = TRUE, summarize_df = TRUE, fun = sum_cover)
# add the names of the protected areas back on to our extraction
names(extract_all) <- prot_areas_all_proj$ORIG_NAME

# convert the list to a data frame
extract_df <- bind_rows(extract_all, .id = 'ORIG_NAME')
# take a look at the first 6 rows
head(extract_df)

# we can now sum all of the area that overlaps with the protected areas for presences (i.e. 1's) and divide this by the total area of all presences
area_under_pas <- extract_df %>% 
  filter(value == 1) %>% 
  summarise(sum(total_area)/(1000^2))

paste(round(area_under_pas/pres_area * 100, 2),'% of the predicted presences are found within protected areas')

# Our final step is to join our IUCN protected area categories onto our presence area data.frame. This will provide us with some information on what percentage of our species area is conserved under different categories. This provides important context on both the quality and quantity of protected areas overlapping with our species range:

iucn_cat <- prot_areas_all_proj %>% 
  st_drop_geometry() %>% 
  dplyr::select(ORIG_NAME, IUCN_CAT)

extract_df %>% 
  left_join(iucn_cat, by = 'ORIG_NAME', relationship = 'many-to-many') %>% 
  filter(value == 1) %>%
  group_by(IUCN_CAT) %>%
  summarise(area = sum(total_area)/(1000^2)) %>%
  mutate(perc = round(area/sum(area) * 100, 2))

#### END ####
renv::snapshot()
