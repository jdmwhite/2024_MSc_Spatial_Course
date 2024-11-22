# RGB Kew MSc Spatial Analysis Course 2024
# Joseph White, Carolina Tovar, Moabe Fernandes, Felix Lim, Jenny Williams
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
library(mapview)
library(blockCV)
source('scripts/helper_functions/mtp.R')

### 1. Load in data ----
#### 1a. Country boundary ----
mad <- vect(ne_countries(country = 'Madagascar', scale = 'medium'))
crs(mad) <- "lonlat"
# we'll use this as our area of interest, but consider changing this to more closely reflect your species geographical area

#### 1b. occurrence data ----
# Name your species. *CHANGE THIS TO YOUR ASSIGNED SPECIES*
species_name <- c("Bauhinia madagascariensis")

# download GBIF occurrence data for this species; this takes time if there are many data points!
gbif_download <- occ_search(scientificName = species_name, hasCoordinate = TRUE, country = 'MG', limit = 10000)

# select only the data (ignore the other metadata for now)
gbif_data <- gbif_download$data
head(gbif_data)

#### 1c. environmental data ----
# load in the WorldClim 2.1 data downloaded from: https://www.worldclim.org/data/worldclim21.html
rast_files <- list.files(paste0(getwd(),'/data/environmental_data/'), full.names = T, recursive = T)
env_vars <- app(sds(rast_files),c)
# give them better names
names(env_vars) <- c("mean_ann_t",'mean_t_warm_q','mean_t_cold_q','ann_p', 'p_wet_m','p_dry_m','p_seas','p_wet_q','p_dry_q','p_warm_q','p_cold_q',"mean_diurnal_t_range","isothermality", "t_seas", 'max_t_warm_m','min_t_cold_m',"t_ann_range",'mean_t_wet_q','mean_t_dry_q', 'elev')
names(env_vars)
env_vars
plot(env_vars[[1]], main = names(env_vars)[[1]])

### 2. Clean occurrence data ----
# keep the necessary columns and rename them
spp_raw <- gbif_data %>% 
              dplyr::select(gbifID, scientificName, latitude = decimalLatitude, longitude = decimalLongitude)

# CoordinateCleaner: remove occurrences with irregular values
# GBIF recommended workflow: https://www.gbif.org/es/data-use/5lkEzTJaUPAZyCGUC3PDKC/coordinatecleaner-fast-and-standardized-cleaning-of-species-occurrence-data
# See also: https://www.nature.com/articles/s41598-024-56158-3
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

# visualise raw vs. clean data
ggplot() +
  geom_spatvector(data = mad) +
  geom_point(data = spp_raw, aes(x = longitude, y = latitude)) +
  labs(title = 'Raw GBIF data') +
  theme_bw() +
  
ggplot() +
  geom_spatvector(data = mad) +
  geom_point(data = spp_clean, aes(x = longitude, y = latitude), col = 'forestgreen') +
  labs(title = 'Cleaned GBIF data') +
  theme_bw()

# Remove records with redundant information and to reduce errors associated with collection bias
# spatial thin
set.seed(1234567)
spp_thin <- thin_by_cell(spp_clean, raster = env_vars[[1]])
nrow(spp_thin)

# convert to sf 
spp_thin <- st_as_sf(spp_thin, coords = c('longitude', 'latitude'), crs = 'EPSG:4326')

# visualise output
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

# Alternative method: background sampling
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
  geom_sf(data = spp_all[spp_all$class == 'presence',], col = '#F8766D') + 
  theme_bw() +
  guides(fill="none") +

ggplot() +
  geom_spatraster(data = env_vars[[1]]) +
  scale_fill_cross_blended_c(palette = 'arid', direction = -1) +
  geom_sf(data = spp_all, aes(col = class)) + 
  theme_bw() +
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
# for effective evaluation and to test model transferability, we set up a several training vs. testing folds using cross-validation

n_blocks <- 5 # if you get warnings that there are no presences in some groups when fitting your models, try using less blocks
set.seed(100)
# cluster points using kmeans into n blocks
spp_cv <- spatial_clustering_cv(data = spp_all, 
                                cluster_function = 'kmeans',
                                v = n_blocks)

# # alternative 1
# n_blocks <- 40
# spp_cv <- spatial_block_cv(data = spp_all, method = 'random', v = n_blocks)

# # alternative 2
# n_blocks <- 5
# cv <- blockCV::cv_spatial(spp_all, column = 'class', r = env_vars, k = n_blocks, selection = 'systematic', plot = F, seed = 1, rows_cols = c(160))
# spp_cv <- blockcv2rsample(cv, spp_all)

# # alternative 3: look up 'targeted group backgrounds' on tidysdm
tidysdm::check_splits_balance(spp_cv, 'class')

# view folds
autoplot(spp_cv, aes(shape = class)) +
  geom_spatvector(data = mad, fill = NA) +
  theme_void()

### 5. Run distribution models ----
# create the model formula
spp_rec <- recipe(spp_all, formula = class ~ .) # the . = all other columns in the dataset, but ignores .geometry

# design the model workflows
spp_models <- workflow_set(
  preproc = list(default = spp_rec),
  models = list(maxent = sdm_spec_maxent(tune = 'sdm'),
                rf = sdm_spec_rf(tune = 'sdm'),
                gbm = sdm_spec_boost_tree(tune = 'sdm')),
  cross = TRUE)  %>%
  # tweak controls to store information needed later to create the ensemble
  option_add(control = control_ensemble_grid())

# alternative: test using more or less models, see https://evolecolgroup.github.io/tidysdm/reference/index.html for list of alternative models

# run the models
# here we specify our combination of splits/folds/resamples, how many parameter combinations we test per algorithm (grid = 5), and which parameters we want to test (sdm_metric_set).
set.seed(1)
spp_models <- spp_models %>%
                workflow_map('tune_grid',
                            resamples = spp_cv,
                            grid = 10,
                            metrics = sdm_metric_set(),
                            verbose = T)

# show mean metric values for each model configuration
spp_models %>% collect_metrics()
autoplot(spp_models)
# We are testing all of these different parameters to optimise our models, based on the available data (occs + env), algorithms (maxent, gbm, rf) and parameters we choose (CV, pseudoabs, hyperparameters), because we want to create an accurate prediction; not hypothesis test.

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

# save hyper-parameters used
write_csv(hyper_parameters_used, 'output/tables/hyper_paramaters_used.csv')

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

# Save parameter outputs
write_csv(species_final_params, 'output/tables/final_parameters.csv')

# Extract the evaluation metrics
(eval_metrics <- spp_ensemble %>% collect_metrics())

# save the evaluation metrics
write_csv(eval_metrics, 'output/tables/evaluation_metrics.csv')

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
  geom_spatvector(data = mad, fill = 'transparent') +
  theme_bw()

# using a species presence specific threshold
# apply the 10th percentile training presence threshold
pred_binary <- sdm_threshold(pred_prob, st_coordinates(spp_thin), type = 'percentile', threshold = 0.1, binary = TRUE)
# this omits the bottom 10th percentile of presences from the binary map
# test different threshold values (e.g. 0, 0.2, 0.5)

# # alternative: use model-specific thresholds
# # select threshold value
# spp_ensemble <- calib_class_thresh(
#                     spp_ensemble,
#                     class_thresh = c('tss_max'),
#                     metric_thresh = c("boyce_cont", 0.25))
# 
# # # create a binary map
# pred_binary <- predict_raster(
#   spp_ensemble,
#   env_vars,
#   type = "class",
#   class_thresh = c('tss_max'),
#   metric_thresh = c("boyce_cont", 0.25),
#   fun = "median")

# view prediction
ggplot() +
  geom_spatraster(data = pred_prob, 
                  aes(fill = median)) +
  scale_fill_whitebox_c(palette = 'high_relief', direction = -1,
                        name = 'Habitat\nsuitability',
                        limits = c(0,1)) +
  geom_sf(data = spp_all[spp_all$class == 'presence',], 
          size = 0.5, alpha = 0.25) +
  geom_spatvector(data = mad, fill = 'transparent') +
  theme_bw() +

# note: if you use tss_max, you need to change the order of the colour values and the order of the labels!
ggplot() +
  geom_spatraster(data = as.factor(pred_binary)) +
  scale_fill_manual(values = c('gray90', 'orange'),
                    na.value = 'transparent',
                    name = 'Binary distribution',
                    labels = c('absent', 'present', '')) +
  geom_sf(data = spp_all[spp_all$class == 'presence',],
          size = 0.5, alpha = 0.25, col = 'black') +
  geom_spatvector(data = mad, fill = 'transparent') +
  theme_bw() &
  plot_annotation(tag_levels = 'A', tag_suffix = ')')

ggsave('output/figures/1_map_predictions.png',
       width = 8.23, height = 4.9, dpi = 320)

### 7. Variable importance ----
explain_models <- explain_tidysdm(spp_ensemble)
vi <- model_parts(explain_models, type = 'variable_importance')

# view variable importance
vi %>%
  filter(variable %in% names(env_vars)) %>%
  ggplot() +
  geom_boxplot(aes(x = reorder(variable, dropout_loss), y = dropout_loss),
               fill = 'lightblue') +
  coord_flip() +
  labs(x = '', y = 'Mean dropout loss (1 - AUC)') +
  theme_bw()
ggsave('output/figures/2_variable_importance.png',
       width = 5.73, height = 3, dpi = 320)

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
  plot_annotation(tag_levels = 'A', tag_suffix = ')')

ggsave('output/figures/3_response_curves.png',
       width = 7.84, height = 5.2, dpi = 320)

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

# What is the area of predicted species presences?
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
  summarise(total_area_under_pas = sum(total_area)/(1000^2))

paste(round(area_under_pas/pres_area * 100, 2),'% of the predicted presences are found within protected areas')

#### join the protected areas vect file onto the area calculations
prot_areas_all_proj_area <- prot_areas_all_proj %>%
  left_join(filter(extract_df, value == 1), by = 'ORIG_NAME')

ggplot() +
  geom_spatvector(data = mad, fill = 'white') +
  geom_spatvector(data = prot_areas_all_proj_area, aes(fill = total_area)) +
  scale_fill_continuous(name = 'Total area of\nspecies covered\n(km2)', type = 'viridis') +
  theme_bw()

mapview::mapview(st_as_sf(prot_areas_all_proj_area), zcol = 'total_area')

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

#### For more information on each function and alternative methods using tidysdm, see:
# https://evolecolgroup.github.io/tidysdm/articles/a0_tidysdm_overview.html
# https://evolecolgroup.github.io/tidysdm/articles/a2_tidymodels_additions.html
