
#################################################################################################################
# 1. Sourcing required packages -------------------------------------------
#################################################################################################################
#################################################################################################################
packages_required <- c("tidyverse", "h2o", "geodata", "terra")

# check and install packages that are not yet installed
installed_packages <- packages_required %in% rownames(installed.packages())
if(any(installed_packages == FALSE)){
  install.packages(packages_required[!installed_packages])}

# load required packages
suppressWarnings(suppressPackageStartupMessages(invisible(lapply(packages_required, library, character.only = TRUE))))

rm(list = ls())

#################################################################################################################
# 2. train the models  
#################################################################################################################
#################################################################################################################
pathIn <- pathOut <- "/Volumes/Data/01UM6P/AI_Papers/proposal/Proposal_versions/Chapters/Chapter1/data_extraction_md/script/downstream_fert_recom"
extracted <- read.csv(paste(pathIn, "extracted_final.csv", sep = "/"), header = T, sep = ",")
gold <- read.csv(paste(pathIn, "gold_final.csv", sep = "/"), header = T, sep = ",")
colnames(extracted)
glimpse(extracted)
colnames(gold)
glimpse(gold)

#select the important attributes
gold <- gold |> dplyr::select(dataset, paper_key, file_name, region, zone, woreda, lat_decimal, long_decimal, soil_type, variety, n_rate, p_rate, grain_yield)
extracted <- extracted |> dplyr::select(dataset, paper_key, file_name, region, zone, woreda, lat_decimal, long_decimal, soil_type, variety, n_rate, p_rate, grain_yield)
colnames(extracted)
glimpse(extracted)
colnames(gold)
glimpse(gold)

# extract the topography, climate and soil
elevation <- geodata::elevation_30s(country = "Ethiopia", path = tempdir())
names(elevation) <- "elevation"
plot(elevation)

aspect <-  elevation |> terra::terrain(v = "aspect")
slope <-  elevation |> terra::terrain(v = "slope")
TPI <-  elevation |> terra::terrain(v = "TPI")
TRI <-  elevation |> terra::terrain(v = "TRI")

topo <- c(elevation, aspect, slope, TPI, TRI)
plot(topo)

#soil
soil_vars <- c("clay", "C.tot","db.od", "eCEC.f", "N.tot", "OC", "P", "pH.H2O", "sand", "silt")
# Create an empty list
soil_list <- list()

# Download each soil variable
for (v in soil_vars) {
  
  cat("Downloading:", v, "\n")
  
  soil_list[[v]] <- geodata::soil_af_isda(
    var = v,
    path = tempdir()
  )
}

# Combine into a terra SpatRaster
soil_stack <- do.call(c, unname(soil_list))

# Check result
soil_stack
names(soil_stack)
names(soil_stack) <- c("clay", "c_tot","dbod", "ecec", "n_tot", "oc", "p_soil", "ph", "sand", "silt")

####climate

# Variables to download
clim_vars <- c("tmin", "tmax", "prec")

# Create an empty list
clim_list <- list()

# Download and process climate variables
for (v in clim_vars) {
  
  cat("Processing:", v, "\n")
  
  # Download WorldClim monthly climate data
  r <- geodata::worldclim_country(country = "Eth",
    var  = v,
    path = tempdir()
  )
  
  # Select June, July, August
  jja <- r[[6:8]]
  
  # Calculate seasonal climate statistics
  if (v == "prec") {
    
    # Total JJA precipitation (mm)
    jja_result <- terra::app(
      jja,
      fun = sum,
      na.rm = FALSE
    )
    
  } else {
    
    # Mean JJA temperature (degrees C)
    jja_result <- terra::app(
      jja,
      fun = mean,
      na.rm = FALSE
    )
  }
  
  # Rename output layer
  names(jja_result) <- paste0(v, "_JJA")
  
  # Store raster in list
  clim_list[[v]] <- jja_result
}

#Combine all three climate rasters
clim_stack <- do.call(c, unname(clim_list))
names(clim_stack)
plot(clim_stack)

#Extracting the climate, soil and topography covariates
# ---------------------------------------------------------------------------------------------
# Step 1. Combine the field data and prepare the point coordinates
# ---------------------------------------------------------------------------------------------
field_data <- rbind(gold, extracted)

# empty text cells are read as "" by read.csv: set them to NA
for (v in c("region", "zone", "woreda")) {
  field_data[[v]][field_data[[v]] == ""] <- NA
}

# recorded coordinates (a record needs both lat and long)
field_data$lon <- field_data$long_decimal
field_data$lat <- field_data$lat_decimal
no_xy <- is.na(field_data$lon) | is.na(field_data$lat)
field_data$lon[no_xy] <- NA
field_data$lat[no_xy] <- NA
field_data$xy_source <- ifelse(no_xy, NA, "recorded")

table(field_data$dataset, no_xy)   # records without coordinates

# ---------------------------------------------------------------------------------------------
# Step 2. Admin boundaries and their centroids (woreda, zone, region)
# The same boundary files were used to fill region/zone/woreda in the EDA notebook,
# so the names (Sidama, Central, South and South West Ethiopia) match the csv.
# ---------------------------------------------------------------------------------------------
adm3 <- terra::vect(file.path(pathIn, "ethiopia_admin", "eth_adm3_regions.geojson"))
adm3 <- adm3[, c("region", "zone", "woreda")]

adm2 <- terra::aggregate(adm3[, c("region", "zone")], by = c("region", "zone"))   # zone polygons
adm1 <- terra::aggregate(adm3[, "region"], by = "region")                          # region polygons

# centroids are forced to fall inside the polygon (inside = TRUE)
cent3 <- data.frame(region = adm3$region, zone = adm3$zone, woreda = adm3$woreda,
                    terra::crds(terra::centroids(adm3, inside = TRUE)))
cent2 <- data.frame(region = adm2$region, zone = adm2$zone,
                    terra::crds(terra::centroids(adm2, inside = TRUE)))
cent1 <- data.frame(region = adm1$region,
                    terra::crds(terra::centroids(adm1, inside = TRUE)))
head(cent3)

# two small helpers for the name matching
clean_name <- function(x) {
  x <- tolower(x)
  x <- gsub("[^a-z0-9 ]", "", x)
  x <- gsub("\\b(woreda|zone|district|town|city)\\b", "", x, perl = TRUE)
  gsub("\\s+", "", x)
}

# index of the most similar name in 'ref' (NA if the similarity is below min_sim)
best_match <- function(x, ref, min_sim = 0.85) {
  if (length(ref) == 0) return(NA)
  sim <- 1 - drop(adist(x, ref)) / pmax(nchar(x), nchar(ref))
  if (max(sim) >= min_sim) which.max(sim) else NA
}

# ---------------------------------------------------------------------------------------------
# Step 3. Records without lat/long: use the centroid of the woreda, then the zone, then the region
# ---------------------------------------------------------------------------------------------
for (i in which(no_xy)) {

  reg <- field_data$region[i]
  zon <- field_data$zone[i]
  wor <- field_data$woreda[i]

  # search only inside the reported region (whole Ethiopia if the region is empty or unknown)
  ref3 <- cent3
  ref2 <- cent2
  if (!is.na(reg) && reg %in% cent1$region) {
    ref3 <- cent3[cent3$region == reg, ]
    ref2 <- cent2[cent2$region == reg, ]
  }

  # 1) woreda / district centroid
  if (!is.na(wor)) {
    j <- best_match(clean_name(wor), clean_name(ref3$woreda))
    if (!is.na(j)) {
      field_data$lon[i] <- ref3$x[j]
      field_data$lat[i] <- ref3$y[j]
      field_data$xy_source[i] <- "woreda_centroid"
      next
    }
  }

  # 2) zone centroid
  if (!is.na(zon)) {
    j <- best_match(clean_name(zon), clean_name(ref2$zone))
    if (!is.na(j)) {
      field_data$lon[i] <- ref2$x[j]
      field_data$lat[i] <- ref2$y[j]
      field_data$xy_source[i] <- "zone_centroid"
      next
    }
  }

  # 3) region centroid
  if (!is.na(reg)) {
    j <- match(reg, cent1$region)
    if (!is.na(j)) {
      field_data$lon[i] <- cent1$x[j]
      field_data$lat[i] <- cent1$y[j]
      field_data$xy_source[i] <- "region_centroid"
    }
  }
}

table(field_data$dataset, field_data$xy_source, useNA = "ifany")

# ---------------------------------------------------------------------------------------------
# Step 4. Extract the topography, soil and climate values at the points
# ---------------------------------------------------------------------------------------------
has_xy <- !is.na(field_data$lon) & !is.na(field_data$lat)
cat("Records without any location (dropped):", sum(!has_xy), "\n")

pts <- terra::vect(field_data[has_xy, ], geom = c("lon", "lat"), crs = "EPSG:4326", keepgeom = TRUE)

topo_vals <- terra::extract(topo,       pts, ID = FALSE)
soil_vals <- terra::extract(soil_stack, pts, ID = FALSE)
clim_vals <- terra::extract(clim_stack, pts, ID = FALSE)

# field data + covariates (one row per record)
field_cov <- cbind(as.data.frame(pts), topo_vals, soil_vals, clim_vals)

colnames(field_cov)
dim(field_cov)
field_cov |> summarise_all(~sum(is.na(.)))   # NA values per column
field_cov_final <- field_cov |> dplyr::select(dataset, soil_type, variety, elevation, aspect, slope, 
                                 TPI, TRI, clay, c_tot, dbod, ecec, n_tot, oc, 
                                 p_soil, ph, sand, silt, tmin_JJA, tmax_JJA, 
                                 prec_JJA, n_rate, p_rate, grain_yield) |>
                                na.omit() |> unique()
field_cov_final |> filter(dataset == "gold") |> summarise_all(~sum(is.na(.)))
field_cov_final |> filter(dataset == "extracted") |> summarise_all(~sum(is.na(.)))
write.csv(field_cov_final, paste(pathOut, "field_cov_final.csv", sep = "/"), row.names = F, col.names = T)

#field_cov_final$variety <- as.factor(field_cov_final$variety)
#field_cov_final$variety <- as.factor(field_cov_final$variety)

# ---------------------------------------------------------------------------------------------
# XGBoost with grid search: gold vs extracted
# Each dataset is split 80/20; the grid is tuned with 5-fold CV on the training part
# (number of trees by early stopping) and the best model is evaluated on the test part.
# XGBoost handles missing values (NA) itself; text columns are coded as integers.
# ---------------------------------------------------------------------------------------------
rm(list = ls())
pathOut <- "/Volumes/Data/01UM6P/AI_Papers/proposal/Proposal_versions/Chapters/Chapter1/data_extraction_md/script/downstream_fert_recom"
field_cov_final <- read.csv(paste(pathOut, "field_cov_final.csv", sep = "/"), header = T, sep = ",") |>
  dplyr::filter(grain_yield <= 7500)
dim(field_cov_final)
colnames(field_cov_final) <- c("dataset", "soil_type", "variety", "elevation", "aspect", "slope",
                                 "TPI", "TRI", "clay", "c_tot", "dbod", "ecec",
                                 "n_tot", "oc", "p_soil", "ph", "sand", "silt",
                                 "tmin", "tmax", "prec",
                                 "n_rate", "p_rate", "grain_yield")
field_cov_final <- field_cov_final |> dplyr::select(-c("soil_type", "variety"))
colnames(field_cov_final)
dim(field_cov_final |> dplyr::filter(dataset == "extracted"))
dim(field_cov_final |> dplyr::filter(dataset == "gold"))
#field_cov_final$variety <- as.factor(field_cov_final$variety)
#field_cov_final$variety <- as.factor(field_cov_final$variety)
if (!requireNamespace("xgboost", quietly = TRUE)) install.packages("xgboost")

xgb_response   <- "grain_yield"
xgb_predictors <- setdiff(names(field_cov_final), c("dataset", "grain_yield"))
xgb_colors     <- c(extracted = "#dca237", gold = "#469c76")

# hyper-parameter grid (2 x 3 x 2 x 2 x 2 = 48 combinations per dataset)
xgb_grid <- expand.grid(eta = c(0.03, 0.1),
                        max_depth = c(3, 5, 7),
                        subsample = c(0.7, 1),
                        colsample_bytree = c(0.6, 1),
                        min_child_weight = c(1, 5))

# every predictor must be numeric: text columns (soil_type, variety) -> integer codes, NA stays NA
# (the coding is done once, so both datasets use the same codes)
# the table is rebuilt from field_cov (it keeps the dataset column) with the variables of field_cov_final
xgb_all <- field_cov_final |>
  dplyr::select(dplyr::all_of(c("dataset", xgb_predictors, xgb_response)))

xgb_complete_cases <- TRUE    # TRUE = drop rows with any NA and duplicates (as for field_cov_final); FALSE = keep the NA
if (xgb_complete_cases) {
  xgb_all <- xgb_all |> na.omit() |> unique()
}
table(xgb_all$dataset)

for (v in xgb_predictors) {
  if (!is.numeric(xgb_all[[v]])) {
    xgb_all[[v]] <- as.integer(factor(as.character(xgb_all[[v]])))
  }
}
sapply(xgb_all[, xgb_predictors], class)   # all should be numeric or integer

xgb_results <- list()
xgb_varimp  <- list()
xgb_pred    <- list()
xgb_shap    <- list()

for (ds in c("gold", "extracted")) {

  cat("\n=== XGBoost:", ds, "===\n")

  # data of one dataset
  xgb_data <- xgb_all |>
    dplyr::filter(dataset == ds) |>
    dplyr::select(-dataset)

  # train / test split
  set.seed(1234)
  train_id  <- sample(nrow(xgb_data), size = round(0.8 * nrow(xgb_data)))
  xgb_train <- xgb_data[train_id, ]
  xgb_test  <- xgb_data[-train_id, ]

  dtrain <- xgboost::xgb.DMatrix(data = as.matrix(xgb_train[, xgb_predictors]), label = xgb_train[[xgb_response]])
  dtest  <- xgboost::xgb.DMatrix(data = as.matrix(xgb_test[, xgb_predictors]),  label = xgb_test[[xgb_response]])

  # grid search: 5-fold CV for every combination, trees chosen by early stopping
  xgb_cv_table <- xgb_grid
  xgb_cv_table$nrounds <- NA
  xgb_cv_table$cv_rmse <- NA

  set.seed(444)
  for (g in 1:nrow(xgb_grid)) {

    params <- list(objective = "reg:squarederror",
                   eta = xgb_grid$eta[g],
                   max_depth = xgb_grid$max_depth[g],
                   subsample = xgb_grid$subsample[g],
                   colsample_bytree = xgb_grid$colsample_bytree[g],
                   min_child_weight = xgb_grid$min_child_weight[g])

    cv <- xgboost::xgb.cv(params = params, data = dtrain, nrounds = 1000, nfold = 5,
                          early_stopping_rounds = 20, verbose = 0)

    cv_log <- as.data.frame(cv$evaluation_log)
    xgb_cv_table$nrounds[g] <- which.min(cv_log$test_rmse_mean)
    xgb_cv_table$cv_rmse[g] <- min(cv_log$test_rmse_mean)
  }

  # best combination = lowest cross-validated RMSE
  best <- xgb_cv_table[which.min(xgb_cv_table$cv_rmse), ]
  print(best)

  best_params <- list(objective = "reg:squarederror",
                      eta = best$eta,
                      max_depth = best$max_depth,
                      subsample = best$subsample,
                      colsample_bytree = best$colsample_bytree,
                      min_child_weight = best$min_child_weight)

  xgb_fit <- xgboost::xgb.train(params = best_params, data = dtrain, nrounds = best$nrounds)

  # performance on the training and the test data
  obs_train  <- xgb_train[[xgb_response]]
  pred_train <- predict(xgb_fit, dtrain)
  obs        <- xgb_test[[xgb_response]]
  pred       <- predict(xgb_fit, dtest)

  xgb_results[[ds]] <- data.frame(dataset = ds,
                                  n_train = nrow(xgb_train),
                                  n_test  = nrow(xgb_test),
                                  nrounds = best$nrounds,
                                  eta = best$eta,
                                  max_depth = best$max_depth,
                                  subsample = best$subsample,
                                  colsample_bytree = best$colsample_bytree,
                                  min_child_weight = best$min_child_weight,
                                  cv_rmse    = best$cv_rmse,
                                  train_MSE  = mean((obs_train - pred_train)^2),
                                  train_MAE  = mean(abs(obs_train - pred_train)),
                                  train_RMSE = sqrt(mean((obs_train - pred_train)^2)),
                                  train_R2   = 1 - sum((obs_train - pred_train)^2) / sum((obs_train - mean(obs_train))^2),
                                  test_MSE   = mean((obs - pred)^2),
                                  test_MAE   = mean(abs(obs - pred)),
                                  test_RMSE  = sqrt(mean((obs - pred)^2)),
                                  test_R2    = 1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2))

  # variable importance (gain, as a share of the total; 0 for variables never used)
  imp <- as.data.frame(xgboost::xgb.importance(model = xgb_fit))
  vi  <- data.frame(variable = xgb_predictors,
                    percentage = imp$Gain[match(xgb_predictors, imp$Feature)],
                    dataset = ds)
  vi$percentage[is.na(vi$percentage)] <- 0
  xgb_varimp[[ds]] <- vi

  # measured vs predicted yield on the test data
  xgb_pred[[ds]] <- data.frame(dataset = ds, observed = obs, predicted = pred)

  # SHAP values of the test data (one value per record and variable) and the scaled variable values (0-1, for the colours)
  shap_values <- predict(xgb_fit, dtest, predcontrib = TRUE)
  feature_values <- as.data.frame(as.matrix(xgb_test[, xgb_predictors]))
  feature_scaled <- as.data.frame(sapply(feature_values, function(x) (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))))

  xgb_shap[[ds]] <- data.frame(dataset  = ds,
                               variable = rep(xgb_predictors, each = nrow(xgb_test)),
                               shap     = as.vector(shap_values[, xgb_predictors]),
                               value    = unlist(feature_scaled, use.names = FALSE))
}

# ---------------------------------------------------------------------------------------------
# Compare the two datasets
# ---------------------------------------------------------------------------------------------
xgb_comparison <- do.call(rbind, xgb_results)
metric_cols <- c("cv_rmse", "train_MSE", "train_MAE", "train_RMSE", "train_R2",
                 "test_MSE", "test_MAE", "test_RMSE", "test_R2")
xgb_comparison[, metric_cols] <- round(xgb_comparison[, metric_cols], 3)
xgb_comparison
write.csv(xgb_comparison, file.path(pathOut, "xgb_gold_vs_extracted_metrics.csv"), row.names = FALSE)

# legend / panel labels of the two datasets
dataset_labels <- c(gold = "Human-annotated", extracted = "Model-extracted")

# theme for the figures (publication size: 7 inch wide, 9 pt axis titles, 8 pt axis and legend text)
theme_pub <- theme_minimal(base_size = 9) +
  theme(plot.title       = element_text(size = 10, face = "bold"),
        axis.title       = element_text(size = 9),
        axis.text        = element_text(size = 8),
        strip.text       = element_text(size = 9),
        legend.title     = element_text(size = 8),
        legend.text      = element_text(size = 8),
        legend.key.size  = unit(0.4, "cm"),
        panel.grid.minor = element_blank())

# MAE, RMSE and R2 for the training and the test data
xgb_metrics_long <- xgb_comparison |>
  dplyr::select(dataset, train_MAE, train_RMSE, train_R2, test_MAE, test_RMSE, test_R2) |>
  tidyr::pivot_longer(cols = -dataset, names_to = c("sample", "metric"), names_sep = "_", values_to = "value")
xgb_metrics_long$sample <- factor(xgb_metrics_long$sample, levels = c("train", "test"))

# bar colour = dataset (dark = test, light = training); the legend uses the same colours
xgb_metrics_long$group <- paste(xgb_metrics_long$dataset, xgb_metrics_long$sample)
metric_colors <- c("gold train"      = "#a3d1bd", "gold test"      = "#469c76",
                   "extracted train" = "#eed0a0", "extracted test" = "#dca237")

p_metrics <- ggplot(xgb_metrics_long, aes(dataset, value, fill = group)) +
  geom_col(aes(group = sample), position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = round(value, 2), group = sample), position = position_dodge(width = 0.8), vjust = -0.3, size = 2.5) +
  facet_wrap(~metric, scales = "free_y") +
  scale_x_discrete(labels = dataset_labels) +
  scale_fill_manual(values = metric_colors, breaks = names(metric_colors),
                    labels = c("Human-annotated (training)", "Human-annotated (test)",
                               "Model-extracted (training)", "Model-extracted (test)")) +
  labs(x = NULL, y = NULL, fill = NULL, title = "XGBoost: training and test performance") +
  theme_pub +
  theme(legend.position = "bottom")
print(p_metrics)
ggsave(file.path(pathOut, "xgb_metrics.png"), p_metrics, width = 7, height = 3.5, dpi = 300)

# measured vs predicted yield (panels without titles, the colour legend identifies the dataset)
xgb_pred_all <- do.call(rbind, xgb_pred)
p_scatter <- ggplot(xgb_pred_all, aes(observed, predicted, colour = dataset)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, colour = "black") +
  facet_wrap(~dataset) +
  scale_colour_manual(values = xgb_colors, labels = dataset_labels) +
  labs(x = "Measured yield (kg/ha)", y = "Predicted yield (kg/ha)", colour = NULL) +
  theme_pub +
  theme(strip.text = element_blank(), legend.position = "bottom")
print(p_scatter)
ggsave(file.path(pathOut, "xgb_pred_vs_obs.png"), p_scatter, width = 7, height = 3.5, dpi = 300)

# variable importance: both datasets side by side, most important variable at the top of each panel
xgb_vi_all <- do.call(rbind, xgb_varimp)
xgb_vi_all <- xgb_vi_all[order(xgb_vi_all$dataset, xgb_vi_all$percentage), ]    # increasing within each dataset
xgb_vi_all$var_ds <- factor(paste(xgb_vi_all$variable, xgb_vi_all$dataset, sep = "___"),
                            levels = paste(xgb_vi_all$variable, xgb_vi_all$dataset, sep = "___"))

p_importance <- ggplot(xgb_vi_all, aes(percentage, var_ds, fill = dataset)) +
  geom_col(width = 0.7, orientation = "y") +
  geom_text(aes(label = paste0(round(100 * percentage, 1), "%")), hjust = -0.1, size = 2.5) +
  facet_wrap(~dataset, scales = "free_y") +
  scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +       # remove the dataset suffix from the labels
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  scale_fill_manual(values = xgb_colors, labels = dataset_labels) +
  labs(x = "Relative importance (gain, share of the total)", y = NULL, fill = NULL,
       title = "XGBoost: variable importance") +
  theme_pub +
  theme(strip.text = element_blank(), legend.position = "bottom")
print(p_importance)
ggsave(file.path(pathOut, "xgb_variable_importance.png"), p_importance, width = 7, height = 5, dpi = 300)

# SHAP summary plot of the test data, gold and extracted in one plot (top 15 variables of each dataset)
# Note: h2o.shap_summary_plot() only works on H2O models; the XGBoost models above are fitted with the
# xgboost package, so the SHAP values come from predict(..., predcontrib = TRUE) and the same kind of plot is drawn with ggplot.
xgb_shap_all <- do.call(rbind, xgb_shap)

shap_rank <- xgb_shap_all |>
  dplyr::group_by(dataset, variable) |>
  dplyr::summarise(mean_abs_shap = mean(abs(shap)), .groups = "drop") |>
  dplyr::group_by(dataset) |>
  dplyr::slice_max(mean_abs_shap, n = 15) |>
  dplyr::ungroup() |>
  dplyr::arrange(dataset, mean_abs_shap)                                        # increasing: the largest ends at the top

shap_plot_data <- xgb_shap_all |>
  dplyr::inner_join(shap_rank, by = c("dataset", "variable"))
shap_plot_data$var_ds <- factor(paste(shap_plot_data$variable, shap_plot_data$dataset, sep = "___"),
                                levels = paste(shap_rank$variable, shap_rank$dataset, sep = "___"))

p_shap <- ggplot(shap_plot_data, aes(shap, var_ds, colour = value)) +
  geom_vline(xintercept = 0, colour = "grey50") +
  geom_jitter(height = 0.2, size = 1, alpha = 0.7) +
  facet_wrap(~dataset, scales = "free_y", labeller = labeller(dataset = dataset_labels)) +
  scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
  scale_colour_gradient(low = "#2166ac", high = "#b2182b", na.value = "grey70",
                        breaks = c(0, 1), labels = c("Low", "High"), name = "Variable value") +
  labs(x = "SHAP value (impact on predicted yield, kg/ha)", y = NULL, title = "XGBoost: SHAP summary plot (test data)") +
  theme_pub
print(p_shap)
ggsave(file.path(pathOut, "xgb_shap_summary.png"), p_shap, width = 7, height = 5, dpi = 300)
