################################################################################
# STEP I - Harmonized TCD and TCA processing on 1km resolution
# A) Aggregating TCD (tree cover density) inside stable forest
# B) Calculating TCA (tree cover area)
#
################################################################################
#TCD change IN STABLE FOREST aggregated to 1000 m
tcd_change_stable_forest <- function(data_vector,
                                     ch_dir = "Z:/ch/",
                                     data_dir = paste0(ch_dir, "data/"),
                                     github_dir = "D:/github/clearing_house/")
{

  ##############################################################################
  # FOR PARALLEL PROCESSING:
  # LOAD PACKAGES
  req_pkg = c("dplyr", "terra", "parallel", "sf")
  lapply(req_pkg, require, character.only = TRUE)
  # LOAD FUNCTIONS
  source(paste0(github_dir, "processing/2_tcd_processing_Step_1_and_2_general_funs.R"))

  ##############################################################################
  # DATA INPUT:
  # COUNTING
  i <- data_vector$n
  # INPUT TCD DATA FOR 2012 AND 2018 (*1 to remove description)
  tcd12raw <- rast(data_vector$tcd12) * 1
  tcd18raw <- rast(data_vector$tcd18) * 1

  ##############################################################################
  # PROCESSING:
  # RETRIEVE FORESTED PIXELS FROM TCD 2012 AND 2018
  tcd12forest <- replaceNonForest(tcd12raw)
  tcd18forest <- replaceNonForest(tcd18raw)
  # aggregate tcd 2018 to same resolution as tcd 2012
  tcd18forest <- aggregate(tcd18forest, fact = 2, fun = mean, na.rm = TRUE)

  # MASK TCD 2018 WITH TCD 2012 AND VICE VERSA
  # (to prevent bias by regrowth or forest loss)
  tcd12forest <- mask(tcd12forest, tcd18forest)
  tcd18forest <- mask(tcd18forest, tcd12forest)

  # AGGREGATE TO 1KM
  tcd12agg <- aggregate(tcd12forest, fact = 50, fun = mean, na.rm = TRUE)
  tcd18agg <- aggregate(tcd18forest, fact = 50, fun = mean, na.rm = TRUE)

  # CALCULATE DELTA
  crs(tcd18agg) <- crs(tcd12agg)
  dTCD <- tcd18agg - tcd12agg

  ##############################################################################
  #OUTPUT
  raster_out(tcd12agg, paste0(data_dir, "TCD12_stable_forest/"), "tcd_2012_", i)
  raster_out(tcd18agg, paste0(data_dir, "TCD18_stable_forest/"), "tcd_2018_", i)
  raster_out(dTCD + 100, paste0(data_dir, "dTCD_stable_forest/"), "dTCD_", i)

  #CLEAN UP RAM
  clean_up(environment())
}


################################################################################
# TCA change aggregated to 1000 m
tca_change <- function(data_vector,
                       ch_dir = "Z:/ch/",
                       data_dir = paste0(ch_dir, "data/"),
                       github_dir = "D:/github/clearing_house/")
{

  ##############################################################################
  # FOR PARALLEL PROCESSING:
  # LOAD PACKAGES
  req_pkg = c("dplyr", "terra", "parallel", "sf")
  lapply(req_pkg, require, character.only = TRUE)
  # LOAD FUNCTIONS
  source(paste0(github_dir, "processing/2_tcd_processing_Step_1_and_2_general_funs.R"))

  ##############################################################################
  #PROCESSING
  i <- data_vector$n
  tca12 <- aggregate(rast(data_vector$tcd12) * 1,
                     fact = 50, fun = rel_rst) * 400
  na12 <- aggregate(rast(data_vector$tcd12) * 1,
                    fact = 50, fun = function(x) {length(x[x > 100])}) * 400
  threshold12 <- aggregate(rast(data_vector$tcd12) * 1,
                    fact = 50, fun = function(x) {length(x[x < 30])}) * 400
  #tca12[tca12 == 0] <- NA
  tca18 <- aggregate(rast(data_vector$tcd18) * 1,
                     fact = 100, fun = rel_rst) * 100
  na18 <- aggregate(rast(data_vector$tcd18) * 1,
                     fact = 100, fun = function(x) {length(x[x > 100])}) * 100
  threshold18 <- aggregate(rast(data_vector$tcd18) * 1,
                           fact = 100, fun = function(x) {length(x[x < 30])}) * 100
  #tca18[tca18 == 0] <- NA
  crs(tca18) <- crs(tca12)
  dTCA <- tca18 - tca12

  ##############################################################################
  #OUTPUT
  raster_out(tca12, paste0(data_dir, "TCA12/"), "tca_2012_", i) # values 0-100
  raster_out(na12, paste0(data_dir, "na12/"), "na_2012_", i) # values 0-100
  raster_out(threshold12, paste0(data_dir, "threshold12/"), "threshold_2012_", i) # values 0-100
  raster_out(tca18, paste0(data_dir, "TCA18/"), "tca_2018_", i) # values 0-100
  raster_out(na18, paste0(data_dir, "na18/"), "na_2018_", i) # values 0-100
  raster_out(threshold18, paste0(data_dir, "threshold18/"), "threshold_2018_", i) # values 0-100
  raster_out(dTCA + 100, paste0(data_dir, "dTCA/"), "dTCA_", i) # values 0-200
  #CLEAN UP
  clean_up(environment())
}


# generate sensitivity layers
sensitivity_lyrs <- function(data_vector,
                        ch_dir = "Z:/ch/",
                        data_dir = "D:/LandOeko/ch/data/",
                        github_dir = "D:/github/clearing_house/")
{

  ##############################################################################
  # FOR PARALLEL PROCESSING:
  # LOAD PACKAGES
  req_pkg = c("dplyr", "terra", "parallel", "sf")
  lapply(req_pkg, require, character.only = TRUE)
  # LOAD FUNCTIONS
  source(paste0(github_dir, "processing/2_tcd_processing_Step_1_and_2_general_funs.R"))

  ##############################################################################
  #PROCESSING
  i <- data_vector$n

  na12 <- aggregate(rast(data_vector$tcd12) * 1,
                    fact = 50, fun = function(x) {length(x[x > 100])}) * 400
  threshold12 <- aggregate(rast(data_vector$tcd12) * 1,
                           fact = 50, fun = function(x) {length(x[x < 30])}) * 400
  #tca12[tca12 == 0] <- NA

  na18 <- aggregate(rast(data_vector$tcd18) * 1,
                    fact = 100, fun = function(x) {length(x[x > 100])}) * 100
  threshold18 <- aggregate(rast(data_vector$tcd18) * 1,
                           fact = 100, fun = function(x) {length(x[x < 30])}) * 100
  #tca18[tca18 == 0] <- NA


  ##############################################################################
  #OUTPUT
  raster_out(na12, paste0(data_dir, "na12/"), "na_2012_", i, dtype = "INT2U") # values 0-100
  raster_out(threshold12, paste0(data_dir, "threshold12/"), "threshold_2012_",
             i, dtype = "INT2U") # values 0-100
  raster_out(na18, paste0(data_dir, "na18/"), "na_2018_", i, dtype = "INT2U") # values 0-100
  raster_out(threshold18, paste0(data_dir, "threshold18/"), "threshold_2018_",
             i, dtype = "INT2U") # values 0-100
  #CLEAN UP
  clean_up(environment())
}

################################################################################
# END
################################################################################

