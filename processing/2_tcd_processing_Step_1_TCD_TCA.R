################################################################################
# STEP I - Harmonized TCD and TCA processing on 1km resolution
# A) Aggregating TCD (tree cover density) inside stable forest
# B) Calculating TCA (tree cover area)
#
################################################################################
# Functions and directories
github_dir <- "D:/github/clearing_house/"
source(paste0(github_dir, "processing/2_tcd_processing_Step_1_and_2_general_funs.R"))
source(paste0(github_dir, "processing/2_tcd_processing_Step_1_funs.R"))

ch_dir <- "Z:/ch/"
data_dir <- paste0(ch_dir, "data/")


################################################################################
# List input files for 2012 and 2018 and gather them in a list
dat_df <- mkDatDf(data_dir)


################################################################################
# DEBUGGING INPUT
#data_list <- split(dat_df, dat_df$n)
#data_vector <- data_list[[221]]


################################################################################
#EXECUTE TCD CHANGE CALCULATION
to_do <- find_missing(paste0(data_dir, "TCD12_stable_forest/"), nrow(dat_df))
if (!is.null(to_do)) data_list <- dat_df[to_do,] %>% split(., .$n)
par_exec(data_list, tcd_change_stable_forest, perc_cores = .3)


################################################################################
#EXECUTE TCA CHANGE CALCULATION
to_do <- find_missing(paste0(data_dir, "TCA12/"), nrow(dat_df))
if (!is.null(to_do)) data_list <- dat_df[to_do,] %>% split(., .$n)
par_exec(data_list, fun = tca_change, perc_cores = .6)


###############################################################################
#EXECUTE SENSITIVITY LaYeR CALCULATION
to_do <- find_missing(paste0(data_dir, "na12/"), nrow(dat_df))
if (!is.null(to_do)) data_list <- dat_df[to_do,] %>% split(., .$n)
par_exec(data_list, fun = sensitivity_lyrs, perc_cores = .6)
################################################################################
# END
################################################################################

