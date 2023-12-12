source("dem_models_timescale2.R")

# Constants
short_sleep_hours <- 6
hrs_in_day <- 24

# Load environment variables from the .env file
dotenv::load_dot_env()
data_dir <- Sys.getenv("DATA_DIR")
output_dir <- Sys.getenv("OUTPUT_DIR")

###  Load data
boot_data <- read_rds(file.path(data_dir, "bootstrap_data.rds"))
# set date variables to strings to avoid errors
boot_data$date_accel <- as.character(boot_data$date_accel)
boot_data$date_acdem2 <- as.character(boot_data$date_acdem2)
boot_data$date_of_death <- as.character(boot_data$date_of_death)


### fit imputation model 

# Matrix of variables to include in imputation model
predmat <- quickpred(
  boot_data,
  mincor = 0,
  exclude = c(
    "date_acdem2",
    "date_accel",
    "date_of_death",
    "avg_sleep",
    "avg_inactivity",
    "avg_light",
    "avg_mvpa",
    "eid"
  )
)

predmat["date_acdem2", ] <- 0
predmat["date_of_death", ] <- 0

# method for each imputed variable
imp_methods <- make.method(boot_data)
# exclude dates from being imputed
imp_methods["date_acdem2"] <- ""
imp_methods["date_of_death"] <- ""

# imputation 
#imp <- mice(boot_data, m = 1, predictorMatrix = predmat, methods = imp_methods)
#imp <- complete(imp)

#write_rds(imp, file.path(data_dir, "imp_timescale.rds"))

imp <- read_rds(file.path(data_dir, "imp_timescale.rds"))


### Create long (person-period) dataset

imp_long <- survSplit(
  Surv(time = time_to_dem, event = dem) ~ .,
  data = imp,
  cut = seq(
    from = min(imp$time_to_dem),
    to = max(imp$time_to_dem),
    length.out = 34
  ),
  episode = "timegroup",
  end = "time_end",
  event = "dem",
  start = "time_start"
)

# datadist

dd <- datadist(imp_long)
options(datadist = "dd")

  
### Fit model 

knots_timegroup <-
  quantile(imp_long[["timegroup"]], c(0.05, 0.35, 0.65, 0.95))
knots_deprivation <-
  quantile(imp_long[["townsend_deprivation_index"]], c(0.1, 0.5, 0.9))

# %ia% excludes higher order product terms 

model <- lrm(
  dem ~
    rcs(timegroup, knots_timegroup) +
    pol(R1) + 
    pol(R2) + 
    pol(R3) +
    pol(R1) %ia% rcs(timegroup, knots_timegroup) +
    pol(R2) %ia% rcs(timegroup, knots_timegroup) +
    pol(R3) %ia% rcs(timegroup, knots_timegroup) +
    sex +
    retired +
    shift +
    apoe_e4 +
    highest_qual +
    rcs(townsend_deprivation_index, knots_deprivation) +
    antidepressant_med +
    antipsychotic_med +
    insomnia_med +
    ethnicity +
    avg_total_household_income +
    smok_status,
  data = imp_long
)


### Estimate substitution effects

get_hr <- function(model, ilr_sub, ilr_ref) {
  out <- contrast(
    model,
    list(
      R1 = ilr_sub[1],
      R2 = ilr_sub[2],
      R3 = ilr_sub[3],
      timegroup = 1:34
    ),
    list(
      R1 = ilr_ref[1],
      R2 = ilr_ref[2],
      R3 = ilr_ref[3],
      timegroup = 1:34
    )
  )
  
  return(tibble(
    timegroup = out$timegroup,
    Contrast = out$Contrast,
    Lower = out$Lower,
    Upper = out$Upper
  ))
}

# test 

get_hr(model, all_comp[1,], all_comp[2,])


## can either wrangle the above function so that the right ILRs are passed in, or ditch it entirely and just adapt calc_substitution 
