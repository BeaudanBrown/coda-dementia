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
predmat <- quickpred(boot_data,
                     mincor = 0,
                     exclude = c(
                       "date_acdem2", "date_accel", "date_of_death",
                       "avg_sleep", "avg_inactivity", "avg_light",
                       "avg_mvpa", "eid"
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

### Fit outcome model 
model <- fit_model_timescale2(imp, get_primary_formula_timescale2)

### Estimate substitution effects

# data for g-computation/standardisation

pred_data <- imp[rep(1,34),]

pred_data <- pred_data |> 
  select(sex,
         R1, R2, R3,
           retired,
           shift,
           apoe_e4,
           highest_qual, townsend_deprivation_index,
           antidepressant_med,
           antipsychotic_med,
           insomnia_med,
           ethnicity,
           avg_total_household_income,
           smok_status)

pred_data$timegroup <- 1:34

setDT(pred_data)

all_comp <- acomp(boot_data[, c("avg_sleep", "avg_inactivity", "avg_light", "avg_mvpa")])


# test

test <- predict_composition_risk(all_comp[1,], pred_data, model)

predict(model, newdata = imp[1,])



# estimate substitutions 

short_sleep_inactive <-
  calc_substitution(short_sleep_geo_mean,
                    imp,
                    model,
                    c("avg_sleep", "avg_inactivity"),
                    timegroup = timegroup)

short_sleep_light <-
  calc_substitution(short_sleep_geo_mean,
                    imp,
                    model,
                    c("avg_sleep", "avg_light"),
                    timegroup = timegroup)

short_sleep_mvpa <-
  calc_substitution(short_sleep_geo_mean,
                    imp,
                    model,
                    c("avg_sleep", "avg_mvpa"),
                    timegroup = timegroup)

avg_sleep_inactive <-
  calc_substitution(avg_sleep_geo_mean,
                    imp,
                    model,
                    c("avg_sleep", "avg_inactivity"),
                    timegroup = timegroup)

avg_sleep_light <-
  calc_substitution(avg_sleep_geo_mean,
                    imp,
                    model,
                    c("avg_sleep", "avg_light"),
                    timegroup = timegroup)

avg_sleep_mvpa <-
  calc_substitution(avg_sleep_geo_mean,
                    imp,
                    model,
                    c("avg_sleep", "avg_mvpa"),
                    timegroup = timegroup)

# pool results 

full_df <- full_join(short_sleep_inactive, short_sleep_light, by = "offset") |>
  full_join(short_sleep_mvpa, by = "offset") |>
  full_join(avg_sleep_inactive, by = "offset") |>
  full_join(avg_sleep_light, by = "offset") |>
  full_join(avg_sleep_mvpa, by = "offset")


