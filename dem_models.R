# load packages
library(mice)
library(tidyverse)
library(survival)
library(rms)
library(compositions)

# Load environment variables from the .env file
dotenv::load_dot_env()
data_dir <- Sys.getenv("DATA_DIR")

## Define SBP 

sbp <- matrix(
  c(
    1, 1, -1, -1,
    1, -1, 0, 0,
    0, 0, 1, -1
  ),
  ncol = 4, byrow = TRUE
)

v <- gsi.buildilrBase(t(sbp))

## Load data
df <- read_rds(file.path(data_dir, "bootstrap_data.rds"))

## Fit imputation model using MICE (single imputation)

# set date variables to strings to avoid errors

df$date_accel <- as.character(df$date_accel)
df$date_acdem2 <- as.character(df$date_acdem2)
df$date_of_death <- as.character(df$date_of_death)

# Matrix of variables to include in imputation model

predmat <- quickpred(df, 
                     mincor = 0, 
                     exclude = c("date_acdem2", "date_accel", "date_of_death", 
                                "avg_sleep", "avg_inactivity", "avg_light", 
                                "avg_mvpa"))

# exclude dates from being imputed
predmat["date_acdem2",] <- 0
predmat["date_of_death",] <- 0

#imp <- mice(df, predictorMatrix = predmat, m = 1)

#write_rds(imp, file.path(data_dir, "24hr_behaviours/imp.rds"))

read_rds(file.path(data_dir, "24hr_behaviours/imp.rds"))

# extract imputed dataset

imp <- complete(imp)

## Create person-period dataset with age as timescale

imp$id <- 1:nrow(imp)

imp_long <- survSplit(Surv(time = age_accel, event = dem, time2 = age_dem) ~ .,
  data = imp,
  cut = seq(
    from = min(imp$age_dem),
    to = max(imp$age_dem),
    length.out = 76
  ),
  episode = "timegroup", end = "age_end", event = "dem",
  start = "age_start"
)

## Fit model ##

m1 <- 
  glm(dem ~ rcs(age_start,5) + pol(R1,2) + pol(R2,2) + pol(R3,2) + 
        rcs(bp_syst_avg,3) + sex + retired + shift + apoe_e4 + highest_qual +
        rcs(townsend_deprivation_index,3) + antidepressant_med +
        antipsychotic_med + insomnia_med + ethnicity + avg_total_household_income + 
        smok_status, family = binomial, data = imp_long)

imp_stacked <- do.call("rbind", replicate(76, imp, simplify = FALSE)) 

imp_stacked$age_start <- as.integer(
  rep(seq(min(imp_long$age_dem),
          max(imp_long$age_dem),
          length.out = 76),
      nrow(imp)))

# reference composition 

ref_comp <-
  acomp(data.frame(
    mean(imp$avg_sleep), mean(imp$avg_inactivity),
    mean(imp$avg_light), mean(imp$avg_mvpa)
  ))

ref_ilr <-
  ilr(ref_comp, V = v) |>
  setNames(c("R1", "R2", "R3"))

ref_data <- 
  mutate(imp_stacked, 
         R1 = ref_ilr[1], R2 = ref_ilr[2], R3 = ref_ilr[3])

ref_data$haz <- 
  predict(m1, newdata = ref_data, type = "response")

ref_data <- ref_data |> 
  group_by(id) |> 
  arrange(age_start) |> 
  mutate(risk = 1 - cumprod(1-haz)) |> 
  ungroup()

risk_ref <- ref_data |> 
  filter(age_start == 75) |> 
  summarise(mean = mean(risk))

risk_ref

# substitution 

sub_comp <-
  acomp(data.frame(
    mean(imp$avg_sleep)+60, mean(imp$avg_inactivity),
    mean(imp$avg_light), mean(imp$avg_mvpa)-60))

sub_ilr <-
  ilr(sub_comp, V = v) |>
  setNames(c("R1", "R2", "R3"))

sub_data <- 
  mutate(imp_stacked, 
         R1 = sub_ilr[1], R2 = sub_ilr[2], R3 = sub_ilr[3])

sub_data$haz <- 
  predict(m1, newdata = sub_data, type = "response")

sub_data <- sub_data |> 
  group_by(id) |> 
  arrange(age_start) |> 
  mutate(risk = 1 - cumprod(1-haz)) |> 
  ungroup()

risk_sub <- sub_data |> 
  filter(age_start == 75) |> 
  summarise(mean = mean(risk))

risk_sub / risk_ref