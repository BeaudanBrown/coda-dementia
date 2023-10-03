# load packages 
library(tidyverse)
library(dotenv)
library(data.table)
library(compositions)
library(mice)

# Load environment variables from the .env file
dotenv::load_dot_env()
data_dir <- Sys.getenv("DATA_DIR")

## Read RDS 



## Fit imputation model using MICE (single imputation)

#imp <- mice(dem_model_data, m = 1)

#write_rds(imp, file.path(data_dir, "24hr_behaviours/imp.rds"))

read_rds(file.path(data_dir, "24hr_behaviours/imp.rds")))

# extract imputed dataset

imp <- complete(imp)

## Constructing risk set ##

# add in dates

imp$date_of_death <- dem_df$date_of_death
imp$date_acdem2 <- dem_df$date_acdem2
imp$date_accel <- dem_df$calendar_date

# Death from other causes remain in risk set until end of FU 
# See Young et al (2019) 

imp <- imp |> 
  mutate(competing = ifelse(!is.na(date_of_death) & is.na(date_acdem2), 1, 0)) |> 
  mutate(time_to_dem = case_when(
    dem == 1 ~ difftime(date_acdem2, date_accel),
    competing == 1 ~ difftime("2022-01-01", date_accel),
    TRUE ~ difftime("2022-01-01", date_accel))) |> 
  mutate(time_to_dem = as.integer(time_to_dem))

# Create age at dementia or censoring/end of follow-up variable 

imp$age_dem <- imp$age_accel + (imp$time_to_dem / 365)

## Create person-period dataset with age as timescale 

imp_long <- survSplit(Surv(time=age_accel, event=dem, time2=age_dem)~., 
                      data = imp, 
                      cut=seq(from=min(imp$age_dem), 
                              to=max(imp$age_dem), 
                              length.out = 76),
                      episode = "timegroup", end = "age_dem", event = "dem",
                      start = "age_accel")

write_rds(imp_long)
