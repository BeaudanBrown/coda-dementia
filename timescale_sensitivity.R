source("dem_models.R")
source("utils.R")

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
# imp <- mice(boot_data, m = 1, predictorMatrix = predmat, methods = imp_methods)
# imp <- complete(imp)

# write_rds(imp, file.path(data_dir, "imp_timescale.rds"))

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

### Function to fit models

fit_model <- function(low, high) {
  df <- imp_long |> filter(timegroup >= low & timegroup < high)

  # datadist
  dd <- datadist(df)
  options(datadist = dd)

  # model
  model <- lrm(
    dem ~
      rcs(timegroup, 4) +
      pol(R1) +
      pol(R2) +
      pol(R3) +
      sex +
      retired +
      shift +
      apoe_e4 +
      highest_qual +
      rcs(townsend_deprivation_index, 3) +
      psych_meds +
      ethnicity +
      avg_total_household_income +
      smok_status,
    data = df
  )

  return(model)
}


### Estimate substitution effects

# reference compositions

all_comp <- acomp(boot_data[, c("avg_sleep", "avg_inactivity", "avg_light", "avg_mvpa")])

short_sleep_comp <-
  all_comp[all_comp$avg_sleep < short_sleep_hours / hrs_in_day, ]

short_sleep_geo_mean <-
  acomp(apply(short_sleep_comp, 2, function(x) exp(mean(log(x)))))

avg_sleep_comp <-
  all_comp[all_comp$avg_sleep >= short_sleep_hours / hrs_in_day, ]

avg_sleep_geo_mean <-
  acomp(apply(avg_sleep_comp, 2, function(x) exp(mean(log(x)))))


## Function for hazard ratios

get_hr <- function(model, ilr_sub, ilr_ref) {
  out <- contrast(
    model,
    list(
      R1 = ilr_sub[1],
      R2 = ilr_sub[2],
      R3 = ilr_sub[3]
    ),
    list(
      R1 = ilr_ref[1],
      R2 = ilr_ref[2],
      R3 = ilr_ref[3]
    )
  )

  return(tibble(
    Contrast = out$Contrast,
    Lower = out$Lower,
    Upper = out$Upper
  ))
}


## function to pass in substitutions to above function

get_sub <- function(model, base_comp, substitution) {
  # increment for substitution
  inc <- c(-1 / 24, 1 / 24)

  # initialise list for sub vectors
  sub_comps_list <- vector("list", length(inc))

  # Loop over inc and create a sub_comps data table for each element
  for (i in seq_along(inc)) {
    # The list of compositions to be fed into the model after applying the substitutions
    sub_comps <- as.data.table(t(base_comp))
    setnames(sub_comps, c("avg_sleep", "avg_inactivity", "avg_light", "avg_mvpa"))

    sub_comps[, (substitution[1]) := .SD[[substitution[1]]] + inc[i]]
    sub_comps[, (substitution[2]) := .SD[[substitution[2]]] - inc[i]]

    # Store the data.table into list
    sub_comps_list[[i]] <- sub_comps
  }

  # Combine all data.tables in the list
  sub_comps <- rbindlist(sub_comps_list)

  # get vector of HRs for each substitution

  sub_hrs <-
    rbindlist(lapply(
      seq_len(nrow(sub_comps)),
      function(i) get_hr(model, ilr(acomp(sub_comps[i]), V = v), base_comp)
    ))

  sub_hrs$Substitution <- substitution[2]
  sub_hrs$Shift <- c(inc[1], inc[2]) * 24

  return(sub_hrs)
}


### Split time group up into 3 chunks and fit modelin each

timegroup_chunks <-
  quantile(imp_long[imp_long$dem == 1, ]$timegroup, c(0.5))

low <- c(0, timegroup_chunks[1])
high <- c(timegroup_chunks[1], 50)

models <- map2(low, high, fit_model)


### estimate hazard ratios in each chunk for each substitution

get_stratified_hr <- function(substitution) {
  sub_res <-
    rbindlist(
      lapply(
        models,
        get_sub,
        base_comp = avg_sleep_geo_mean,
        substitution = substitution
      )
    )

  sub_res$time_chunk <- rep(1:2, each = 2)

  return(sub_res)
}


# sleep modvig

all_subs <- list(
  c("avg_sleep", "avg_mvpa"),
  c("avg_sleep", "avg_light"),
  c("avg_sleep", "avg_inactivity")
)

out <- rbindlist(
  lapply(
    all_subs, get_stratified_hr
  )
)


### Plot

out |>
  mutate(across(c(Contrast, Lower, Upper), exp)) |>
  mutate(Shift = ifelse(Shift == -1, "Add 1 hr sleep", "Remove 1 hr sleep")) |>
  ggplot(aes(x = time_chunk, y = Contrast, colour = Shift)) +
  geom_pointrange(aes(ymin = Lower, ymax = Upper),
    position = position_dodge(0.15)
  ) +
  geom_hline(aes(yintercept = 1), linetype = "dashed") +
  facet_wrap(~Substitution, ncol = 1) +
  labs(x = "Time since accelerometry", y = "HR") +
  scale_x_continuous(
    breaks = c(1, 2),
    labels = c("<5 years", "\u2265 5 years")
  ) +
  theme_cowplot()

ggsave(file.path(data_dir, "stratified_hrs.png"),
  device = "png", bg = "white", width = 6,
  height = 8
)


# #### Bootstrap
#
# ### Create long (person-period) dataset
#
# boot_sub <- function(data,indices){
#
#   df_sample <- imp[indices,]
#
#   imp_long <- survSplit(
#     Surv(time = time_to_dem, event = dem) ~ .,
#     data = df_sample,
#     cut = seq(
#       from = min(df_sample$time_to_dem),
#       to = max(df_sample$time_to_dem),
#       length.out = 34
#     ),
#     episode = "timegroup",
#     end = "time_end",
#     event = "dem",
#     start = "time_start"
#   )
#
#   # fit stratified models
#
#   models <- map2(low, high, fit_model)
#
#   out <- rbindlist(
#     lapply(
#       all_subs, get_stratified_hr
#     )
#   )
#
#   return(out$Contrast)
# }
#
#
# ## run bootstrap locally
#
# boot_out <- boot(imp, boot_sub, R = 50)
