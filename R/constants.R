m <- 1 # Number of imputed datasets
maxit <- 10 # Number of MICE iterations
n_boots <- 500

hrs_in_day <- 24
mins_in_day <- 1440

short_sleep_hours <- 6
long_sleep_min_hours <- 8
long_sleep_hours <- 9
intervention_threshold <- 0.75

synth_threshold <- 0.05
uni_threshold <- 0.025

no_filter_fn <- function(df) {
  df
}

not_retired_filter_fn <- function(df) {
  df |> filter(retired == 0)
}

retired_filter_fn <- function(df) {
  df |> filter(retired == 1)
}


short_sleeper_filter_fn <- function(df) {
  df |> filter(avg_sleep < short_sleep_hours * 60)
}

avg_sleeper_filter_fn <- function(df) {
  df |>
    filter(
      avg_sleep >= short_sleep_hours * 60 & avg_sleep <= long_sleep_hours * 60
    )
}

long_sleeper_filter_fn <- function(df) {
  df |> filter(avg_sleep > long_sleep_min_hours * 60)
}

substitutions <- expand.grid(
  from_var = c("avg_mvpa", "avg_light", "avg_inactivity"),
  to_var = "avg_sleep",
  stringsAsFactors = FALSE
)

durations <- data.frame(
  duration = seq(from = -60, to = 60, by = 15)[
    seq(from = -60, to = 60, by = 15) != 0
  ]
)
all_subs <- merge(substitutions, durations)


retired_cohorts <- list(
  filter_fn = rlang::syms(c(
    "not_retired_filter_fn",
    "retired_filter_fn"
  )),
  cohort = c("not_retired", "retired"),
  colour = c("#ff747b", "#6ed853")
)

cohorts <- list(
  filter_fn = rlang::syms(c(
    "short_sleeper_filter_fn",
    "avg_sleeper_filter_fn",
    "long_sleeper_filter_fn",
    "no_filter_fn"
  )),
  cohort = c("short_sleeper", "avg_sleeper", "long_sleeper", "full_cohort"),
  colour = c("#ff747b", "#6ed853", "#708ff9", "#f39c12")
)
