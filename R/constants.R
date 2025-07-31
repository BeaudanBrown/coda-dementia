m <- 1 # Number of imputed datasets
maxit <- 10 # Number of MICE iterations
n_boots <- 500

hrs_in_day <- 24
mins_in_day <- 1440

short_sleep_hours <- 6
long_sleep_min_hours <- 8
long_sleep_hours <- 9

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

cohorts <- list(
  filter_fn = rlang::syms(c(
    "short_sleeper_filter_fn",
    "avg_sleeper_filter_fn",
    "long_sleeper_filter_fn"
  )),
  cohort = c("short_sleeper", "avg_sleeper", "long_sleeper"),
  colour = c("#ff747b", "#6ed853", "#708ff9")
)
