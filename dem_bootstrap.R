source("dem_models.R")

## Load data
boot_data <- read_rds(file.path(data_dir, "bootstrap_data.rds"))

# set date variables to strings to avoid errors
boot_data$date_accel <- as.character(boot_data$date_accel)
boot_data$date_acdem2 <- as.character(boot_data$date_acdem2)
boot_data$date_of_death <- as.character(boot_data$date_of_death)

## UNCOMMENT FOR TEST RUNS
# boot_data <- boot_data[1:1000, ]

# Run bootstrap for primary model
run_primary_bootstrap <- function(boot_data) {
  run_bootstrap(
    boot_data = boot_data,
    timegroup = 55,
    create_formula_fn = get_primary_formula,
    output_name = "boot_primary.rds"
  )
}

# Run bootstrap for sensitivity analysis model 1
run_s1_bootstrap <- function(boot_data) {
  run_bootstrap(
    boot_data = boot_data,
    timegroup = 55,
    create_formula_fn = get_s1_formula,
    output_name = "boot_s1.rds"
  )
}

# Run bootstrap for sensitivity analysis model 2
run_s2_bootstrap <- function(boot_data) {
  run_bootstrap(
    boot_data = boot_data,
    timegroup = 55,
    create_formula_fn = get_s2_formula,
    output_name = "boot_s2.rds"
  )
}

## Run bootstrap for sensitivity analysis model 3
# excluding first four years of follow-up

run_s3_bootstrap <- function(boot_data) {
  subset_data <- boot_data[boot_data$dem == 0 | boot_data$time_to_dem > (4 * 365), ]

  run_bootstrap(
    boot_data = subset_data,
    timegroup = 55,
    create_formula_fn = get_primary_formula,
    output_name = "boot_s3_excludefirst4.rds"
  )
}

# Run bootstrap for a particular model formula and timegroup target, outputting to a file
run_bootstrap <- function(boot_data, timegroup, create_formula_fn, output_name) {
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

  ## reference compositions
  all_comp <- acomp(boot_data[, c("avg_sleep", "avg_inactivity", "avg_light", "avg_mvpa")])

  short_sleep_comp <-
    all_comp[all_comp$avg_sleep < short_sleep_hours / hrs_in_day, ]
  short_sleep_geo_mean <-
    acomp(apply(short_sleep_comp, 2, function(x) exp(mean(log(x)))))

  avg_sleep_comp <-
    all_comp[all_comp$avg_sleep >= short_sleep_hours / hrs_in_day, ]
  avg_sleep_geo_mean <-
    acomp(apply(avg_sleep_comp, 2, function(x) exp(mean(log(x)))))

  result <- boot(
    data = boot_data,
    statistic = bootstrap_substitutions_fn,
    create_formula_fn = create_formula_fn,
    timegroup = timegroup,
    predmat = predmat,
    imp_methods = imp_methods,
    short_sleep_geo_mean = short_sleep_geo_mean,
    avg_sleep_geo_mean = avg_sleep_geo_mean,
    R = bootstrap_iterations,
    parallel = "multicore",
    ncpus = ncpus
  )

  # Prepend timestamp to avoid accidental data loss
  timestamp <- format(Sys.time(), "%Y-%m-%d_%H:%M")
  output_name_with_timestamp <- paste0(output_name, "_", timestamp, ".rds")
  saveRDS(result, file.path(output_dir, output_name_with_timestamp))
}

# Ran each bootstrap iteration
# Imputes the data, fits the model for the given formula and predicts the risks for a
# variety of substitutions and cohorts
bootstrap_substitutions_fn <- function(
  data,
  indices,
  create_formula_fn,
  timegroup,
  predmat,
  imp_methods,
  short_sleep_geo_mean,
  avg_sleep_geo_mean
) {
  this_sample <- data[indices, ]

  print("Imputing")
  print(format(Sys.time(), "%H:%M:%S"))
  imp <- mice(this_sample, m = 1, predictorMatrix = predmat, methods = imp_methods)
  imp <- complete(imp)
  setDT(imp)
  imp[, id := .I]
  imp_len <- nrow(imp)
  # fit model
  print("Fitting model")
  print(format(Sys.time(), "%H:%M:%S"))
  print(gc())
  model <- fit_model(imp, create_formula_fn)

  # data for g-computation/standardisation
  imp <- imp[rep(seq_len(imp_len), each = timegroup)]
  imp[, timegroup := rep(1:timegroup, imp_len)]

  print("Generating predictions")
  print(format(Sys.time(), "%H:%M:%S"))
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

  full_df <- full_join(short_sleep_inactive, short_sleep_light, by = "offset") |>
    full_join(short_sleep_mvpa, by = "offset") |>
    full_join(avg_sleep_inactive, by = "offset") |>
    full_join(avg_sleep_light, by = "offset") |>
    full_join(avg_sleep_mvpa, by = "offset")

  print("Returning results")
  print(format(Sys.time(), "%H:%M:%S"))
  return(as.matrix(full_df))
}

process_boot_output <- function(rds_path) {
  data <- readRDS(file.path(data_dir, rds_path))
  num_subs <- sub_steps * 2 + 1
  sub_col_names <- data$t[1, 1:num_subs]

  plot_data <- pivot_longer(as.data.frame(data$t0), -offset, values_to = "risk", names_to = "Substitution") |>
    group_by(Substitution) |>
    mutate(ref_risk = ifelse(offset == 0, risk, NA_real_)) |>
    fill(ref_risk, .direction = "downup") |>
    mutate(risk_dif = risk - ref_risk,
           risk_ratio = risk / ref_risk) |>
    mutate(Reference = ifelse(str_detect(Substitution, "short_sleep"),
                              "Short sleepers", "Normal sleepers")) |>
    mutate(Substitution = str_remove(Substitution, "_short_sleep_geo_mean")) |>
    mutate(Substitution = str_remove(Substitution, "_avg_sleep_geo_mean")) |>
    mutate(Substitution = ifelse(Substitution == "avg_inactivity", "Inactivity",
                                 ifelse(Substitution == "avg_light", "Light activity", "MVPA")))

  get_quantiles <- function(start, substitution, reference) {
    slice <- data$t[, start:(start + num_subs - 1)]
    middle_col <- ceiling(ncol(slice) / 2)

    slice <- t(apply(slice, 1, function(row, idx) {
      zero_offset <- row[middle_col]
      return(row / zero_offset)
    }, idx = middle_col))

    quantiles <-
      as.data.frame(t(apply(slice, 2, function(column) quantile(column, probs = c(0.025, 0.975)))))
    colnames(quantiles) <- c("lower", "upper")
    quantiles$Substitution <- substitution
    quantiles$Reference <- reference
    quantiles$offset <- sub_col_names
    return(quantiles)
  }

  sub_start <- num_subs + 1

  inactivity_short_sleep <- get_quantiles(sub_start, "Inactivity", "Short sleepers")
  sub_start <- sub_start + num_subs

  light_short_sleep <- get_quantiles(sub_start, "Light activity", "Short sleepers")
  sub_start <- sub_start + num_subs

  mvpa_short_sleep <- get_quantiles(sub_start, "MVPA", "Short sleepers")
  sub_start <- sub_start + num_subs

  inactivity_avg_sleep <- get_quantiles(sub_start, "Inactivity", "Normal sleepers")
  sub_start <- sub_start + num_subs

  light_avg_sleep <- get_quantiles(sub_start, "Light activity", "Normal sleepers")
  sub_start <- sub_start + num_subs

  mvpa_avg_sleep <- get_quantiles(sub_start, "MVPA", "Normal sleepers")

  all_quantiles <- rbind(inactivity_short_sleep,
                         light_short_sleep,
                         mvpa_short_sleep,
                         inactivity_avg_sleep,
                         light_avg_sleep,
                         mvpa_avg_sleep)


  plot_data <- full_join(plot_data, all_quantiles, by = c("offset", "Substitution", "Reference"))

  rr_plot <- function(sub, refcomp, colour) {
    
    plot_data$Substitution <- 
      ifelse(plot_data$Substitution == "Inactivity","inactivity", 
             ifelse(plot_data$Substitution == "Light activity",
                               "light activity", "MVPA"))

    plot_data2 <- plot_data[plot_data$Substitution == sub &
                plot_data$Reference == refcomp,]
    

    plot_data2 |>
      ggplot(aes(x = offset, y = risk_ratio)) +
      geom_line(colour = colour) +
      geom_ribbon(aes(ymin = lower, ymax = upper),
                  alpha = 0.25, fill = colour) +
      facet_wrap(~ Substitution, nrow = 2) +
      geom_hline(yintercept = 1, linetype = "dotted") +
      xlab("") +
      ylab("Risk ratio") +
      annotate(geom = "text", x=0, y = -0.15,
               hjust = 0.5,fontface = 1,size = 14/.pt,
               label = "Minutes", family = "serif") +
      annotate(geom = "text", x=-20, y = -0.5,
               hjust = 1,fontface = 1,size = 12/.pt,
               label = "Less sleep", family = "serif") +
      annotate(geom = "text", x=20, y = -0.5,
               hjust = 0,fontface = 1, size = 12/.pt,
               label = "More sleep", family = "serif") +
      geom_segment(aes(x = 1, y = -0.625, 
                       xend = 15, yend = -0.625),
                   arrow = arrow(length = unit(0.15, "cm"))) +
      geom_segment(aes(x = -1, y = -0.625, 
                       xend = -15, yend = -0.625),
                   arrow = arrow(length = unit(0.15, "cm"))) +
      annotate(geom = "text", x=-20, y = -0.75,
               hjust = 1, size = 12/.pt,
               label = paste("More", plot_data2$Substitution),
               family = "serif", fontface = 1, size = 12/.pt) +
      annotate(geom = "text", x=20, y = -0.75,
               hjust = 0,
               label = paste("Less", plot_data2$Substitution),
               family = "serif", fontface = 1, size = 12/.pt) +
      coord_cartesian(ylim = c(0.33, 3), expand = FALSE, clip = "off") +
      cowplot::theme_cowplot() +
      theme(text=element_text(size = 12, family="serif"),
            plot.margin = unit(c(1, 1, 4, 1), "lines"),
            strip.background = element_blank(),
            strip.text.x = element_blank())
  }

  # normal sleepers
  p1 <- rr_plot("inactivity", "Normal sleepers", "#fc020f")
  p2 <- rr_plot("light activity", "Normal sleepers", "#145e01")
  p3 <- rr_plot("MVPA", "Normal sleepers", "#011869")

  pnorm <-
    plot_grid(p1 + labs(x = "",title = "Normal sleepers") + 
                theme(legend.position = "none", plot.title = element_text(hjust = 0.5)),
              p2 + labs(x = "", ) + theme(legend.position = "none"),
              p3 + labs(x = "") + theme(legend.position = "none"),
              align = "vh",
              nrow = 3)

  # short sleepers
  p4 <- rr_plot("inactivity", "Short sleepers", "#ff747b")
  p5 <- rr_plot("light activity", "Short sleepers", "#6ed853")
  p6 <- rr_plot("MVPA", "Short sleepers", "#708ff9")

  pshort <-
    plot_grid(p4 + labs(title = "Short sleepers", y = "") + 
                theme(legend.position="none", plot.title = element_text(hjust = 0.5)),
              p5 + labs(y = "") +
                theme(legend.position="none"),
              p6 + labs(y = "") + theme(legend.position = "none"),
              align = "vh",
              nrow = 3)

  plot <- plot_grid(pnorm,
            pshort,
            nrow = 1)
  
  return(list(plot, plot_data))
}

#### Results 

# ## Primary model
# 
process_boot_output("boot_primary.rds")[[1]]

# save
ggsave(
  file.path(
    data_dir,
    "../../Papers/Substitution Analysis/Main_figures/Substitutions.png"
  ),
  device = "png",
  bg = "white",
  width = 10,
  height = 12,
  dpi = 500
)
# 
# # risk ratios 
# 
# process_boot_output("boot_primary.rds")[[2]] |> 
#   filter(abs(offset) == 60) |> 
#   filter(Reference == "Short sleepers")
# 
# 
# ## Sensitivity 1
# 
# process_boot_output("boot_s1.rds")[[1]]
# 
# ggsave(
#   file.path(
#     data_dir,
#     "../../Papers/Substitution Analysis/Appendix_figures/Sensitivity_1.png"
#   ),
#   device = "png",
#   bg = "white",
#   width = 10,
#   height = 12
# )

# 
# ## Sensitivity 2
# 
# process_boot_output("boot_s2.rds")[[1]]
# 
# ggsave(
#   file.path(
#     data_dir,
#     "../../Papers/Substitution Analysis/Appendix_figures/Sensitivity_2.png"
#   ),
#   device = "png",
#   bg = "white",
#   width = 10,
#   height = 12
# )
