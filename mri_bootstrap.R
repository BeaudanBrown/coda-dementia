source("mri_models.R")
library(tidyverse)

## Read MRI data

mri_df <-
  fread(file.path(data_dir, "mri_data.csv"),
    stringsAsFactors = TRUE
  ) |>
  as_tibble()

mri_df$assessment_centre_mri1 <- as.factor(mri_df$assessment_centre_mri1)

## Read main data
boot_df <- read_rds(file.path(data_dir, "bootstrap_data.rds"))

## Merge

mri_df <- left_join(mri_df, boot_df, by = "eid")

# Age at MRI scan
mri_df <- mri_df |>
  mutate(age_mri = (as.numeric(
    as.Date(date_mri1) - as.Date(calendar_date)
  ) / 365) + age_accel)

# create datadist
options(datadist = datadist(mri_df))

## function for bootstrapping substitution effects

make_ilr <- function(comp) {
  return(ilr(acomp(comp), V = v))
}

best_and_worst <- get_best_and_worst_comp(boot_df)
best_and_worst <- as.data.frame(apply(best_and_worst, 2, make_ilr))

run_mri_subs_bootstrap <- function(
    boot_data, comp_df, create_formula_fn, output_name) {
  # Matrix of variables to include in imputation model
  predmat <- quickpred(boot_data,
    mincor = 0,
    exclude = c(
      "eid",
      "calendar_date",
      "avg_sleep",
      "avg_inactivity",
      "avg_light",
      "avg_mvpa"
    )
  )

  ## reference compositions
  all_comp <- acomp(comp_df[
    , c("avg_sleep", "avg_inactivity", "avg_light", "avg_mvpa")
  ])
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
    statistic = bootstrap_mri_subs_fn,
    create_formula_fn = create_formula_fn,
    predmat = predmat,
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

run_mri_bootstrap <- function(boot_data, create_formula_fn, output_name) {
  # Matrix of variables to include in imputation model
  predmat <- quickpred(boot_data,
    mincor = 0,
    exclude = c(
      "eid",
      "calendar_date",
      "avg_sleep",
      "avg_inactivity",
      "avg_light",
      "avg_mvpa"
    )
  )

  result <- boot(
    data = boot_data,
    statistic = bootstrap_mri_fn,
    create_formula_fn = create_formula_fn,
    predmat = predmat,
    best_and_worst = best_and_worst,
    R = bootstrap_iterations,
    parallel = "multicore",
    ncpus = ncpus
  )

  # Prepend timestamp to avoid accidental data loss
  timestamp <- format(Sys.time(), "%Y-%m-%d_%H:%M")
  output_name_with_timestamp <- paste0(output_name, "_", timestamp, ".rds")
  saveRDS(result, file.path(output_dir, output_name_with_timestamp))
}

bootstrap_mri_fn <- function(
    data,
    indices,
    create_formula_fn,
    predmat,
    best_and_worst) {
  this_sample <- data[indices, ]

  imp <- mice(
    this_sample,
    m = 1, maxit = maxit, predictorMatrix = predmat
  )
  imp <- complete(imp)

  result_df <- predict_mri_results(imp, best_and_worst)

  return(as.matrix(result_df))
}

bootstrap_mri_subs_fn <- function(
    data,
    indices,
    create_formula_fn,
    predmat,
    short_sleep_geo_mean,
    avg_sleep_geo_mean) {
  this_sample <- data[indices, ]

  imp <- mice(
    this_sample,
    m = 1, maxit = maxit, predictorMatrix = predmat
  )
  imp <- complete(imp)

  result_df <- predict_all_substitutions(
    imp, short_sleep_geo_mean, avg_sleep_geo_mean
  )

  return(as.matrix(result_df))
}


process_boot_output <- function(directory, rds_path) {
  data <- readRDS(file.path(directory, rds_path))

  # tidy bootstrap output

  ## prepare data for plotting
  # Define the phenos and comps
  phenos <- rep(c("tbv", "wmv", "gmv", "hip", "log_wmh"), each = 3)
  comps <- rep(c("worst", "common", "best"), times = 5)

  col_names <- paste(comps, phenos, sep = "_")
  boot_reps <- as_tibble(data$t)
  colnames(boot_reps) <- col_names

  # Combine the phenos, comps, and values into a new data frame
  plot_data <- data.frame(comp = comps, pheno = phenos, value = data$t0)

  num_comps <- 3

  get_quantiles <- function(start, pheno) {
    slice <- data$t[, start:(start + num_comps - 1)]

    quantiles <-
      as.data.frame(t(apply(
        slice, 2, function(column) quantile(column, probs = c(0.025, 0.975))
      )))
    colnames(quantiles) <- c("lower", "upper")
    comps <- c("worst", "common", "best")
    quantiles$comp <- comps
    quantiles$pheno <- pheno
    return(quantiles)
  }

  pheno_start <- 1

  tbv_quants <- get_quantiles(pheno_start, "tbv")
  pheno_start <- pheno_start + num_comps

  wmv_quants <- get_quantiles(pheno_start, "wmv")
  pheno_start <- pheno_start + num_comps

  gmv_quants <- get_quantiles(pheno_start, "gmv")
  pheno_start <- pheno_start + num_comps

  hip_quants <- get_quantiles(pheno_start, "hip")
  pheno_start <- pheno_start + num_comps

  log_wmh_quants <- get_quantiles(pheno_start, "log_wmh")

  all_quantiles <- rbind(
    tbv_quants,
    wmv_quants,
    gmv_quants,
    hip_quants,
    log_wmh_quants
  )


  plot_data <- full_join(plot_data, all_quantiles, by = c("comp", "pheno"))

  return(list(plot_data, boot_reps))
}


### Plot

plot_mri <- function() {
  plot_data <- result_list[[1]]

  plot_data$comp <- str_to_title(plot_data$comp)
  plot_data$comp <- ifelse(
    plot_data$comp == "Common", "Average", plot_data$comp
  )
  plot_data$comp <- factor(
    plot_data$comp,
    levels = c("Best", "Average", "Worst")
  )

  single_plot <- function(pheno) {
    plot_data[plot_data$pheno == pheno, ] |>
      ggplot(aes(x = comp, y = value)) +
      geom_pointrange(aes(ymin = lower, ymax = upper),
        size = 0.1
      ) +
      labs(x = "Composition") +
      theme_cowplot()
  }

  tbv_plot <- single_plot("tbv")
  gmv_plot <- single_plot("gmv")
  wmv_plot <- single_plot("wmv")
  hip_plot <- single_plot("hip")
  wmh_plot <- single_plot("log_wmh")

  plot_grid(
    tbv_plot +
      labs(y = expression(paste("Total brain volume ", (cm^{
        3
      })))) +
      labs(x = "") +
      ylim(c(
        mean(mri_model_data$tbv, na.rm = T) - 0.5 *
          sd(mri_model_data$tbv, na.rm = T),
        mean(mri_model_data$tbv, na.rm = T) + 0.5 *
          sd(mri_model_data$tbv, na.rm = T)
      )),
    gmv_plot +
      labs(y = expression(paste("Grey matter volume ", (cm^{
        3
      })))) +
      labs(x = "") +
      ylim(c(
        mean(mri_model_data$gmv, na.rm = T) - 0.5 *
          sd(mri_model_data$gmv, na.rm = T),
        mean(mri_model_data$gmv, na.rm = T) + 0.5 *
          sd(mri_model_data$gmv, na.rm = T)
      )),
    wmv_plot +
      labs(
        y = expression(paste("White matter volume ", (cm^{
          3
        }))),
        x = ""
      ) +
      ylim(c(
        mean(mri_model_data$wmv, na.rm = T) - 0.5 *
          sd(mri_model_data$wmv, na.rm = T),
        mean(mri_model_data$wmv, na.rm = T) + 0.5 *
          sd(mri_model_data$wmv, na.rm = T)
      )),
    hip_plot +
      labs(y = expression(paste("Hippocampal volume ", (cm^{
        3
      })))) +
      ylim(c(
        mean(mri_model_data$hip, na.rm = T) - 0.5 *
          sd(mri_model_data$hip, na.rm = T),
        mean(mri_model_data$hip, na.rm = T) + 0.5 *
          sd(mri_model_data$hip, na.rm = T)
      )),
    wmh_plot +
      labs(y = "Log white matter hyperintensities") +
      ylim(c(
        mean(mri_model_data$log_wmh, na.rm = T) - 0.5 *
          sd(mri_model_data$log_wmh, na.rm = T),
        mean(mri_model_data$log_wmh, na.rm = T) + 0.5 *
          sd(mri_model_data$log_wmh, na.rm = T)
      )),
    nrow = 3
  )
}

# plot_mri()

# save plot

# ggsave(
#   file.path(
#     data_dir,
#     "../../Papers/Substitution Analysis/Main_figures/MRI_compositions.png"
#   ),
#   device = "png",
#   bg = "white",
#   width = 8,
#   height = 12,
#   dpi = 500
# )

### Estimated differences between compositions

# estimates <- result_list[[1]]
# boot_reps <- result_list[[2]]

get_contrasts <- function(pheno) {
  # estimates
  estimates <- estimates[estimates$pheno == pheno, ]
  estimates <- estimates |>
    select(pheno, comp, value) |>
    pivot_wider(names_from = comp, values_from = value)

  # CIs
  boot_reps <- boot_reps |> select(ends_with(pheno))
  names(boot_reps) <- str_remove(names(boot_reps), "_(.*)")

  out <-
    c(
      paste(
        pheno,
        ": worst - common ",
        round(estimates$worst - estimates$common, 2),
        " (",
        round(quantile(boot_reps$worst - boot_reps$common, 0.025), 2),
        ",",
        round(quantile(boot_reps$worst - boot_reps$common, 0.975), 2),
        ")",
        sep = ""
      ),
      paste(
        pheno,
        ": best - common ",
        round(estimates$best - estimates$common, 3),
        " (",
        round(quantile(boot_reps$best - boot_reps$common, 0.025), 3),
        ",",
        round(quantile(boot_reps$best - boot_reps$common, 0.975), 3),
        ")",
        sep = ""
      )
    )

  return(out)
}

process_boot_subs_output <- function(directory, rds_path) {
  data <- readRDS(file.path(directory, rds_path))

  # tidy bootstrap output

  num_subs <- sub_steps * 2 + 1
  sub_len <- num_subs * 5

  sub_col_names <- data$t[1, 1:num_subs]

  outcomes <- c("tbv", "wmv", "gmv", "hip", "log_wmh")
  num_outcomes <- length(outcomes)

  # Initialize indices
  t_ref_idx <- sub_len + 1
  t0_sub_idx <- 2
  sub_idx <- 1
  t0_outcome_idx <- 1

  # Initialize a list to store the plot data for each outcome
  get_plot_data <- function(
      t0_sub_idx, t_ref_idx, t0_outcome_idx, sub_idx) {
    whole_sample_values <- data$t0[
      t0_outcome_idx:(t0_outcome_idx + num_subs - 1), t0_sub_idx
    ]
    ref_values <- data$t[, t_ref_idx:(t_ref_idx + sub_len - 1)]
    sub_values <- ref_values[, sub_idx:(sub_idx + num_subs - 1)]

    return(data.frame(
      offset = sub_col_names,
      value = whole_sample_values,
      lower = apply(sub_values, 2, quantile, probs = 0.025),
      upper = apply(sub_values, 2, quantile, probs = 0.975)
    ))
  }

  # Initialize empty lists
  tbv_plot_list <- list()
  gmv_plot_list <- list()
  wmv_plot_list <- list()
  hip_plot_list <- list()
  log_wmh_plot_list <- list()

  # Loop over the outcomes
  for (reference in c("avg", "short")) {
    for (activity_level in c("inactive", "light", "mvpa")) {
      for (i in 1:num_outcomes) {
        # Construct the plot data
        plot_data <- get_plot_data(
          t0_sub_idx, t_ref_idx, t0_outcome_idx, sub_idx
        )

        # Store the plot data in the corresponding list
        plot_name <- paste0(reference, "_", activity_level)
        if (outcomes[i] == "tbv") tbv_plot_list[[plot_name]] <- plot_data
        if (outcomes[i] == "gmv") gmv_plot_list[[plot_name]] <- plot_data
        if (outcomes[i] == "wmv") wmv_plot_list[[plot_name]] <- plot_data
        if (outcomes[i] == "hip") hip_plot_list[[plot_name]] <- plot_data
        if (outcomes[i] == "log_wmh") {
          log_wmh_plot_list[[plot_name]] <- plot_data
        }

        # Update the indices for the next iteration
        t0_outcome_idx <- t0_outcome_idx + num_subs
        sub_idx <- sub_idx + num_subs
      }

      t_ref_idx <- t_ref_idx + sub_len
      t0_sub_idx <- t0_sub_idx + 1
      sub_idx <- 1
      t0_outcome_idx <- 1
    }
  }

  tbv_plot_df <- bind_rows(tbv_plot_list, .id = "Type")
  wmv_plot_df <- bind_rows(wmv_plot_list, .id = "Type")
  gmv_plot_df <- bind_rows(gmv_plot_list, .id = "Type")
  hip_plot_df <- bind_rows(hip_plot_list, .id = "Type")
  log_wmh_plot_df <- bind_rows(log_wmh_plot_list, .id = "Type")

  tbv_plot_df <- tbv_plot_df %>%
    separate(Type, into = c("Reference", "Substitution"), sep = "_")
  wmv_plot_df <- wmv_plot_df %>%
    separate(Type, into = c("Reference", "Substitution"), sep = "_")
  gmv_plot_df <- gmv_plot_df %>%
    separate(Type, into = c("Reference", "Substitution"), sep = "_")
  hip_plot_df <- hip_plot_df %>%
    separate(Type, into = c("Reference", "Substitution"), sep = "_")
  log_wmh_plot_df <- log_wmh_plot_df %>%
    separate(Type, into = c("Reference", "Substitution"), sep = "_")


  mri_plots <- function(outcome_df, outcome, ylabel) {
    get_plot <- function(outcome_df, outcome, sub, refcomp, colour, ylabel) {
      out_df <- outcome_df
      out_df$Substitution <-
        ifelse(out_df$Substitution == "inactive", "inactivity",
          ifelse(out_df$Substitution == "light",
            "light activity", "MVPA"
          )
        )

      out_df <- out_df[out_df$Substitution == sub &
        out_df$Reference == refcomp, ]

      # plot limits
      lower_lim <- mean(mri_model_data[[outcome]], na.rm = T) - 0.25 *
        sd(mri_model_data[[outcome]], na.rm = T)
      upper_lim <- mean(mri_model_data[[outcome]], na.rm = T) + 0.25 *
        sd(mri_model_data[[outcome]], na.rm = T)

      minutes_offset <- 0.2 * (upper_lim - lower_lim)
      sleep_offset <- 0.275 * (upper_lim - lower_lim)
      arrow_offset <- 0.325 * (upper_lim - lower_lim)
      sub_offset <- 0.375 * (upper_lim - lower_lim)

      out_df |>
        ggplot(aes(x = offset, y = value)) +
        geom_line(colour = colour) +
        geom_ribbon(aes(ymin = lower, ymax = upper),
          alpha = 0.2, fill = colour
        ) +
        xlab("") +
        ylab(ylabel) +
        facet_wrap(~Substitution, nrow = 2) +
        annotate(
          geom = "text", x = 0, y = lower_lim - minutes_offset,
          hjust = 0.5, fontface = 1, size = 14 / .pt,
          label = "Minutes", family = "serif"
        ) +
        annotate(
          geom = "text", x = -20, y = lower_lim - sleep_offset,
          hjust = 1, fontface = 1, size = 12 / .pt,
          label = "Less sleep", family = "serif"
        ) +
        annotate(
          geom = "text", x = 20, y = lower_lim - sleep_offset,
          hjust = 0, fontface = 1, size = 12 / .pt,
          label = "More sleep", family = "serif"
        ) +
        geom_segment(
          aes(
            x = 1, y = lower_lim - arrow_offset,
            xend = 15, yend = lower_lim - arrow_offset
          ),
          arrow = arrow(length = unit(0.15, "cm"))
        ) +
        geom_segment(
          aes(
            x = -1, y = lower_lim - arrow_offset,
            xend = -15, yend = lower_lim - arrow_offset
          ),
          arrow = arrow(length = unit(0.15, "cm"))
        ) +
        annotate(
          geom = "text", x = -20, y = lower_lim - sub_offset,
          hjust = 1, size = 12 / .pt,
          label = paste("More", out_df$Substitution),
          family = "serif", fontface = 1, size = 12 / .pt
        ) +
        annotate(
          geom = "text", x = 20, y = lower_lim - sub_offset,
          hjust = 0,
          label = paste("Less", out_df$Substitution),
          family = "serif", fontface = 1, size = 12 / .pt
        ) +
        coord_cartesian(
          ylim = c(lower_lim, upper_lim),
          expand = FALSE, clip = "off"
        ) +
        cowplot::theme_cowplot() +
        theme(
          text = element_text(size = 12, family = "serif"),
          plot.margin = unit(c(1, 1, 4, 1), "lines"),
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          axis.title.y = element_text(margin = margin(
            t = 0, r = 10, b = 0, l = 0
          ))
        )
    }

    # normal sleepers
    p1 <- get_plot(
      outcome_df, outcome, "inactivity", "avg", "#fc020f",
      ylabel = ylabel
    )
    p2 <- get_plot(
      outcome_df, outcome, "light activity", "avg", "#145e01",
      ylabel = ylabel
    )
    p3 <- get_plot(
      outcome_df, outcome, "MVPA", "avg", "#011869",
      ylabel = ylabel
    )

    pnorm <-
      plot_grid(NULL,
        p1 + labs(x = "", title = "A") + theme(
          legend.position = "none", plot.title.position = "plot",
          plot.title = element_text(size = 16)
        ),
        p2 + labs(x = "", title = "C") + theme(
          legend.position = "none", plot.title.position = "plot",
          plot.title = element_text(size = 16)
        ),
        p3 + labs(x = "", title = "D") + theme(
          legend.position = "none", plot.title.position = "plot",
          plot.title = element_text(size = 16)
        ),
        align = "vh",
        rel_heights = c(0.05, 1, 1, 1),
        nrow = 4,
        labels = "Normal sleepers",
        hjust = -1
      )

    # short sleepers
    p4 <- get_plot(
      outcome_df, outcome, "inactivity", "short", "#ff747b",
      ylabel = ylabel
    )
    p5 <- get_plot(
      outcome_df, outcome, "light activity", "short", "#6ed853",
      ylabel = ylabel
    )
    p6 <- get_plot(
      outcome_df, outcome, "MVPA", "short", "#708ff9",
      ylabel = ylabel
    )

    pshort <-
      plot_grid(NULL,
        p4 + labs(x = "", y = "", title = "B") + theme(
          legend.position = "none", plot.title.position = "plot",
          plot.title = element_text(size = 16)
        ),
        p5 + labs(x = "", y = "", title = "D") + theme(
          legend.position = "none", plot.title.position = "plot",
          plot.title = element_text(size = 16)
        ),
        p6 + labs(x = "", y = "", title = "F") + theme(
          legend.position = "none", plot.title.position = "plot",
          plot.title = element_text(size = 16)
        ),
        align = "vh",
        rel_heights = c(0.05, 1, 1, 1),
        nrow = 4,
        labels = "Short sleepers",
        hjust = -1
      )

    plot <- plot_grid(pnorm,
      pshort,
      nrow = 1
    )
    return(plot)
  }

  tbv_plot <- mri_plots(
    tbv_plot_df, "tbv",
    ylabel = expression(
      paste("Total brain volume ", cm^{
        3
      })
    )
  )
  wmv_plot <- mri_plots(
    wmv_plot_df, "wmv",
    ylabel = expression(
      paste("White matter volume ", cm^{
        3
      })
    )
  )
  gmv_plot <- mri_plots(
    gmv_plot_df, "gmv",
    ylabel = expression(
      paste("Grey matter volume ", cm^{
        3
      })
    )
  )
  hip_plot <- mri_plots(
    hip_plot_df, "hip",
    ylabel = expression(
      paste("Hippocampal volume ", cm^{
        3
      })
    )
  )
  wmh_plot <- mri_plots(log_wmh_plot_df, "log_wmh", ylabel = "Log WMH")

  return(list(
    tbv = tbv_plot, wmv = wmv_plot, gmv = gmv_plot,
    hip = hip_plot, wmh = wmh_plot
  ))
}


result_list <- process_boot_subs_output(data_dir, "boot_mri_subs.rds")

# save plots

for (i in names(result_list)) {
  ggsave(
    file.path(
      data_dir,
      paste("../../Papers/Substitution Analysis/Main_figures/Substitutions_",
        str_to_upper(i), ".png",
        sep = ""
      )
    ),
    plot = result_list[[i]],
    device = "png",
    bg = "white",
    width = 10,
    height = 12,
    dpi = 500
  )
}


#### Boot contrasts

get_boot_contrasts <- function(offset) {
  data <- readRDS(file.path(data_dir, "boot_mri_subs.rds"))

  # tidy bootstrap output

  num_subs <- sub_steps * 2 + 1
  sub_len <- num_subs * 5

  sub_col_names <- data$t[1, 1:num_subs]

  outcomes <- c("tbv", "wmv", "gmv", "hip", "log_wmh")
  num_outcomes <- length(outcomes)

  # Initialize indices
  t_ref_idx <- sub_len + 1
  t0_sub_idx <- 2
  sub_idx <- 1
  t0_outcome_idx <- 1

  get_contrast_data <- function(
      t0_sub_idx, t_ref_idx, t0_outcome_idx, sub_idx) {
    whole_sample_values <-
      data$t0[t0_outcome_idx:(t0_outcome_idx + num_subs - 1), t0_sub_idx]
    whole_sample_values <- cbind(
      whole_sample_values, sub_col_names
    ) |> as_tibble()


    # whole sample
    ref_value <- as.numeric(whole_sample_values[
      whole_sample_values$sub_col_names == 0, 1
    ])
    offset_value <- as.numeric(whole_sample_values[
      whole_sample_values$sub_col_names == offset, 1
    ])
    dif_full_sample <- offset_value - ref_value

    # bootstrap samples
    ref_values <- data$t[, t_ref_idx:(
      t_ref_idx + sub_len - 1)]
    sub_values <- ref_values[, sub_idx:(
      sub_idx + num_subs - 1)]
    dif_samples <- sub_values[
      , sub_col_names == offset
    ] - sub_values[, sub_col_names == 0]

    return(tibble(
      offset = offset,
      ref_value = ref_value,
      offset_value = offset_value,
      mean_dif = dif_full_sample,
      lower = quantile(dif_samples, probs = 0.025),
      upper = quantile(dif_samples, probs = 0.975),
    ))
  }

  # Initialize empty lists
  tbv_list <- list()
  gmv_list <- list()
  wmv_list <- list()
  hip_list <- list()
  log_wmh_list <- list()

  # Loop over the outcomes
  for (reference in c("avg", "short")) {
    for (activity_level in c("inactive", "light", "mvpa")) {
      for (i in 1:num_outcomes) {
        # Construct the contrast data
        contrast_data <- get_contrast_data(
          t0_sub_idx, t_ref_idx, t0_outcome_idx, sub_idx
        )

        # Store the contrast data in the corresponding list
        contrast_name <- paste0(reference, "_", activity_level)
        if (outcomes[i] == "tbv") tbv_list[[contrast_name]] <- contrast_data
        if (outcomes[i] == "gmv") gmv_list[[contrast_name]] <- contrast_data
        if (outcomes[i] == "wmv") wmv_list[[contrast_name]] <- contrast_data
        if (outcomes[i] == "hip") hip_list[[contrast_name]] <- contrast_data
        if (outcomes[i] == "log_wmh") {
          log_wmh_list[[contrast_name]] <- contrast_data
        }

        # Update the indices for the next iteration
        t0_outcome_idx <- t0_outcome_idx + num_subs
        sub_idx <- sub_idx + num_subs
      }

      t_ref_idx <- t_ref_idx + sub_len
      t0_sub_idx <- t0_sub_idx + 1
      sub_idx <- 1
      t0_outcome_idx <- 1
    }
  }

  tbv_df <- bind_rows(tbv_list, .id = "Type")
  wmv_df <- bind_rows(wmv_list, .id = "Type")
  gmv_df <- bind_rows(gmv_list, .id = "Type")
  hip_df <- bind_rows(hip_list, .id = "Type")
  log_wmh_df <- bind_rows(log_wmh_list, .id = "Type")

  return(bind_rows(list(
    tbv = tbv_df, wmv = wmv_df, gmv = gmv_df,
    hip = hip_df, wmh = log_wmh_df
  ), .id = "outcome"))
}


# dataframe of contrasts

map(c(-60, 60), get_boot_contrasts) |>
  bind_rows() |>
  filter(outcome == "hip" & str_detect(Type, "light"))
