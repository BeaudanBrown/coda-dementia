source("mri_models.R")

mri_df <-
  fread(file.path(data_dir, "../MRI/mri_full_trimmed_v3.csv"), stringsAsFactors = TRUE) |>
  as_tibble()
mri_df <- mri_df |> rename("insomnia_scale_sr" = "insomnia_sr")
mri_df$assessment_centre_mri1 <- as.factor(mri_df$assessment_centre_mri1)
mri_df <- mri_df |> select(-starts_with("dur_day_total_"))

# MRI QC data
#########################################################################################
mri_qc_df <-
  fread(file.path(data_dir, "../MRI/mri_qc.csv"))

mri_qc_df <- mri_qc_df |>
  select(eid, ends_with("-2.0")) |>
  as_tibble()

mri_qc_df <- mri_qc_df |>
  rename(
    "discrep_t1_stand_lin" = `25731-2.0`,
    "discrep_t1_stand_nonlin" = `25732-2.0`,
    "warping_t1" = `25733-2.0`,
    "inv_sig_noise_t1" = `25734-2.0`,
    "inv_contr_noise_t1" = `25735-2.0`,
    "discrep_t2_t1" = `25736-2.0`,
    "discrep_swi_t1" = `25738-2.0`,
    "discrep_rfmri_t1" = `25739-2.0`,
    "discrep_tfmri_t1" = `25740-2.0`,
    "mean_tfmri_headmot" = `25742-2.0`,
    "inv_temp_signoise_pp_rfmri" = `25743-2.0`,
    "inv_temp_signoise_art_rfmri" = `25744-2.0`,
    "num_dmri_outslices" = `25746-2.0`,
    "scan_lat_bpos" = `25756-2.0`,
    "scan_trans_bpos" = `25757-2.0`,
    "scan_long_bpos" = `25758-2.0`,
    "scan_tabpos" = `25759-2.0`,
    "acq_prot_phase" = `25780-2.0`
  )

mri_qc_df <- mri_qc_df |> select(-starts_with("2"))

boot_df <- read_rds(file.path(data_dir, "bootstrap_data.rds"))

mri_df <- mri_df |>
  select(-any_of(names(boot_df)[names(boot_df) != "eid"])) |>
  left_join(boot_df, by = "eid")

mri_df <- mri_df |>
  left_join(mri_qc_df |> select(-any_of(names(mri_df)[names(boot_df) != "eid"])), by = "eid")

# Age at MRI scan
mri_df <- mri_df %>%
  mutate(
    date_mri1 = parse_date_time(date_mri1, "ymd"),
    date_birth = paste(year_birth, month_birth, sep = "-")
  ) %>%
  mutate(date_birth = parse_date_time(date_birth, "ym")) %>%
  mutate(age_mri = as.numeric(date_mri1 - date_birth) / 365)

## Create derived variables
mri_df$tbv <-
  (mri_df$vol_grey_matter_norm + mri_df$vol_white_matter_norm) / 1000
mri_df$wmv <-
  mri_df$vol_white_matter_norm / 1000
mri_df$gmv <-
  mri_df$vol_grey_matter_norm / 1000
mri_df$hip <-
  (mri_df$vol_hippocampus_l + mri_df$vol_hippocampus_r) *
  mri_df$volumetric_scaling_from_t1_head_image_to_standard_space / 1000
mri_df$tot_wmh <- mri_df$total_vol_white_matter_hyperintensities_from_t1_and_t2_flair_images
mri_df$log_wmh <-
  log(mri_df$tot_wmh * mri_df$volumetric_scaling_from_t1_head_image_to_standard_space)

# Model data for mri model
mri_model_data <- select(
  mri_df,
  R1, R2, R3,
  avg_sleep, avg_inactivity, avg_light, avg_mvpa,
  assessment_centre_mri1,
  tbv, wmv, gmv, hip, log_wmh,
  mean_tfmri_headmot,
  scan_lat_bpos,
  scan_trans_bpos,
  scan_long_bpos,
  scan_tabpos,
  smok_status,
  dem,
  time_to_dem,
  bp_syst_avg,
  avg_sri,
  age_assessment,
  bp_med,
  any_cvd,
  highest_qual,
  apoe_e4,
  BMI,
  antidepressant_med,
  antipsychotic_med,
  avg_WASO,
  retired,
  shift,
  pa_modvig,
  insomnia_med,
  ethnicity,
  age_mri,
  sex,
  parent_dementa,
  avg_total_household_income,
  townsend_deprivation_index,
  vol_white_matter_norm,
  vol_grey_matter_norm,
  vol_hippocampus_l,
  vol_hippocampus_r,
  volumetric_scaling_from_t1_head_image_to_standard_space,
  mri_accel_time_dif,
  tot_wmh,
  dem,
  avg_sleep_duration,
  time_to_dem,
  overall_health_rating
)

options(datadist = datadist(mri_model_data))

run_mri_bootstrap <- function(boot_data, create_formula_fn, output_name) {
  # Matrix of variables to include in imputation model
  predmat <- quickpred(boot_data,
    mincor = 0,
    exclude = c(
      "date_acdem2",
      "date_accel",
      "date_of_death",
      "avg_sleep",
      "avg_inactivity",
      "avg_light",
      "avg_mvpa"
    )
  )

  # method for each imputed variable
  imp_methods <- make.method(boot_data)
  # exclude dates from being imputed
  imp_methods["date_acdem2"] <- ""
  imp_methods["date_of_death"] <- ""

  make_ilr <- function(comp) {
    return(ilr(acomp(comp), V = v))
  }

  print("Calculating best and worst")
  print(format(Sys.time(), "%H:%M:%S"))
  best_and_worst <- get_best_and_worst_comp()
  best_and_worst <- as.data.frame(apply(best_and_worst, 2, make_ilr))

  result <- boot(
    data = boot_data,
    statistic = bootstrap_mri_fn,
    create_formula_fn = create_formula_fn,
    predmat = predmat,
    imp_methods = imp_methods,
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
  imp_methods,
  best_and_worst
) {
  this_sample <- data[indices, ]

  imp <- mice(this_sample, m = 1, maxit = 1, predictorMatrix = predmat, methods = imp_methods)
  imp <- complete(imp)

  result_df <- predict_mri_results(imp, best_and_worst)

  return(as.matrix(result_df))
}

process_boot_output <- function(directory, rds_path) {
  
  data <- readRDS(file.path(directory, rds_path))
  
  # tidy bootstrap output 
  
  col_names <- paste(comps, phenos, sep = "_")
  boot_reps <- as_tibble(data$t)
  colnames(boot_reps) <- col_names

  ## prepare data for plotting 
  # Define the phenos and comps
  phenos <- rep(c("tbv", "wmv", "gmv", "hip", "log_wmh"), each = 3)
  comps <- rep(c("worst", "common", "best"), times = 5)

  # Combine the phenos, comps, and values into a new data frame
  plot_data <- data.frame(comp = comps, pheno = phenos, value = data$t0)

  num_comps <- 3

  get_quantiles <- function(start, pheno) {
    slice <- data$t[, start:(start + num_comps - 1)]

    quantiles <-
      as.data.frame(t(apply(slice, 2, function(column) quantile(column, probs = c(0.025, 0.975)))))
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

  all_quantiles <- rbind(tbv_quants,
                         wmv_quants,
                         gmv_quants,
                         hip_quants,
                         log_wmh_quants)


  plot_data <- full_join(plot_data, all_quantiles, by = c("comp", "pheno"))
  
  return(list(plot_data, boot_reps))
}


# run_mri_bootstrap(
#   boot_data = mri_model_data,
#   create_formula_fn = get_mri_formula,
#   output_name = "boot_mri"
# )

result_list <- process_boot_output(data_dir, "boot_mri.rds")

### Plot 

plot_mri <- function(){
  
  plot_data <- result_list[[1]]
  
  plot_data$comp <- str_to_title(plot_data$comp)
  plot_data$comp <- ifelse(plot_data$comp == "Common", "Average", plot_data$comp)
  plot_data$comp <- factor(plot_data$comp, levels = c("Best","Average","Worst"))
  
  single_plot <- function(pheno){
    plot_data[plot_data$pheno == pheno,] |> 
      ggplot(aes(x = comp, y = value)) +
      geom_pointrange(aes(ymin = lower, ymax = upper),
                      size = 0.1) +
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
      labs(y = expression(paste("Total brain volume ", (cm^{3})))) +
      labs(x = "") +
      ylim(c(mean(mri_model_data$tbv,na.rm = T) - 0.5*sd(mri_model_data$tbv,na.rm = T),
             mean(mri_model_data$tbv,na.rm = T) + 0.5*sd(mri_model_data$tbv,na.rm = T))),
    gmv_plot + 
      labs(y = expression(paste("Grey matter volume ", (cm^{3})))) +
      labs(x = "") +
      ylim(c(mean(mri_model_data$gmv,na.rm = T) - 0.5*sd(mri_model_data$gmv,na.rm = T),
             mean(mri_model_data$gmv,na.rm = T) + 0.5*sd(mri_model_data$gmv,na.rm = T))),
    wmv_plot + 
      labs(y = expression(paste("White matter volume ", (cm^{3}))),
           x = "") +
      ylim(c(mean(mri_model_data$wmv,na.rm = T) - 0.5*sd(mri_model_data$wmv,na.rm = T),
             mean(mri_model_data$wmv,na.rm = T) + 0.5*sd(mri_model_data$wmv,na.rm = T))),
    hip_plot + 
      labs(y = expression(paste("Hippocampal volume ", (cm^{3})))) +
      ylim(c(mean(mri_model_data$hip,na.rm = T) - 0.5*sd(mri_model_data$hip,na.rm = T),
             mean(mri_model_data$hip,na.rm = T) + 0.5*sd(mri_model_data$hip,na.rm = T))),
    wmh_plot + 
      labs(y = "Log white matter hyperintensities") +
      ylim(c(mean(mri_model_data$log_wmh,na.rm = T) - 0.5*sd(mri_model_data$log_wmh,na.rm = T),
             mean(mri_model_data$log_wmh,na.rm = T) + 0.5*sd(mri_model_data$log_wmh,na.rm = T))),
    nrow = 3
  )  
}

#plot_mri()

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

estimates <- result_list[[1]]
boot_reps <- result_list[[2]]

get_contrasts <- function(pheno){
  
  # estimates
  estimates <- estimates[estimates$pheno == pheno,]
  estimates <- estimates |> select(pheno, comp, value) |> 
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

get_contrasts("wmv")

