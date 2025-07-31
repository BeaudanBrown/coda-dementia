representative_targets <- list(
  ### MODELS ###
  tar_target(
    imp_rep,
    # shift covariates to match mean (continuous vars) or
    # probability (categorical vars)
    # of Schoeler et al pseudo-pop (see paper)
    full_imp[, `:=`(
      sex = sample(
        c("female", "male"),
        .N,
        replace = TRUE,
        prob = c(0.504, 0.496)
      ),
      retired = rbinom(.N, 1, prob = 0.193),
      avg_total_household_income = sample(
        c("<18", "18-30", "31-50", "52-100", ">100"),
        .N,
        replace = TRUE,
        prob = c(0.264, 0.141, 0.205, 0.145, 0.245)
      ),
      smok_status = sample(
        c("current", "former", "never"),
        .N,
        replace = TRUE,
        prob = c(0.208, 0.359, 0.433)
      )
    )]
  ),
  tar_target(
    representative_models,
    fit_models(full_imp, timegroup_cuts, get_primary_formula)
  ),
  ### REF RISK ###
  tar_target(
    representative_ref_risk,
    get_ref_risk(imp_rep, representative_models, final_time)
  ),
  ### SUB RISK ###
  tar_target(
    all_subs_target,
    all_subs
  ),
  tar_target(
    representative_sub_risk,
    {
      get_sub_risk(
        imp_rep,
        all_subs_target$from_var,
        all_subs_target$to_var,
        as.numeric(all_subs_target$duration),
        representative_models,
        final_time,
        comp_limits
      )
    },
    pattern = map(all_subs_target)
  ),
  ### ANALYSIS ###
  tar_target(
    representative_ref_avg_risks,
    average_sub_results(representative_ref_risk, df, no_filter_fn)
  ),
  tar_target(
    representative_sub_avg_risks,
    average_sub_results(representative_sub_risk, df, no_filter_fn)
  ),
  tar_target(
    representative_risk_ratios,
    merge_risks(
      representative_sub_avg_risks,
      representative_ref_avg_risks,
      "full_cohort"
    )
  ),
  tar_target(
    representative_plots,
    make_plot(representative_risk_ratios, "#6ed853")
  ),
  # ### PLOTS ###
  tar_target(
    representative_plots_grid,
    {
      cohort_order <- c("full_cohort")
      subtype_order <- c("avg_inactivity", "avg_light", "avg_mvpa")
      make_plot_grid(representative_plots, cohort_order, subtype_order)
    }
  )
)
