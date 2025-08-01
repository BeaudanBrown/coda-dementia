representative_targets <- list(
  ### MODELS ###
  tar_target(
    imp_rep,
    # shift covariates to match mean (continuous vars) or
    # probability (categorical vars)
    # of Schoeler et al pseudo-pop (see paper)
    half_imp |>
      mutate(
        sex = sample(
          c("female", "male"),
          n(),
          replace = TRUE,
          prob = c(0.504, 0.496)
        ),
        sex = as.factor(sex),
        retired = rbinom(n(), 1, prob = 0.193),
        avg_total_household_income = sample(
          c("<18", "18-30", "31-50", "52-100", ">100"),
          n(),
          replace = TRUE,
          prob = c(0.264, 0.141, 0.205, 0.145, 0.245)
        ),
        avg_total_household_income = as.factor(avg_total_household_income),
        smok_status = sample(
          c("current", "former", "never"),
          n(),
          replace = TRUE,
          prob = c(0.208, 0.359, 0.433)
        ),
        smok_status = as.factor(smok_status),
      ),
    pattern = map(half_imp)
  ),
  tar_target(
    representative_models,
    fit_models(imp_rep, timegroup_cuts, get_primary_formula),
    pattern = map(imp_rep)
  ),
  ### REF RISK ###
  tar_target(
    representative_ref_risk,
    get_ref_risk(imp_rep, representative_models, final_time),
    pattern = map(imp_rep, representative_models)
  ),
  ### SUB RISK ###
  tar_target(
    all_subs_target,
    all_subs
  ),
  tar_target(
    representative_sub_risk,
    bind_rows(apply(all_subs, 1, function(sub) {
      get_sub_risk(
        imp_rep,
        sub["from_var"],
        sub["to_var"],
        as.numeric(sub["duration"]),
        representative_models,
        final_time,
        comp_limits
      )
    }))
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
    make_plot(representative_risk_ratios, "full_cohort")
  ),
  # ### PLOTS ###
  tar_target(
    representative_plots_grid,
    make_plot_grid(
      representative_plots,
      list(
        cohort_order = c("full_cohort"),
        color_order = c("#708ff9"),
        subtype_order = c("avg_inactivity", "avg_light", "avg_mvpa")
      )
    )
  )
)
