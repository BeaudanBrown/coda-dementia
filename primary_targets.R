primary_targets <- list(
  ### MODELS ###
  tar_target(
    primary_models,
    fit_models(imp, timegroup_cuts, get_primary_formula),
    pattern = map(imp)
  ),
  ### REF RISK ###
  tar_target(
    primary_ref_risk,
    get_ref_risk(imp, primary_models, final_time),
    pattern = map(imp, primary_models)
  ),
  ### SUB RISK ###
  tar_target(
    primary_sub_risk,
    bind_rows(apply(all_subs, 1, function(sub) {
      get_sub_risk(
        imp,
        sub["from_var"],
        sub["to_var"],
        as.numeric(sub["duration"]),
        primary_models,
        final_time,
        comp_limits
      )
    })),
    pattern = map(imp, primary_models),
  ),
  ### COHORT ANALYSES ###
  tar_map(
    values = cohorts,
    names = cohort,
    tar_target(
      primary_ref_avg_risks,
      average_sub_results(primary_ref_risk, df, filter_fn)
    ),
    tar_target(
      primary_sub_avg_risks,
      average_sub_results(primary_sub_risk, df, filter_fn)
    ),
    tar_target(
      primary_risk_ratios,
      merge_risks(primary_sub_avg_risks, primary_ref_avg_risks, cohort)
    ),
    tar_target(
      primary_plots,
      make_plot(primary_risk_ratios, cohort)
    )
  ),
  ### PLOTS ###
  tar_target(
    all_primary_plots,
    rbind(
      primary_plots_short_sleeper,
      primary_plots_avg_sleeper
    )
  ),
  tar_target(
    primary_plot_grid,
    make_plot_grid(
      all_primary_plots,
      list(
        cohort_order = c("avg_sleeper", "short_sleeper"),
        color_order = c("#708ff9", "#ff747b"),
        subtype_order = c("avg_inactivity", "avg_light", "avg_mvpa")
      )
    )
  )
)
