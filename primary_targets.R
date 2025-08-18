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
  ### RETIRED ANALYSES ###
  tar_map(
    values = retired_cohorts,
    names = cohort,
    tar_target(
      retired_ref_avg_risks,
      average_sub_results(primary_ref_risk, df, filter_fn)
    ),
    tar_target(
      retired_sub_avg_risks,
      average_sub_results(primary_sub_risk, df, filter_fn)
    ),
    tar_target(
      retired_risk_ratios,
      merge_risks(retired_sub_avg_risks, retired_ref_avg_risks, cohort)
    ),
    tar_target(
      primary_plots,
      make_plot(retired_risk_ratios, cohort)
    )
  ),
  ### PLOTS ###
  tar_target(
    all_primary_plots,
    rbind(
      primary_plots_avg_sleeper,
      primary_plots_short_sleeper
    )
  ),
  tar_target(
    retired_short_plot_grid,
    make_plot_grid(
      rbind(
        primary_plots_retired_short,
        primary_plots_not_retired_short
      ),
      list(
        cohort_order = c("retired", "not_retired"),
        subtype_order = c("avg_inactivity", "avg_light", "avg_mvpa")
      )
    )
  ),
  tar_target(
    retired_avg_plot_grid,
    make_plot_grid(
      rbind(
        primary_plots_retired_avg,
        primary_plots_not_retired_avg
      ),
      list(
        cohort_order = c("retired", "not_retired"),
        subtype_order = c("avg_inactivity", "avg_light", "avg_mvpa")
      )
    )
  ),
  tar_target(
    retired_long_plot_grid,
    make_plot_grid(
      rbind(
        primary_plots_retired_long,
        primary_plots_not_retired_long
      ),
      list(
        cohort_order = c("retired", "not_retired"),
        subtype_order = c("avg_inactivity", "avg_light", "avg_mvpa")
      )
    )
  ),
  tar_target(
    primary_plot_grid,
    make_plot_grid(
      all_primary_plots,
      list(
        cohort_order = c("avg_sleeper", "short_sleeper"),
        subtype_order = c("avg_inactivity", "avg_light", "avg_mvpa")
      )
    )
  )
)
