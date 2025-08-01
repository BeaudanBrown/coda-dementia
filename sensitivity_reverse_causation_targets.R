reverse_causation_targets <- list(
  ### MODELS ###
  tar_target(
    imp_rc,
    {
      half_imp[half_imp$time_to_dem > (3 * 365), ]
    },
    pattern = map(half_imp)
  ),
  tar_target(
    reverse_causation_models,
    fit_models(imp_rc, timegroup_cuts, get_primary_formula),
    pattern = map(imp_rc)
  ),
  ### REF RISK ###
  tar_target(
    reverse_causation_ref_risk,
    get_ref_risk(imp_rc, reverse_causation_models, final_time),
    pattern = map(imp_rc, reverse_causation_models)
  ),
  ### SUB RISK ###
  tar_target(
    reverse_causation_sub_risk,
    bind_rows(apply(all_subs, 1, function(sub) {
      get_sub_risk(
        imp_rc,
        sub["from_var"],
        sub["to_var"],
        as.numeric(sub["duration"]),
        reverse_causation_models,
        final_time,
        comp_limits
      )
    })),
    pattern = map(imp_rc, reverse_causation_models)
  ),
  ### ANALYSIS ###
  tar_target(
    reverse_causation_ref_avg_risks,
    average_sub_results(reverse_causation_ref_risk, df, no_filter_fn)
  ),
  tar_target(
    reverse_causation_sub_avg_risks,
    average_sub_results(reverse_causation_sub_risk, df, no_filter_fn)
  ),
  tar_target(
    reverse_causation_risk_ratios,
    merge_risks(
      reverse_causation_sub_avg_risks,
      reverse_causation_ref_avg_risks,
      "full_cohort"
    )
  ),
  tar_target(
    reverse_causation_plots,
    make_plot(reverse_causation_risk_ratios, "full_cohort")
  ),
  # ### PLOTS ###
  tar_target(
    reverse_causation_plots_grid,
    make_plot_grid(
      reverse_causation_plots,
      list(
        cohort_order = c("full_cohort"),
        color_order = c("#708ff9"),
        subtype_order = c("avg_inactivity", "avg_light", "avg_mvpa")
      )
    )
  )
)
