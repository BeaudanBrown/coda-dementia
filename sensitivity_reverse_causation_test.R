reverse_causation_targets_test <- list(
  ### MODELS ###
  tar_target(
    imp_rc_test,
    {
      half_imp[half_imp$time_to_dem > (5 * 365), ]
    },
    pattern = slice(half_imp, index = c(1))
  ),
  tar_target(
    reverse_causation_models_test,
    fit_models(imp_rc_test, timegroup_cuts, get_primary_formula),
    pattern = map(imp_rc_test)
  ),
  ### REF RISK ###
  tar_target(
    reverse_causation_ref_risk_test,
    get_ref_risk(imp_rc_test, reverse_causation_models_test, final_time),
    pattern = map(imp_rc_test, reverse_causation_models_test)
  ),
  ### SUB RISK ###
  tar_target(
    reverse_causation_sub_risk_test,
    bind_rows(apply(all_subs, 1, function(sub) {
      get_sub_risk(
        imp_rc_test,
        sub["from_var"],
        sub["to_var"],
        as.numeric(sub["duration"]),
        reverse_causation_models_test,
        final_time,
        comp_limits
      )
    })),
    pattern = map(imp_rc_test, reverse_causation_models_test)
  ),
  ### ANALYSIS ###
  tar_map(
    values = cohorts,
    names = cohort,
    tar_target(
      reverse_causation_ref_avg_risks_test,
      average_sub_results(reverse_causation_ref_risk_test, df, filter_fn)
    ),
    tar_target(
      reverse_causation_sub_avg_risks_test,
      average_sub_results(reverse_causation_sub_risk_test, df, filter_fn)
    ),
    tar_target(
      reverse_causation_risk_ratios_test,
      merge_risks(
        reverse_causation_sub_avg_risks_test,
        reverse_causation_ref_avg_risks_test,
        "reverse_causation"
      )
    ),
    tar_target(
      reverse_causation_plots_test,
      make_plot(
        reverse_causation_risk_ratios_test,
        "reverse_causation",
        sleep_cohort
      )
    )
  ),
  tar_target(
    sensitivity_plot_grid_test,
    make_plot_grid(
      rbind(
        reverse_causation_plots_test_short_sleeper,
        reverse_causation_plots_test_avg_sleeper
      ),
      list(
        cohort_order = c(
          "short_sleeper",
          "long_sleeper"
        ),
        subtype_order = c("avg_inactivity", "avg_light", "avg_mvpa")
      )
    )
  )
)
