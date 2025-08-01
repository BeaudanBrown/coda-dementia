sensitivities <- list(
  create_formula_fn = rlang::syms(c(
    "get_s1_formula",
    "get_s2_formula"
  )),
  sens_name = c("s1", "s2"),
  colour = c("#ff747b", "#708ff9")
)

covar_sensitivity_targets <- list(
  ### MODELS ###
  tar_map(
    values = sensitivities,
    names = sens_name,
    tar_target(
      sensitivity_models,
      fit_models(half_imp, timegroup_cuts, create_formula_fn),
      pattern = map(half_imp)
    ),

    ### REF RISK ###
    tar_target(
      sensitivity_ref_risk,
      get_ref_risk(full_imp, sensitivity_models, final_time),
      pattern = map(half_imp, sensitivity_models)
    ),
    ### SUB RISK ###
    tar_target(
      sensitivity_sub_risk,
      bind_rows(apply(all_subs, 1, function(sub) {
        get_sub_risk(
          full_imp,
          sub["from_var"],
          sub["to_var"],
          as.numeric(sub["duration"]),
          sensitivity_models,
          final_time,
          comp_limits
        )
      })),
      pattern = map(half_imp, sensitivity_models)
    ),
    ### ANALYSIS ###
    tar_target(
      sensitivity_ref_avg_risks,
      average_sub_results(sensitivity_ref_risk, df, no_filter_fn)
    ),
    tar_target(
      sensitivity_sub_avg_risks,
      average_sub_results(sensitivity_sub_risk, df, no_filter_fn)
    ),
    tar_target(
      sensitivity_risk_ratios,
      merge_risks(
        sensitivity_sub_avg_risks,
        sensitivity_ref_avg_risks,
        sens_name
      )
    ),
    tar_target(
      sensitivity_plots,
      make_plot(sensitivity_risk_ratios, sens_name)
    )
  ),
  tar_target(
    all_sensitivity_plots,
    rbind(
      sensitivity_plots_s1,
      sensitivity_plots_s2,
      representative_plots,
      reverse_causation_plots
    )
  ),
  # ### PLOTS ###
  tar_target(
    sensitivity_plot_grid,
    make_plot_grid(
      all_sensitivity_plots,
      list(
        cohort_order = c("s1", "s2", "representative", "reverse_causation"),
        subtype_order = c("avg_inactivity", "avg_light", "avg_mvpa")
      )
    )
  )
)
