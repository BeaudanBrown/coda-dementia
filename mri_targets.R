mri_targets <- list(
  tar_map(
    values = list(
      outcome = c("tbv", "wmv", "gmv", "hip", "log_wmh")
    ),
    names = outcome,
    ### MODELS ###
    tar_target(
      mri_models,
      get_mri_model(mri_imp, outcome),
      pattern = map(mri_imp),
      iteration = "list"
    ),
    ### REF ESTIMATE ###
    tar_target(
      mri_ref_results,
      get_mri_ref(mri_imp, outcome, mri_models),
      pattern = map(mri_imp, mri_models)
    ),
    ### SUB ESTIMATE ###
    tar_target(
      mri_sub_results,
      bind_rows(apply(all_subs, 1, function(sub) {
        get_mri_subs(
          mri_imp,
          outcome,
          mri_models,
          sub["from_var"],
          sub["to_var"],
          as.numeric(sub["duration"]),
          comp_limits
        )
      })),
      pattern = map(mri_imp, mri_models)
    ),
    tar_map(
      values = cohorts,
      names = cohort,
      ### FILTERED REF ESTIMATES ###
      tar_target(
        mri_ref_avg_estimate,
        average_sub_results(
          mri_ref_results,
          df,
          filter_fn,
          result_name = "estimate"
        )
      ),
      ### FILTERED SUB ESTIMATES ###
      tar_target(
        mri_sub_avg_estimate,
        average_sub_results(
          mri_sub_results,
          df,
          filter_fn,
          result_name = "estimate"
        )
      ),
      ### CONTRASTED ESTIMATES ###
      tar_target(
        mri_mean_diffs,
        merge_estimates(
          mri_sub_avg_estimate,
          mri_ref_avg_estimate,
          outcome,
          cohort
        )
      )
    ),
    tar_target(
      mri_plot_data,
      rbind(
        mri_mean_diffs_short_sleeper,
        mri_mean_diffs_avg_sleeper,
        mri_mean_diffs_full_cohort
      )
    ),
    tar_target(
      mri_plots,
      make_mri_plots(mri_plot_data)
    )
  ),
  ### PLOTS ###
  tar_target(
    all_mri_plots,
    rbind(
      mri_plots_log_wmh,
      mri_plots_tbv,
      mri_plots_hip,
      mri_plots_gmv,
      mri_plots_wmv
    )
  ),
  tar_map(
    values = list(
      outcome_name = c("tbv", "wmv", "gmv", "hip", "log_wmh")
    ),
    names = outcome_name,
    tar_target(
      mri_plot_grid,
      make_plot_grid(
        all_mri_plots[outcome == outcome_name, ],
        list(
          cohort_order = c("full_cohort", "avg_sleeper", "short_sleeper"),
          subtype_order = c("avg_inactivity", "avg_light", "avg_mvpa")
        )
      )
    )
  )
)
