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
    ### SYNTH TARGETS ###
    tar_map(
      values = list(
        comp = c("Best", "Worst", "Typical")
      ),
      names = comp,
      tar_target(
        mri_synth_results,
        {
          mri_imp$R1 <- reference_comps[[comp]]$R1
          mri_imp$R2 <- reference_comps[[comp]]$R2
          mri_imp$R3 <- reference_comps[[comp]]$R3
          get_mri_ref(mri_imp, outcome, mri_models)
        },
        pattern = map(mri_imp, mri_models)
      ),
      tar_target(
        mri_synth_avg_estimate,
        {
          est_data <- average_sub_results(
            mri_synth_results,
            df,
            no_filter_fn,
            result_name = "estimate"
          )
          data.table(
            comp = comp,
            outcome = outcome,
            estimate = mean(est_data$results),
            lower = quantile(est_data$results, 0.025, names = FALSE),
            upper = quantile(est_data$results, 0.975, names = FALSE)
          )
        }
      )
    ),
    tar_target(
      mri_synth_plot_data,
      rbind(
        mri_synth_avg_estimate_Best,
        mri_synth_avg_estimate_Worst,
        mri_synth_avg_estimate_Typical
      )
    ),
    tar_target(
      mri_synth_plot,
      make_mri_synth_plot(mri_synth_plot_data)
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
      mri_long_plot_data,
      rbind(
        mri_mean_diffs_long_sleeper
      )
    ),
    tar_target(
      mri_long_plots,
      make_mri_plots(mri_long_plot_data)
    ),
    tar_target(
      mri_plot_data,
      rbind(
        mri_mean_diffs_short_sleeper,
        mri_mean_diffs_avg_sleeper
      )
    ),
    tar_target(
      mri_plots,
      make_mri_plots(mri_plot_data)
    ),
    tar_target(
      mri_tables,
      list(
        shorts = process_mri_table(mri_mean_diffs_short_sleeper),
        avgs = process_mri_table(mri_mean_diffs_avg_sleeper),
        longs = process_mri_table(mri_mean_diffs_long_sleeper)
      )
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
  tar_target(
    long_mri_plots,
    rbind(
      mri_long_plots_log_wmh,
      mri_long_plots_tbv,
      mri_long_plots_hip,
      mri_long_plots_gmv,
      mri_long_plots_wmv
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
          cohort_order = c("avg_sleeper", "short_sleeper"),
          subtype_order = c("avg_inactivity", "avg_light", "avg_mvpa")
        )
      )
    ),
    tar_target(
      mri_long_plot_grid,
      make_plot_grid(
        long_mri_plots[outcome == outcome_name, ],
        list(
          cohort_order = c("long_sleeper"),
          subtype_order = c("avg_inactivity", "avg_light", "avg_mvpa")
        )
      )
    )
  ),
  ### PLOTS ###
  tar_target(
    mri_synth_plots,
    {
      plots <- list(
        mri_synth_plot_tbv,
        mri_synth_plot_gmv,
        mri_synth_plot_wmv,
        mri_synth_plot_hip,
        mri_synth_plot_log_wmh
      )
      patchwork::wrap_plots(c(plots), nrow = 3, ncol = 2)
    }
  )
)
