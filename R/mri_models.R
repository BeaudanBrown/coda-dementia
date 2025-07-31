get_mri_formula <- function(data) {
  knots_age_mri <- quantile(
    data[["age_mri"]],
    c(0.1, 0.5, 0.9)
  )
  knots_deprivation <- quantile(
    data[["townsend_deprivation_index"]],
    c(0.1, 0.5, 0.9)
  )
  knots_fruit_veg <- quantile(
    data[["fruit_veg"]],
    c(0.1, 0.5, 0.9)
  )
  knots_tfmri <- quantile(
    data[["mean_tfmri_headmot"]],
    c(0.1, 0.5, 0.9)
  )
  knots_head_scale <- quantile(
    data[["head_scale"]],
    c(0.1, 0.5, 0.9)
  )
  knots_lat_bpos <- quantile(
    data[["scan_lat_bpos"]],
    c(0.1, 0.5, 0.9)
  )
  knots_trans_bpos <- quantile(
    data[["scan_trans_bpos"]],
    c(0.1, 0.5, 0.9)
  )
  knots_long_bpos <- quantile(
    data[["scan_long_bpos"]],
    c(0.1, 0.5, 0.9)
  )

  primary_formula <- as.formula(
    ~ poly(R1, 2) *
      as.numeric(avg_total_household_income) +
      poly(R2, 2) * as.numeric(avg_total_household_income) +
      poly(R3, 2) * as.numeric(avg_total_household_income) +
      poly(R1, 2) * sex +
      poly(R2, 2) * sex +
      poly(R3, 2) * sex +
      poly(R1, 2) * retired +
      poly(R2, 2) * retired +
      poly(R3, 2) * retired +
      poly(R1, 2) * smok_status +
      poly(R2, 2) * smok_status +
      poly(R3, 2) * smok_status +
      poly(R1, 2) * rcs(age_mri, knots_age_mri) +
      poly(R2, 2) * rcs(age_mri, knots_age_mri) +
      poly(R3, 2) * rcs(age_mri, knots_age_mri) +
      as.numeric(apoe_e4) +
      highest_qual +
      rcs(fruit_veg, knots_fruit_veg) +
      alc_freq +
      shift +
      rcs(townsend_deprivation_index, knots_deprivation) +
      psych_meds +
      ethnicity +
      rcs(age_mri, knots_age_mri) * sex +
      rcs(age_mri, knots_age_mri) * as.numeric(apoe_e4) +
      rcs(mean_tfmri_headmot, knots_tfmri) +
      rcs(head_scale, knots_head_scale) +
      rcs(scan_lat_bpos, knots_lat_bpos) +
      rcs(scan_trans_bpos, knots_trans_bpos) +
      scan_tabpos +
      assessment_centre_mri1
  )

  return(primary_formula)
}
