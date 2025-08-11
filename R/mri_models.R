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
  knots_R1 <- quantile(
    data[["R1"]],
    c(0.1, 0.5, 0.9)
  )
  knots_R2 <- quantile(
    data[["R2"]],
    c(0.1, 0.5, 0.9)
  )
  knots_R3 <- quantile(
    data[["R3"]],
    c(0.1, 0.5, 0.9)
  )

  primary_formula <- as.formula(
    ~ rcs(R1, knots_R1) *
      as.numeric(avg_total_household_income) +
      rcs(R2, knots_R2) * as.numeric(avg_total_household_income) +
      rcs(R3, knots_R3) * as.numeric(avg_total_household_income) +
      rcs(R1, knots_R1) * sex +
      rcs(R2, knots_R2) * sex +
      rcs(R3, knots_R3) * sex +
      rcs(R1, knots_R1) * retired +
      rcs(R2, knots_R2) * retired +
      rcs(R3, knots_R3) * retired +
      rcs(R1, knots_R1) * smok_status +
      rcs(R2, knots_R2) * smok_status +
      rcs(R3, knots_R3) * smok_status +
      rcs(R1, knots_R1) * rcs(age_mri, knots_age_mri) +
      rcs(R2, knots_R2) * rcs(age_mri, knots_age_mri) +
      rcs(R3, knots_R3) * rcs(age_mri, knots_age_mri) +
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
