testthat::context("meta")
library(testthat)
library(TARGETdata)
mat = load_miR_RPM()
clin_dt = load_clinical_data()

suppressWarnings({
  mat.cov = filter_expression_to_valid(mat, clin_dt)
  dt = convert_matrix_2_data.table(mat)
  dt.cov = filter_expression_to_valid(dt, clin_dt)
  meta_dt = make_meta_dt(mat, clin_dt)
  mat = filter_expression_to_valid(mat, meta_dt)

})

setdiff(
  meta_dt$sample_id,
  colnames(mat)
)

setdiff(
  colnames(mat),
  meta_dt$sample_id
)

exp_cn = c(
  "sample_code",
  "patient_id",
  "sample_id",
  "Phase",
  "vital_status",
  "days_to_last_follow_up",
  "days_to_death",
  "sample_type",
  "sample_type_short",
  "cluster_id",
  "ik_status"
)

test_that("make_meta_dt mat", {
  meta_dt = make_meta_dt(mat)
  testthat::expect_setequal(
    colnames(meta_dt),
    exp_cn
  )

  meta_dt.full = make_meta_dt.full(mat)
  testthat::expect_true(
    all(exp_cn %in% colnames(meta_dt.full))

  )
  testthat::expect_gt(
    ncol(meta_dt.full),
    length(exp_cn)
  )
})


test_that("make_meta_dt dt", {
  dt = convert_matrix_2_data.table(mat)
  meta_dt = make_meta_dt(dt)
  testthat::expect_setequal(
    colnames(meta_dt),
    exp_cn
  )

  meta_dt.full = make_meta_dt.full(dt)
  testthat::expect_true(
    all(exp_cn %in% colnames(meta_dt.full))

  )
  testthat::expect_gt(
    ncol(meta_dt.full),
    length(exp_cn)
  )
})

