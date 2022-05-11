testthat::context("meta")

mat = load_miR_RPM()
exp_cn = c(
  "sample_code",
  "patient_id",
  "sample_id",
  "Phase",
  "vital_status",
  "days_to_last_follow_up",
  "days_to_death",
  "sample_type",
  "sample_type_short"
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
  dt = matrix_2_data.table(mat)
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

