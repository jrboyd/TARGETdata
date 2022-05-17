testthat::context("surival")
library(testthat)
library(TARGETdata)

out_dir = tempdir()
set_out_dir(out_dir)

suppressWarnings({
  mat_rpm = load_RNA_RPM()
  clin_dt = load_clinical_data()
  meta_dt = make_meta_dt(mat_rpm, clin_dt)
  meta_dt = meta_dt[Phase %in% c("Phase_II_Discovery", "Phase_II_Validation")]
  meta_dt[, rnk := frank(runif(.N), ties.method = "first"), .(Phase)]
  meta_dt = meta_dt[rnk <= 10]
  mat_rpm = filter_expression_to_valid(mat_rpm, meta_dt)
})


test_that("run_survival", {
  res1 = run_survival(meta_dt, group_var = "ik_status")
  expect_is(res1$pval, "numeric")
  expect_is(res1$plot, "ggplot")
})

test_that("surv_fun.z", {
  res2 = surv_fun.z(mat_rpm, meta_dt, goi = rownames(mat_rpm)[1])
  expect_is(res2$result$pval, "numeric")
  expect_is(res2$result$plot, "ggplot")
  expect_is(res2$plots, "ggplot")
  expect_is(res2$expression_plot, "ggplot")
  expect_is(res2$expression_data, "data.frame")
})


test_that("run_survival_scan.phaseII", {
  res3 = run_survival_scan.phaseII(goi_todo = rownames(mat_rpm)[1:5],
                                   rpm_mat = mat_rpm,
                                   meta_dt = meta_dt,
                                   analysis_name = "test",
                                   n_cores = 1)
  expect_is(res3, "data.frame")

  res5 = plot_survival_discovery_vs_validation(res3)
  expect_is(res5, "ggplot")
})



test_that("plot_survival_goi", {
  res4 = plot_survival_goi(rownames(mat_rpm)[1], mat_rpm, meta_dt)
  expect_is(res4$plot, "ggplot")
  expect_is(res4$plot_parts, "list")
  expect_is(res4$pvals, "data.frame")
  expect_equal(ncol(res4$pvals), 3)
  expect_equal(colnames(res4$pvals)[-1], c("discovery", "validation"))
})

unlink(out_dir, recursive = TRUE)
