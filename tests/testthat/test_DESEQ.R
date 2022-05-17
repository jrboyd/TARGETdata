testthat::context("DESeq")
library(testthat)
library(TARGETdata)

out_dir = tempdir()
set_out_dir(out_dir)

suppressWarnings({
  mat_cnt = load_RNA_counts()
  clin_dt = load_clinical_data()
  meta_dt = make_meta_dt(mat_cnt, clin_dt)
  meta_dt = meta_dt[Phase %in% c("Phase_II_Discovery", "Phase_II_Validation")]
  meta_dt[, rnk := frank(runif(.N), ties.method = "first"), .(Phase)]
  meta_dt = meta_dt[rnk <= 10]
  mat_cnt = filter_expression_to_valid(mat_cnt, meta_dt)
})


test_that("run_DESeq", {
  des = run_DESeq(mat_cnt[1:100,], meta_dt)
  expect_is(des, "DESeqDataSet")

  write_out = write_DESeq_results(des, "test")
  expect_is(write_out, "list")
  expect_setequal(names(write_out), c("full", "significant"))
  expect_true(all(file.exists(write_out$full)))
  expect_true(all(file.exists(write_out$significant)))
})

unlink(out_dir, recursive = TRUE)
