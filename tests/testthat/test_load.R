testthat::context("load")

test_that("load and convert", {
    ### load
    mir_cnt_mat = load_miR_counts()
    testthat::expect_is(mir_cnt_mat, "matrix")

    mir_rpm_mat = load_miR_RPM()
    testthat::expect_is(mir_rpm_mat, "matrix")

    rna_cnt_mat = load_RNA_counts()
    testthat::expect_is(rna_cnt_mat, "matrix")

    rna_rpm_mat = load_RNA_RPM()
    testthat::expect_is(rna_rpm_mat, "matrix")

    ### to tidy data.table
    mir_cnt_dt = matrix_2_data.table(mir_cnt_mat)
    testthat::expect_is(mir_cnt_dt, "data.table")

    mir_rpm_dt = matrix_2_data.table(mir_rpm_mat)
    testthat::expect_is(mir_rpm_dt, "data.table")

    rna_cnt_dt = matrix_2_data.table(rna_cnt_mat)
    testthat::expect_is(rna_cnt_dt, "data.table")

    rna_rpm_dt = matrix_2_data.table(rna_rpm_mat)
    testthat::expect_is(rna_rpm_dt, "data.table")

    ### back to matrix
    mir_cnt_mat2 = data.table_2_matrix(mir_cnt_dt)
    testthat::expect_is(mir_cnt_mat2, "matrix")

    mir_rpm_mat2 = data.table_2_matrix(mir_rpm_dt)
    testthat::expect_is(mir_rpm_mat2, "matrix")

    rna_cnt_mat2 = data.table_2_matrix(rna_cnt_dt)
    testthat::expect_is(rna_cnt_mat2, "matrix")

    rna_rpm_mat2 = data.table_2_matrix(rna_rpm_dt)
    testthat::expect_is(rna_rpm_mat2, "matrix")

})

