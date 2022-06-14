testthat::context("app_functions")
library(testthat)
library(TARGETdata)

parse_breaks = TARGETdata:::parse_breaks
capitalize = TARGETdata:::capitalize

test_that("parse_breaks", {
  expect_equal(parse_breaks("-1,1"), c(-1, 1))
  expect_equal(parse_breaks("-1,"), c(-1))
  expect_equal(parse_breaks("   -1  ,   1  "), c(-1, 1))
  expect_equal(parse_breaks("1,-1"), c(-1, 1))
  expect_equal(parse_breaks("-1,-1,1,1"), c(-1, 1))
  expect_equal(parse_breaks("0"), c(0))
  expect_equal(parse_breaks("-1.5,-.5,.5,1.5"), c(-1.5, -.5, .5, 1.5))

  expect_error(parse_breaks(","), regexp = "invalid breaks input")
  expect_error(parse_breaks(""), regexp = "invalid breaks input")
  expect_error(parse_breaks("-1.5,-.5,.5 1.5"), regexp = "invalid breaks input")
  expect_error(parse_breaks("-1,1,asdf"), regexp = "invalid breaks input")
})

test_that("capitalize", {
  expect_equal(capitalize("asdf"), "Asdf")
  expect_equal(capitalize("foo bar"), "Foo Bar")
  expect_equal(capitalize("ASDF"), "ASDF")
})
