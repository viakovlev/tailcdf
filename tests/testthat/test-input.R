test_that("invalid event times are rejected", {
  expect_error(tailfit(c(1, 2, NA, 4)))
  expect_error(tailfit(c(-1, 2, 3, 4)))
  expect_error(tailfit(c(0, 2, 3, 4)))
  expect_error(tailfit(c(1, 2, Inf, 4)))
})
