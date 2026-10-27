test_that("numerically saturated bounds are handled", {
  fit <- tailfit(c(1, 2, 4))
  expect_true(is.infinite(fit$maxFnInf_bounds["upper"]))
  expect_equal(unname(fit$pct_complete_bounds["lower"]), 0)
})
