test_that("tailfit is reproducible", {
  set.seed(123)
  x <- rweibull(100, shape = 2, scale = 10)
  x <- x[x <= 12]
  a <- tailfit(x, dists = "weibull", n_starts_bounds = 10)
  b <- tailfit(x, dists = "weibull", n_starts_bounds = 10)
  expect_equal(a$maxFnInf_best, b$maxFnInf_best)
  expect_equal(a$pct_complete_best, b$pct_complete_best)
  expect_equal(a$maxFnInf_bounds, b$maxFnInf_bounds)
})
