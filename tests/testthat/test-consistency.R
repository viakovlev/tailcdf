test_that("total and completeness are consistent", {
  set.seed(123)
  x <- rweibull(100, shape = 2, scale = 10)
  x <- x[x <= 12]
  fit <- tailfit(x, dists = "weibull", n_starts_bounds = 10)
  expect_equal(fit$pct_complete_best,
               100 * length(x) / fit$maxFnInf_best)
  expect_lte(fit$maxFnInf_bounds["lower"],
             fit$maxFnInf_bounds["upper"])
})
