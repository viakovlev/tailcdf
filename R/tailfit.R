# tailfit.R
# Final MLE-only tail extrapolation for right-truncated cumulative event data.
# The fitted family, parameters, Max/Ninf, completeness, and bounds are all based on truncated MLE.
# CDF values are calculated only from the MLE-fitted distribution for plotting and diagnostics.

#' Fit an MLE tail model to right-truncated event times
#'
#' @param x Numeric vector of strictly positive event times.
#' @param dists Candidate distributions. Default: lnorm, gamma, weibull, chisq, f.
#' @param envelope_expansion Bound width multiplier. 1.0 is the standard E100 bound.
#' @param n_starts Multistarts for best-fit optimization.
#' @param n_starts_bounds Multistarts for bound optimization.
#' @param seed Random seed.
#' @param pcts Probabilities for reporting times to fractions of total expected events.
#' @param ... Ignored compatibility arguments.
#' @return Object of class tailfit.
#' @examples
#' set.seed(123)
#' x <- rweibull(100, shape = 2, scale = 10)
#' x <- x[x <= 12]
#' fit <- tailfit(x, dists = 'weibull', n_starts_bounds = 10)
#' fit$maxFnInf_best
#' fit$pct_complete_best
#' @export
tailfit <- function(
  x,
  dists = c("lnorm", "gamma", "weibull", "chisq", "f"),
  envelope_expansion = 1.0,
  n_starts = 10,
  n_starts_bounds = max(80, 8 * n_starts),
  seed = 1,
  pcts = c(0.6, 0.7, 0.8, 0.9),
  ...
) {

  if (any(!is.finite(x) | x <= 0)) {
  stop("x must contain only positive finite values.", call. = FALSE)
}

  x <- sort(as.numeric(x))
  x <- x[is.finite(x) & x > 0]

  if (length(x) < 3L) {
    message("No fit possible: x must contain at least 3 positive finite values.")
    return(invisible(NULL))
  }

  n <- length(x)
  Xt <- max(x)
  y_emp <- seq_len(n)

  E <- sort(unique(as.numeric(envelope_expansion)))
  E <- E[is.finite(E) & E >= 0]
  if (!length(E)) E <- 1.0
  names_E <- paste0("E", formatC(round(100 * E), width = 2, flag = "0"))

  pcts <- sort(unique(as.numeric(pcts)))
  pcts <- pcts[is.finite(pcts) & pcts > 0 & pcts < 1]

  # ----------------------------
  # Distribution specifications
  # ----------------------------
  make_specs <- function(xvec) {
    list(
      lnorm = list(
        name = "lnorm",
        pfun = stats::plnorm,
        dfun = stats::dlnorm,
        qfun = stats::qlnorm,
        theta_to_par = function(th) c(meanlog = th[1], sdlog = exp(th[2])),
        start = {
          lx <- log(xvec)
          c(meanlog = mean(lx), lsdlog = log(max(stats::sd(lx), 1e-8)))
        }
      ),
      gamma = list(
        name = "gamma",
        pfun = stats::pgamma,
        dfun = stats::dgamma,
        qfun = stats::qgamma,
        theta_to_par = function(th) c(shape = exp(th[1]), rate = exp(th[2])),
        start = {
          m <- mean(xvec)
          v <- stats::var(xvec)
          if (!is.finite(v) || v <= 0) v <- m^2
          c(lshape = log(max(m^2 / v, 1e-6)), lrate = log(max(m / v, 1e-6)))
        }
      ),
      weibull = list(
        name = "weibull",
        pfun = stats::pweibull,
        dfun = stats::dweibull,
        qfun = stats::qweibull,
        theta_to_par = function(th) c(shape = exp(th[1]), scale = exp(th[2])),
        start = c(lshape = log(2), lscale = log(mean(xvec)))
      ),
      chisq = list(
        name = "chisq",
        pfun = stats::pchisq,
        dfun = stats::dchisq,
        qfun = stats::qchisq,
        theta_to_par = function(th) c(df = exp(th[1])),
        start = c(ldf = log(max(mean(xvec), 1e-6)))
      ),
      f = list(
        name = "f",
        pfun = stats::pf,
        dfun = stats::df,
        qfun = stats::qf,
        theta_to_par = function(th) c(df1 = exp(th[1]), df2 = exp(th[2])),
        start = c(ldf1 = log(5), ldf2 = log(10))
      )
    )
  }

  specs <- make_specs(x)
  dists <- intersect(dists, names(specs))
  if (!length(dists)) stop("No valid distributions in dists.", call. = FALSE)

  safe_par <- function(spec, theta) {
    par <- tryCatch(spec$theta_to_par(theta), error = function(e) NULL)
    if (is.null(par) || any(!is.finite(par))) return(NULL)
    par
  }

  safe_p <- function(spec, q, par) {
    tryCatch(do.call(spec$pfun, c(list(q = q), as.list(par))), error = function(e) NA_real_)
  }

  safe_d <- function(spec, xx, par, log = FALSE) {
    tryCatch(
      do.call(spec$dfun, c(list(x = xx), as.list(par), list(log = log))),
      error = function(e) rep(NA_real_, length(xx))
    )
  }

  safe_q <- function(spec, pp, par) {
    tryCatch(
      do.call(spec$qfun, c(list(p = pp), as.list(par))),
      error = function(e) rep(NA_real_, length(pp))
    )
  }

  make_starts <- function(theta0, n_start, seed_offset = 0L, wide = FALSE) {
    set.seed(seed + seed_offset)
    k <- length(theta0)
    starts <- matrix(NA_real_, nrow = n_start, ncol = k)
    starts[1, ] <- theta0

    if (n_start >= 2) {
      for (i in 2:n_start) {
        sdj <- if (wide) sample(c(0.4, 0.8, 1.5, 2.5), 1) else sample(c(0.25, 0.5, 0.9), 1)
        starts[i, ] <- theta0 + stats::rnorm(k, 0, sdj)
      }
    }

    starts
  }

  optimize_fit <- function(spec, objective_fun, starts, maxit = 6000) {
    best_val <- Inf
    best_theta <- NULL

    for (i in seq_len(nrow(starts))) {
      th0 <- as.numeric(starts[i, ])
      val0 <- objective_fun(th0)
      if (!is.finite(val0)) next

      # stats::optim(method = "Nelder-Mead") gives warnings and is unreliable
      # for one-parameter families such as chisq. Use one-dimensional Brent
      # optimization in that case, and Nelder-Mead only for multi-parameter
      # families.
      if (length(th0) == 1L) {
        interval_width <- 8
        opt <- try(
          stats::optimize(
            f = function(z) objective_fun(c(z)),
            interval = c(th0 - interval_width, th0 + interval_width),
            tol = 1e-10
          ),
          silent = TRUE
        )

        if (inherits(opt, "try-error") || !is.finite(opt$objective)) next

        val <- opt$objective
        theta <- c(opt$minimum)
      } else {
        opt <- try(
          stats::optim(
            th0,
            objective_fun,
            method = "Nelder-Mead",
            control = list(maxit = maxit, reltol = 1e-12)
          ),
          silent = TRUE
        )

        if (inherits(opt, "try-error") || !is.finite(opt$value)) next

        val <- opt$value
        theta <- opt$par
      }

      if (val < best_val) {
        best_val <- val
        best_theta <- theta
      }
    }

    if (is.null(best_theta)) return(NULL)
    list(theta = best_theta, value = best_val)
  }

  # ----------------------------
  # Truncated MLE fit
  # ----------------------------
  mle_stats <- function(spec, theta) {
    par <- safe_par(spec, theta)
    if (is.null(par)) return(NULL)

    Ft <- safe_p(spec, Xt, par)
    if (!is.finite(Ft) || Ft <= 0 || Ft >= 1) return(NULL)

    Fv <- safe_p(spec, x, par)
    if (length(Fv) != n || any(!is.finite(Fv)) || any(Fv < 0) || any(Fv > 1)) return(NULL)

    logfx <- safe_d(spec, x, par, log = TRUE)
    if (length(logfx) != n || any(!is.finite(logfx))) return(NULL)

    ll_i <- logfx - log(Ft)
    sum_ll <- sum(ll_i)
    mean_ll <- mean(ll_i)
    Ninf <- n / Ft
    pred <- Ninf * Fv
    absD <- abs(y_emp - pred)

    if (!is.finite(sum_ll) || !is.finite(Ninf) || any(!is.finite(absD))) return(NULL)

    list(
      par = par,
      Ft = Ft,
      Fv = Fv,
      logfx = logfx,
      ll_i = ll_i,
      sum_ll = sum_ll,
      mean_ll = mean_ll,
      Ninf = Ninf,
      comp_pct = 100 * n / Ninf,
      pred = pred,
      absD = absD,
      cdf_sum_absD = sum(absD),
      cdf_mean_absD = mean(absD),
      cdf_median_absD = stats::median(absD)
    )
  }

  # Excel-style spacing-to-log noise at a parameter vector.
  # Uses rows i=2..n:
  # empirical spacing = x_i - x_(i-1)
  # theoretical truncated spacing = F(Xt) / (n * f(x_i))
  # shift = |theoretical spacing - empirical spacing|
  # noise_i = |log f(x_i + shift) - log f(x_i)|
  mle_spacing_log_noise_sum <- function(spec, theta) {
    st <- mle_stats(spec, theta)
    if (is.null(st)) return(NA_real_)

    xx <- x[-1]
    emp_spacing <- diff(x)

    pdf1 <- safe_d(spec, xx, st$par, log = FALSE)
    if (length(pdf1) != length(xx) || any(!is.finite(pdf1)) || any(pdf1 <= 0)) return(NA_real_)

    theor_spacing <- st$Ft / (n * pdf1)
    shift <- abs(theor_spacing - emp_spacing)

    pdf2 <- safe_d(spec, xx + shift, st$par, log = FALSE)
    if (length(pdf2) != length(xx) || any(!is.finite(pdf2)) || any(pdf2 <= 0)) return(NA_real_)

    noise <- abs(log(pdf2) - log(pdf1))
    if (any(!is.finite(noise))) return(NA_real_)

    sum(noise)
  }

  mle_fits <- list()
  mle_skipped <- character(0)

  for (dn in dists) {
    spec <- specs[[dn]]
    starts <- make_starts(spec$start, n_starts, seed_offset = 500 + match(dn, dists), wide = FALSE)

    obj <- function(th) {
      st <- mle_stats(spec, th)
      if (is.null(st) || !is.finite(st$sum_ll)) return(1e100)
      -st$sum_ll
    }

    opt <- optimize_fit(spec, obj, starts, maxit = 9000)

    if (is.null(opt)) {
      mle_skipped <- c(mle_skipped, paste0(dn, ": MLE optimization failed"))
      next
    }

    st <- mle_stats(spec, opt$theta)

    if (is.null(st)) {
      mle_skipped <- c(mle_skipped, paste0(dn, ": MLE final evaluation failed"))
      next
    }

    mle_fits[[dn]] <- c(list(name = dn, spec = spec, theta = opt$theta), st)
  }

  if (!length(mle_fits)) {
    message("No MLE fit possible.")
    return(invisible(NULL))
  }

  mle_best <- mle_fits[[which.max(vapply(mle_fits, function(z) z$sum_ll, numeric(1)))]]

  mle_best_dist <- mle_best$name
  mle_params <- mle_best$par
  mle_best_sum_ll <- mle_best$sum_ll
  mle_best_mean_ll <- mle_best$mean_ll
  mle_maxFnInf_best <- mle_best$Ninf
  mle_pct_complete_best <- mle_best$comp_pct
  mle_noise_sum <- mle_spacing_log_noise_sum(mle_best$spec, mle_best$theta)

  # ----------------------------
  # MLE bounds
  # ----------------------------
  mle_E_floor_sum_ll <- stats::setNames(rep(NA_real_, length(E)), names_E)
  mle_E_maxFnInf_bounds <- list()
  mle_E_pct_complete_bounds <- list()
  mle_E_params_bounds <- list()
  mle_E_J_bounds <- list()

  solve_mle_bound <- function(direction = c("lower", "upper"), floor_ll) {
    direction <- match.arg(direction)
    spec <- mle_best$spec
    theta0 <- mle_best$theta

    starts <- make_starts(
      theta0,
      n_starts_bounds,
      seed_offset = 8000 + round(abs(floor_ll)) + ifelse(direction == "upper", 1, 2),
      wide = TRUE
    )

    pen <- 1e14

    obj <- function(th) {
      st <- mle_stats(spec, th)
      if (is.null(st) || !is.finite(st$Ninf) || !is.finite(st$sum_ll)) return(1e100)

      viol <- max(0, floor_ll - st$sum_ll)
      primary <- if (direction == "upper") -log(st$Ninf) else log(st$Ninf)

      primary + pen * viol^2
    }

    opt <- optimize_fit(spec, obj, starts, maxit = 16000)
    if (is.null(opt)) return(NULL)

    st <- mle_stats(spec, opt$theta)
    if (is.null(st) || !is.finite(st$sum_ll) || st$sum_ll < floor_ll - 1e-5) return(NULL)

    c(theta = list(opt$theta), list(par = st$par, Ninf = st$Ninf, sum_ll = st$sum_ll, comp_pct = st$comp_pct))
  }

  if (is.finite(mle_noise_sum)) {
    mle_E_floor_sum_ll <- stats::setNames(mle_best_sum_ll - E * mle_noise_sum, names_E)

    for (ii in seq_along(E)) {
      nm <- names_E[ii]
      floor_ll <- mle_E_floor_sum_ll[ii]

      lo <- solve_mle_bound("lower", floor_ll)
      hi <- solve_mle_bound("upper", floor_ll)

      lower_N <- if (!is.null(lo)) lo$Ninf else NA_real_
      upper_N <- if (!is.null(hi)) hi$Ninf else NA_real_
      if (is.finite(upper_N) && upper_N >= .Machine$double.xmax * 0.99) {
        upper_N <- Inf
      }

      mle_E_maxFnInf_bounds[[nm]] <- c(lower = lower_N, upper = upper_N)
      mle_E_pct_complete_bounds[[nm]] <- c(
        lower = if (!is.na(upper_N)) 100 * n / upper_N else NA_real_,
        upper = if (is.finite(lower_N)) 100 * n / lower_N else NA_real_
      )
      mle_E_params_bounds[[nm]] <- list(
        lower = if (!is.null(lo)) lo$par else NULL,
        upper = if (!is.null(hi)) hi$par else NULL
      )
      mle_E_J_bounds[[nm]] <- c(
        lower = if (!is.null(lo)) lo$sum_ll else NA_real_,
        upper = if (!is.null(hi)) hi$sum_ll else NA_real_
      )
    }
  }

  # ----------------------------
  # Times to event fractions from MLE fit
  # ----------------------------
  mle_times_to_pct_best <- NULL
  if (!is.null(mle_best$spec$qfun) && is.function(mle_best$spec$qfun)) {
    tbest <- safe_q(mle_best$spec, pcts, mle_best$par)
    mle_times_to_pct_best <- stats::setNames(as.numeric(tbest), paste0(round(100 * pcts), "%"))
  }

  # ----------------------------
  # Output
  # ----------------------------
  out <- list(
    x = x,
    n_events = n,
    t_max = Xt,
    dists = dists,
    envelope_expansion = E,
    envelope_labels = names_E,

    # MLE fit
    mle_fit = TRUE,
    mle_best_dist = mle_best_dist,
    mle_params = mle_params,
    mle_best_sum_ll = mle_best_sum_ll,
    mle_best_mean_ll = mle_best_mean_ll,
    mle_maxFnInf_best = mle_maxFnInf_best,
    mle_pct_complete_best = mle_pct_complete_best,
    mle_spacing_log_noise_sum = mle_noise_sum,
    mle_E_floor_sum_ll = mle_E_floor_sum_ll,
    mle_E_maxFnInf_bounds = mle_E_maxFnInf_bounds,
    mle_E_pct_complete_bounds = mle_E_pct_complete_bounds,
    mle_E_params_bounds = mle_E_params_bounds,
    mle_E_J_bounds = mle_E_J_bounds,
    mle_env_maxFnInf_bounds = mle_E_maxFnInf_bounds,
    mle_env_pct_complete_bounds = mle_E_pct_complete_bounds,
    mle_maxFnInf_bounds =
      if ("E100" %in% names(mle_E_maxFnInf_bounds)) mle_E_maxFnInf_bounds[["E100"]]
      else if (length(mle_E_maxFnInf_bounds)) mle_E_maxFnInf_bounds[[1]]
      else c(lower = NA_real_, upper = NA_real_),
    mle_pct_complete_bounds =
      if ("E100" %in% names(mle_E_pct_complete_bounds)) mle_E_pct_complete_bounds[["E100"]]
      else if (length(mle_E_pct_complete_bounds)) mle_E_pct_complete_bounds[[1]]
      else c(lower = NA_real_, upper = NA_real_),
    mle_times_to_pct_best = mle_times_to_pct_best,
    mle_skipped = mle_skipped,

    # MLE-derived CDF diagnostics and compatibility aliases
    best_dist = mle_best_dist,
    params = mle_params,
    best_abs_D = mle_best$cdf_mean_absD,
    best_mean_abs_D = mle_best$cdf_mean_absD,
    best_median_abs_D = mle_best$cdf_median_absD,
    best_sum_abs_D = mle_best$cdf_sum_absD,
    maxFnInf_best = mle_maxFnInf_best,
    pct_complete_best = mle_pct_complete_best,
    times_to_pct_best = mle_times_to_pct_best,
    maxFnInf_bounds =
      if ("E100" %in% names(mle_E_maxFnInf_bounds)) mle_E_maxFnInf_bounds[["E100"]]
      else if (length(mle_E_maxFnInf_bounds)) mle_E_maxFnInf_bounds[[1]]
      else c(lower = NA_real_, upper = NA_real_),
    pct_complete_bounds =
      if ("E100" %in% names(mle_E_pct_complete_bounds)) mle_E_pct_complete_bounds[["E100"]]
      else if (length(mle_E_pct_complete_bounds)) mle_E_pct_complete_bounds[[1]]
      else c(lower = NA_real_, upper = NA_real_),

    skipped = mle_skipped
  )

  class(out) <- "tailfit"
  out
}


#' @export
print.tailfit <- function(x, ...) {
  fmt <- function(z, digits = 3) {
    if (length(z) == 0 || all(is.na(z))) return("NA")
    z <- as.numeric(z)
    if (length(z) == 1 && is.finite(z)) {
      return(format(round(z, digits), nsmall = ifelse(abs(z) < 10, digits, 0), scientific = FALSE, trim = TRUE))
    }
    format(round(z, digits), scientific = FALSE, trim = TRUE)
  }

  fmt_int <- function(z) {
    if (length(z) == 0 || all(is.na(z)) || !is.finite(as.numeric(z)[1])) return("NA")
    format(round(as.numeric(z)[1]), scientific = FALSE, trim = TRUE)
  }

  fmt_pct <- function(z) {
    if (length(z) == 0 || all(is.na(z)) || !is.finite(as.numeric(z)[1])) return("NA")
    format(round(as.numeric(z)[1], 1), nsmall = 1, scientific = FALSE, trim = TRUE)
  }

  fmt_par <- function(par) {
    if (is.null(par) || length(par) == 0 || all(is.na(par))) return("NA")
    paste(paste0(names(par), " = ", fmt(unname(par), 3)), collapse = " ")
  }

  cat("tailfit\n")
  cat("-------\n")
  cat("n events: ", x$n_events, "\n", sep = "")

  cat("\nBest fit function:\n")
  cat("family: ", x$mle_best_dist, "\n", sep = "")
  cat("parameters: ", fmt_par(x$mle_params), "\n", sep = "")

  cat("\nMLE fit (truncated log-likelihood):\n")
  cat("best sum logLik: ", fmt(x$mle_best_sum_ll, 3), "\n", sep = "")
  cat("best mean logLik: ", fmt(x$mle_best_mean_ll, 3), "\n", sep = "")
  cat("spacing-log noise sum: ", fmt(x$mle_spacing_log_noise_sum, 3), "\n", sep = "")

  cat("\nCDF estimators based on function parameters of the best MLE fit:\n")
  cat("Estimated total events (N\u221e): ", fmt_int(x$mle_maxFnInf_best), "\n", sep = "")
  cat("Estimated completeness: ", fmt_pct(x$mle_pct_complete_best), "%\n", sep = "")

  if (length(x$mle_E_maxFnInf_bounds)) {
    cat("\nCDF envelope bounds based on MLE noise:\n")
    for (nm in names(x$mle_E_maxFnInf_bounds)) {
      b <- x$mle_E_maxFnInf_bounds[[nm]]
      cb <- x$mle_E_pct_complete_bounds[[nm]]
      cat("  ", nm, ":\n", sep = "")
      cat("    Estimated total events (N\u221e): ", fmt_int(b["lower"]), " to ", fmt_int(b["upper"]), "\n", sep = "")
      cat("    Estimated completeness: ", fmt_pct(cb["lower"]), "% to ", fmt_pct(cb["upper"]), "%\n", sep = "")
    }
  }

  if (!is.null(x$mle_times_to_pct_best) && length(x$mle_times_to_pct_best)) {
    cat("\nTime to reach % of total expected events:\n")
    for (nm in names(x$mle_times_to_pct_best)) {
      cat("  ", nm, ": ", fmt(x$mle_times_to_pct_best[[nm]], 2), "\n", sep = "")
    }
  }

  invisible(x)
}


#' @export
plot.tailfit <- function(x, envelope = NULL, extend_prob = 0.99, ...) {
  xx <- x$x

  if (is.null(xx) || !length(xx)) stop("No x data in tailfit object.")
  if (!isTRUE(x$mle_fit) || !is.finite(x$mle_maxFnInf_best)) stop("No MLE fit in object.")

  spec <- make_plot_spec(x$mle_best_dist)
  if (is.null(spec)) stop("No plot spec for MLE best distribution.")

  par <- as.list(x$mle_params)

  if (is.null(envelope)) {
    if ("E100" %in% names(x$mle_E_maxFnInf_bounds)) {
      envelope <- "E100"
    } else if (length(x$mle_E_maxFnInf_bounds)) {
      envelope <- names(x$mle_E_maxFnInf_bounds)[1]
    }
  }

  b <- NULL
  if (!is.null(envelope) && envelope %in% names(x$mle_E_maxFnInf_bounds)) {
    b <- x$mle_E_maxFnInf_bounds[[envelope]]
  }

  extend_prob <- suppressWarnings(as.numeric(extend_prob)[1])
  if (!is.finite(extend_prob) || extend_prob <= 0 || extend_prob >= 1) {
    extend_prob <- 0.99
  }

  xmax <- max(xx)
  x_ext <- tryCatch(
    do.call(spec$qfun, c(list(p = extend_prob), par)),
    error = function(e) NA_real_
  )
  if (is.finite(x_ext) && x_ext > xmax) xmax <- x_ext
  xmax <- max(xmax, max(xx) * 1.05)

  grid <- seq(min(xx), xmax, length.out = 600)
  Fg <- do.call(spec$pfun, c(list(q = grid), par))
  pred <- x$mle_maxFnInf_best * Fg

  y_upper <- max(seq_along(xx), pred, na.rm = TRUE)
  if (!is.null(b) && all(is.finite(b))) {
    y_upper <- max(y_upper, b["upper"] * Fg, na.rm = TRUE)
  }
  y_upper <- y_upper * 1.05

  graphics::plot(
    xx,
    seq_along(xx),
    pch = 16,
    xlab = "Time to event x",
    ylab = "Cumulative events",
    main = paste0("CDF plot, based on MLE best fit: ", x$mle_best_dist),
    xlim = c(min(xx), xmax),
    ylim = c(0, y_upper),
    ...
  )

  graphics::lines(grid, pred, lwd = 2)

  if (!is.null(b) && all(is.finite(b))) {
    graphics::lines(grid, b["lower"] * Fg, lty = 2)
    graphics::lines(grid, b["upper"] * Fg, lty = 2)
  }

  graphics::abline(h = length(xx), lty = 3)
  graphics::abline(h = x$mle_maxFnInf_best, lty = 3)

  graphics::legend(
    "bottomright",
    inset = c(0.03, 0.03),
    legend = c("Observed events", "Best fit", "Bounds"),
    pch = c(16, NA, NA),
    lty = c(NA, 1, 2),
    lwd = c(NA, 2, 1),
    bty = "n",
    cex = 0.8
  )

  # Summary labels placed under the legend to avoid right-edge clipping.
  usr <- graphics::par("usr")
  x_info <- usr[1] + 0.08 * diff(usr[1:2])
  y_info <- usr[4] - 0.06 * diff(usr[3:4])
  dy_info <- 0.045 * diff(usr[3:4])

  graphics::text(
    x_info,
    y_info,
    labels = paste0("Observed = ", length(xx)),
    adj = c(0, 1),
    cex = 0.8
  )

  graphics::text(
    x_info,
    y_info - dy_info,
    labels = paste0("Estimated total = ", round(x$mle_maxFnInf_best)),
    adj = c(0, 1),
    cex = 0.8
  )

  invisible(x)
}

make_plot_spec <- function(dist) {
  switch(
    dist,
    lnorm = list(pfun = stats::plnorm, qfun = stats::qlnorm),
    gamma = list(pfun = stats::pgamma, qfun = stats::qgamma),
    weibull = list(pfun = stats::pweibull, qfun = stats::qweibull),
    chisq = list(pfun = stats::pchisq, qfun = stats::qchisq),
    f = list(pfun = stats::pf, qfun = stats::qf),
    NULL
  )
}
