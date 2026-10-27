# tailcdf

`tailcdf` fits parametric distributions to right-truncated event times
to estimate the unobserved tail of a cumulative event process.

It estimates:

- Eventual total number of events
- Percentage of events already observed
- Empirical likelihood-envelope bounds
- Times to reach 60%, 70%, 80%, and 90% completeness

The bounds are not conventional confidence intervals.

## Installation

After CRAN release:

```r
install.packages("tailcdf")
```

## Example

```r
library(tailcdf)

set.seed(123)
x <- rweibull(100, shape = 2, scale = 10)
x <- x[x <= 12]

fit <- tailfit(x, dists = "weibull")
print(fit)
plot(fit)
```

## Input template

A customizable input template is included in the package.
It provides distribution selection, bound-width adjustment,
and completeness reporting options.
