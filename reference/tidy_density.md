# Tidy the log density from a set of monty samples

Tidy the log density from a set of monty samples

## Usage

``` r
tidy_density(samples)
```

## Arguments

- samples:

  A `monty_samples` object, from
  [`monty::monty_sample()`](https://mrc-ide.github.io/monty/reference/monty_sample.html)
  or
  [`run_parameter_sets()`](https://mrc-ide.github.io/tidyodin/reference/run_parameter_sets.md).

## Value

A tibble with columns `.chain`, `.iteration`, `.draw` and `density`,
which can be joined to the output of the other `tidy_` functions by
`.draw`.
