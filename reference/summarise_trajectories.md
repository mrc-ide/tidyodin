# Summarise tidy trajectories across draws

Collapse a tidy tibble of trajectories to a mean, median and central
interval at each time, ready for
[`ggplot2::geom_ribbon()`](https://ggplot2.tidyverse.org/reference/geom_ribbon.html).
The summary is taken over draws, chains and particles; every other
column (variable, time, group, array dimensions and any columns you have
added, such as a scenario) is kept as a grouping.

## Usage

``` r
summarise_trajectories(data, interval = 0.95)
```

## Arguments

- data:

  A tibble from
  [`tidy_trajectories()`](https://mrc-ide.github.io/tidyodin/reference/tidy_trajectories.md),
  [`tidy_simulation()`](https://mrc-ide.github.io/tidyodin/reference/tidy_simulation.md)
  or
  [`tidy_filter()`](https://mrc-ide.github.io/tidyodin/reference/tidy_filter.md).

- interval:

  Width of the central interval, by default 0.95.

## Value

A tibble with the grouping columns and `mean`, `median`, `lower` and
`upper`.
