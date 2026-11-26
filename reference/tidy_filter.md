# Tidy the trajectories from a particle filter

Convert the trajectories held by a filter or unfilter after
[`dust2::dust_likelihood_run()`](https://mrc-ide.github.io/dust2/reference/dust_likelihood_run.html)
into a long-format tibble. Times and group names are read from the
filter.

## Usage

``` r
tidy_filter(
  obj,
  select_random_particle = FALSE,
  index_state = NULL,
  vars = NULL,
  dims = NULL,
  labels = NULL
)
```

## Arguments

- obj:

  A `dust_likelihood` object, created by
  [`dust2::dust_filter_create()`](https://mrc-ide.github.io/dust2/reference/dust_filter_create.html)
  or
  [`dust2::dust_unfilter_create()`](https://mrc-ide.github.io/dust2/reference/dust_unfilter_create.html),
  that has been run with `save_trajectories = TRUE`.

- select_random_particle:

  Logical, passed to
  [`dust2::dust_likelihood_last_trajectories()`](https://mrc-ide.github.io/dust2/reference/dust_likelihood_last_trajectories.html):
  return a single randomly selected particle rather than all of them.

- index_state:

  The `index_state` passed to
  [`dust2::dust_likelihood_run()`](https://mrc-ide.github.io/dust2/reference/dust_likelihood_run.html),
  if any.

- vars:

  Optional character vector of variables to keep, for example
  `c("I", "incidence")`. Selecting variables before reshaping is much
  cheaper than filtering afterwards.

- dims:

  Optional names for the dimensions of array variables. Either an
  unnamed character vector applied by position to every array variable
  (`c("age", "vax")`), or a named list giving the dimension names of
  each variable (`list(S = c("age", "vax"), R = "age")`). Dimensions
  that are not named produce columns `index1`, `index2`, and so on.

- labels:

  Optional named list of labels for dimensions, for example
  `list(age = c("child", "adult"))`. Labelled dimensions become factors.
  Use the name `group` to label parameter groups.

## Value

A tibble with columns `group` and `particle` (where present), `time`,
`variable`, one column per array dimension and `value`. Times are the
times of the data.
