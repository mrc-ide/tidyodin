# Tidy trajectories saved alongside monty samples

Convert the trajectories, final state or snapshots saved with a set of
samples into a long-format tibble. These are saved when the likelihood
is created with `save_trajectories`, `save_state` or `save_snapshots`
(see
[`dust2::dust_likelihood_monty()`](https://mrc-ide.github.io/dust2/reference/dust_likelihood_monty.html)),
and are also present in the output of
[`run_parameter_sets()`](https://mrc-ide.github.io/tidyodin/reference/run_parameter_sets.md).

## Usage

``` r
tidy_trajectories(
  samples,
  obj = NULL,
  times = NULL,
  save_trajectories = NULL,
  n_draws = NULL,
  vars = NULL,
  dims = NULL,
  labels = NULL
)

tidy_state(
  samples,
  obj = NULL,
  n_draws = NULL,
  vars = NULL,
  dims = NULL,
  labels = NULL
)

tidy_snapshots(
  samples,
  obj = NULL,
  times,
  n_draws = NULL,
  vars = NULL,
  dims = NULL,
  labels = NULL
)
```

## Arguments

- samples:

  A `monty_samples` object, from
  [`monty::monty_sample()`](https://mrc-ide.github.io/monty/reference/monty_sample.html)
  or
  [`run_parameter_sets()`](https://mrc-ide.github.io/tidyodin/reference/run_parameter_sets.md).

- obj:

  The filter or unfilter (a `dust_likelihood` object) that the samples
  were created with. Not needed for output from
  [`run_parameter_sets()`](https://mrc-ide.github.io/tidyodin/reference/run_parameter_sets.md).

- times:

  Times for the trajectories. By default these are taken from `obj` (the
  times of the data) or from the samples. For `tidy_snapshots()` these
  are the snapshot times and are required.

- save_trajectories:

  If the likelihood was created with `save_trajectories` as a character
  vector of variable names, that same vector.

- n_draws:

  Optional number of draws to keep, selected at random from all chains.

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

A tibble with columns `chain`, `iteration`, `sample`, `group` and
`particle` (where present), `time` (except for `tidy_state()`),
`variable`, one column per array dimension and `value`.

## Details

Samples from
[`monty::monty_sample()`](https://mrc-ide.github.io/monty/reference/monty_sample.html)
do not record which state variables were saved or at what times, so the
filter used to create the likelihood is needed as `obj`. Output from
[`run_parameter_sets()`](https://mrc-ide.github.io/tidyodin/reference/run_parameter_sets.md)
carries this information with it and `obj` can be omitted.

Trajectories can be very large (variables x times x draws). Use `vars`
and `n_draws` to reduce them before reshaping.
