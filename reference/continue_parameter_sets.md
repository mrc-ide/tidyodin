# Continue a simulation from the end of a fit or a previous run

Run a model forward from the final state and time of a set of samples,
with each sample continuing from its own state using its own parameters.
The samples can be a fit from
[`monty::monty_sample()`](https://mrc-ide.github.io/monty/reference/monty_sample.html),
for projecting forward from the end of the data, or the result of
[`run_parameter_sets()`](https://mrc-ide.github.io/tidyodin/reference/run_parameter_sets.md)
or of this function.

## Usage

``` r
continue_parameter_sets(
  samples,
  times,
  obj = NULL,
  packer = NULL,
  fixed = NULL,
  modify = NULL,
  n_particles = NULL,
  deterministic = NULL,
  n_draws = NULL,
  state = NULL,
  time = NULL,
  save_trajectories = NULL,
  dt = NULL,
  seed = NULL,
  n_threads = 1
)
```

## Arguments

- samples:

  A `monty_samples` object, from
  [`monty::monty_sample()`](https://mrc-ide.github.io/monty/reference/monty_sample.html),
  [`run_parameter_sets()`](https://mrc-ide.github.io/tidyodin/reference/run_parameter_sets.md)
  or `continue_parameter_sets()`.

- times:

  Times to save output at, increasing and no earlier than the time being
  continued from.

- obj:

  For a fit, the filter or unfilter the likelihood was created from. Not
  needed for output from
  [`run_parameter_sets()`](https://mrc-ide.github.io/tidyodin/reference/run_parameter_sets.md).

- packer:

  For a fit, the packer the likelihood was created with. Not needed for
  output from
  [`run_parameter_sets()`](https://mrc-ide.github.io/tidyodin/reference/run_parameter_sets.md).

- fixed:

  Optional named list of parameters to hold constant from this point on,
  replacing fixed values used so far.

- modify:

  Optional function applied to each parameter list from this point on,
  in addition to any applied so far.

- n_particles:

  Number of particles to run for each parameter set, each starting from
  that set's final state. If the samples already hold several particles
  per set, each particle continues from its own state and this must not
  be changed.

- deterministic:

  Logical, indicating that a stochastic model should be run
  deterministically. By default this follows the fit (`TRUE` for an
  unfilter, `FALSE` for a filter) or the previous run.

- n_draws:

  Optional number of parameter sets to select at random from `pars`. If
  given, the result has a single chain.

- state, time:

  Optional state and time to continue from, instead of the final state.
  `state` has one column per sample, so a saved snapshot can be passed
  directly: `samples$observations$snapshots[, 1, , ]`.

- save_trajectories:

  Either `TRUE`, to save all state variables, or a character vector of
  variable names. By default the same variables as the previous run, or
  all of them for a fit.

- dt:

  Time step, by default that of the fit or previous run.

- seed, n_threads:

  Passed to
  [`dust2::dust_system_create()`](https://mrc-ide.github.io/dust2/reference/dust_system_create.html).

## Value

An object of class `monty_samples`, as for
[`run_parameter_sets()`](https://mrc-ide.github.io/tidyodin/reference/run_parameter_sets.md).

## Details

For a fit, the final state must have been saved by creating the
likelihood with `save_state = TRUE` (see
[`dust2::dust_likelihood_monty()`](https://mrc-ide.github.io/dust2/reference/dust_likelihood_monty.html));
it cannot be recovered afterwards. The filter is needed too, as it holds
the model, its time step and the time of the last data point.

To start from somewhere other than the end, pass `state` and `time`
explicitly, for example a snapshot saved with `save_snapshots`.
