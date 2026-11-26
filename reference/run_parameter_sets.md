# Run a model over many parameter sets

Simulate an odin2 model from its initial conditions, once for each of a
collection of parameter sets such as posterior draws or scenarios. The
result has the same structure as the output of
[`monty::monty_sample()`](https://mrc-ide.github.io/monty/reference/monty_sample.html),
so it can be passed to
[`tidy_trajectories()`](https://mrc-ide.github.io/tidyodin/reference/tidy_trajectories.md),
[`tidy_pars()`](https://mrc-ide.github.io/tidyodin/reference/tidy_pars.md)
and to monty's own tools such as
[`monty::monty_samples_thin()`](https://mrc-ide.github.io/monty/reference/monty_samples_thin.html).
Use
[`continue_parameter_sets()`](https://mrc-ide.github.io/tidyodin/reference/continue_parameter_sets.md)
to carry on from the end of a fit or of a previous run.

## Usage

``` r
run_parameter_sets(
  generator,
  pars,
  times,
  packer = NULL,
  fixed = NULL,
  n_particles = 1,
  deterministic = NULL,
  n_draws = NULL,
  time_start = NULL,
  save_trajectories = TRUE,
  dt = NULL,
  ode_control = NULL,
  seed = NULL,
  n_threads = 1
)
```

## Arguments

- generator:

  An odin2 model generator, as returned by
  [`odin2::odin()`](https://mrc-ide.github.io/odin2/reference/odin.html).
  When rerunning a fitted model, pass the filter or unfilter used for
  the fit instead: the model, its time step and start time, and whether
  it is deterministic are then taken from the fit rather than given
  again.

- pars:

  The parameter sets. One of:

  - a `monty_samples` object (requires `packer`); chains are kept

  - a data frame with one row per parameter set and one numeric column
    per parameter. Columns `.chain`, `.iteration` and `.draw`, as
    produced by
    [`posterior::as_draws_df()`](https://mc-stan.org/posterior/reference/draws_df.html),
    are not treated as parameters, and `.chain` is used to keep chains

  - a list of parameter lists; without a `packer`, elements that are the
    same in every set are treated as fixed

  - a matrix with one column per parameter set (requires `packer`)

- times:

  Times to save output at, increasing and no earlier than `time_start`.

- packer:

  A
  [`monty::monty_packer()`](https://mrc-ide.github.io/monty/reference/monty_packer.html)
  used to convert between parameter vectors and the list of parameters
  the model needs. Required if `pars` is a `monty_samples` object or a
  matrix.

- fixed:

  Optional named list of parameters held constant across all parameter
  sets, for example `list(N = 1000)`. These are added to every set, and
  replace any fixed values in `packer`. They may not name a parameter
  that varies in `pars`.

- n_particles:

  Number of particles (stochastic realisations) to run for each
  parameter set.

- deterministic:

  Logical, indicating that a stochastic model should be run
  deterministically, replacing random draws with their expectations. Has
  no effect on ODE models, which are always deterministic. By default
  `FALSE`, or as the fit if `generator` is a filter or unfilter.

- n_draws:

  Optional number of parameter sets to select at random from `pars`. If
  given, the result has a single chain.

- time_start:

  The time at which the model starts. By default 0, or as the fit if
  `generator` is a filter or unfilter.

- save_trajectories:

  Either `TRUE`, to save all state variables, or a character vector of
  variable names to save.

- dt, ode_control:

  Time step for discrete-time models and control for ODE models, passed
  to
  [`dust2::dust_system_create()`](https://mrc-ide.github.io/dust2/reference/dust_system_create.html).
  By default those of the fit if `generator` is a filter or unfilter,
  and otherwise the model's defaults.

- seed, n_threads:

  Passed to
  [`dust2::dust_system_create()`](https://mrc-ide.github.io/dust2/reference/dust_system_create.html).

## Value

An object of class `monty_samples`, with elements:

- `pars`: parameters x samples x chains

- `density`: samples x chains, all `NA` as no density is computed

- `observations$trajectories`: state x time x samples x chains, or state
  x particle x time x samples x chains if `n_particles > 1`

- `observations$state`: the final state, laid out as trajectories
  without the time dimension

Everything needed to tidy the result or continue the run is stored as an
attribute, so
[`tidy_trajectories()`](https://mrc-ide.github.io/tidyodin/reference/tidy_trajectories.md)
and
[`continue_parameter_sets()`](https://mrc-ide.github.io/tidyodin/reference/continue_parameter_sets.md)
need no further information.

## Details

dust2 provides
[`dust2::dust_system_simulate()`](https://mrc-ide.github.io/dust2/reference/dust_system_simulate.html),
which runs a system that has already been created. This function does
the work around that: turning the parameter sets into the form the model
needs, creating a system with one group per parameter set, running it,
and arranging the output so that each trajectory stays with the
parameters that produced it.

All parameter sets are run together in a single dust2 system, so every
set must produce a model with the same state (for example, the same
array sizes).

To compare scenarios, call this once per scenario and combine the tidied
results with `dplyr::bind_rows(..., .id = "scenario")`. For a
counterfactual based on a fit, convert the samples to a data frame with
[`posterior::as_draws_df()`](https://mc-stan.org/posterior/reference/draws_df.html),
change the parameters there, and pass the data frame as `pars`.
