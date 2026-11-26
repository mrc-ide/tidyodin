# Tidy the output of a dust2 simulation

Convert the array returned by
[`dust2::dust_system_simulate()`](https://mrc-ide.github.io/dust2/reference/dust_system_simulate.html)
(or the state returned by
[`dust2::dust_system_state()`](https://mrc-ide.github.io/dust2/reference/dust_system_state.html))
into a long-format tibble with one row per variable, time, particle and
group.

## Usage

``` r
tidy_simulation(
  y,
  sys,
  times,
  index_state = NULL,
  vars = NULL,
  dims = NULL,
  labels = NULL
)
```

## Arguments

- y:

  An array returned by
  [`dust2::dust_system_simulate()`](https://mrc-ide.github.io/dust2/reference/dust_system_simulate.html),
  or by
  [`dust2::dust_system_state()`](https://mrc-ide.github.io/dust2/reference/dust_system_state.html)
  if `times` is `NULL`.

- sys:

  The `dust_system` that produced `y`.

- times:

  The times passed to
  [`dust2::dust_system_simulate()`](https://mrc-ide.github.io/dust2/reference/dust_system_simulate.html).
  Use `NULL` if `y` is a state with no time dimension.

- index_state:

  The `index_state` used to produce `y`, if any.

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

A tibble with columns `group` and `particle` (where the system has
them), `time`, `variable` (a factor, with levels in the order the model
stores its state), one column per array dimension (`NA` for variables
without that dimension) and `value`.

## Details

The array on its own does not say whether its middle dimensions are
particles or parameter groups, nor what its times were, so the system
that produced it and the times are needed too.

## Examples

``` r
if (FALSE) { # requireNamespace("odin2", quietly = TRUE) && interactive()
sir <- odin2::odin({
  deriv(S) <- -beta * S * I / N
  deriv(I) <- beta * S * I / N - gamma * I
  initial(S) <- N - 1
  initial(I) <- 1
  N <- parameter(1000)
  beta <- parameter(0.4)
  gamma <- parameter(0.1)
})
sys <- dust2::dust_system_create(sir, list())
dust2::dust_system_set_state_initial(sys)
t <- seq(0, 100)
y <- dust2::dust_system_simulate(sys, t)
tidy_simulation(y, sys, t)
}
```
