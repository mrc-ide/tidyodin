# Tidy the parameters from a set of monty samples

Convert the parameters in a `monty_samples` object into a long-format
tibble with one row per parameter, draw and chain. Array parameters
(`beta[1]`, `beta[2]`) are split into a parameter name and index
columns, and parameters from a grouped packer (`beta<a>`) into a
parameter name and a `group` column.

## Usage

``` r
tidy_pars(samples, vars = NULL, dims = NULL, labels = NULL)
```

## Arguments

- samples:

  A `monty_samples` object, from
  [`monty::monty_sample()`](https://mrc-ide.github.io/monty/reference/monty_sample.html)
  or
  [`run_parameter_sets()`](https://mrc-ide.github.io/tidyodin/reference/run_parameter_sets.md).

- vars:

  Optional character vector of parameters to keep.

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

A tibble with columns `.chain`, `.iteration`, `.draw`, `parameter`,
index columns for array parameters, `group` for grouped parameters, and
`value`.

## Details

For a wide format with one column per parameter, use
[`posterior::as_draws_df()`](https://mc-stan.org/posterior/reference/draws_df.html),
which monty supports directly. The `.chain`, `.iteration` and `.draw`
columns here match those.
