#' Tidy the output of a dust2 simulation
#'
#' Convert the array returned by [dust2::dust_system_simulate()] (or the
#' state returned by [dust2::dust_system_state()]) into a long-format
#' tibble with one row per variable, time, particle and group.
#'
#' The array on its own does not say whether its middle dimensions are
#' particles or parameter groups, nor what its times were, so the system
#' that produced it and the times are needed too.
#'
#' @param y An array returned by [dust2::dust_system_simulate()], or by
#'   [dust2::dust_system_state()] if `times` is `NULL`.
#'
#' @param sys The `dust_system` that produced `y`.
#'
#' @param times The times passed to [dust2::dust_system_simulate()].
#'   Use `NULL` if `y` is a state with no time dimension.
#'
#' @param index_state The `index_state` used to produce `y`, if any.
#'
#' @param vars Optional character vector of variables to keep, for
#'   example `c("I", "incidence")`. Selecting variables before reshaping
#'   is much cheaper than filtering afterwards.
#'
#' @param dims Optional names for the dimensions of array variables.
#'   Either an unnamed character vector applied by position to every
#'   array variable (`c("age", "vax")`), or a named list giving the
#'   dimension names of each variable
#'   (`list(S = c("age", "vax"), R = "age")`). Dimensions that are not
#'   named produce columns `index1`, `index2`, and so on.
#'
#' @param labels Optional named list of labels for dimensions, for
#'   example `list(age = c("child", "adult"))`. Labelled dimensions
#'   become factors. Use the name `group` to label parameter groups.
#'
#' @return A tibble with columns `group` and `particle` (where the
#'   system has them), `time`, `variable`, one column per array
#'   dimension (`NA` for variables without that dimension) and `value`.
#'
#' @export
#' @examplesIf requireNamespace("odin2", quietly = TRUE) && interactive()
#' sir <- odin2::odin({
#'   deriv(S) <- -beta * S * I / N
#'   deriv(I) <- beta * S * I / N - gamma * I
#'   initial(S) <- N - 1
#'   initial(I) <- 1
#'   N <- parameter(1000)
#'   beta <- parameter(0.4)
#'   gamma <- parameter(0.1)
#' })
#' sys <- dust2::dust_system_create(sir, list())
#' dust2::dust_system_set_state_initial(sys)
#' t <- seq(0, 100)
#' y <- dust2::dust_system_simulate(sys, t)
#' tidy_simulation(y, sys, t)
tidy_simulation <- function(y, sys, times, index_state = NULL, vars = NULL,
                            dims = NULL, labels = NULL) {
  if (!inherits(sys, "dust_system")) {
    cli::cli_abort("{.arg sys} must be a {.cls dust_system}")
  }
  names <- sys$packer_state$names()
  if (!is.null(index_state)) {
    names <- names[index_state]
  }
  trailing <- list()
  if (sys$preserve_particle_dimension) {
    trailing$particle <- seq_len(sys$n_particles)
  }
  if (sys$preserve_group_dimension) {
    trailing$group <- seq_len(sys$n_groups)
  }
  if (!is.null(times)) {
    trailing$time <- times
  }
  melt_state(y, names, trailing, vars = vars, dims = dims, labels = labels)
}


#' Tidy the trajectories from a particle filter
#'
#' Convert the trajectories held by a filter or unfilter after
#' [dust2::dust_likelihood_run()] into a long-format tibble. Times and
#' group names are read from the filter.
#'
#' @param obj A `dust_likelihood` object, created by
#'   [dust2::dust_filter_create()] or [dust2::dust_unfilter_create()],
#'   that has been run with `save_trajectories = TRUE`.
#'
#' @param select_random_particle Logical, passed to
#'   [dust2::dust_likelihood_last_trajectories()]: return a single
#'   randomly selected particle rather than all of them.
#'
#' @param index_state The `index_state` passed to
#'   [dust2::dust_likelihood_run()], if any.
#'
#' @inheritParams tidy_simulation
#'
#' @return A tibble with columns `group` and `particle` (where present),
#'   `time`, `variable`, one column per array dimension and `value`.
#'   Times are the times of the data.
#'
#' @export
tidy_filter <- function(obj, select_random_particle = FALSE,
                        index_state = NULL, vars = NULL, dims = NULL,
                        labels = NULL) {
  if (!inherits(obj, "dust_likelihood")) {
    cli::cli_abort("{.arg obj} must be a {.cls dust_likelihood}")
  }
  y <- dust2::dust_likelihood_last_trajectories(
    obj, select_random_particle = select_random_particle)
  names <- likelihood_state_names(obj)
  if (!is.null(index_state)) {
    names <- names[index_state]
  }
  trailing <- likelihood_trailing(
    obj, particle = !(select_random_particle && inherits(obj, "dust_filter")))
  trailing$time <- obj$inputs$time
  melt_state(y, names, trailing, vars = vars, dims = dims, labels = labels)
}


likelihood_state_names <- function(obj, call = rlang::caller_env()) {
  packer <- obj$packer_state
  if (is.null(packer)) {
    cli::cli_abort(
      c("The likelihood has not been run yet, so its state is not known",
        i = "Run it with {.fn dust2::dust_likelihood_run} first"),
      call = call)
  }
  packer$names()
}


## The particle and group dimensions of trajectories from a likelihood
likelihood_trailing <- function(obj, particle = TRUE) {
  trailing <- list()
  if (particle && isTRUE(obj$preserve_particle_dimension)) {
    trailing$particle <- seq_len(obj$n_particles)
  }
  if (isTRUE(obj$preserve_group_dimension)) {
    trailing$group <- obj$groups %||% seq_len(obj$n_groups)
  }
  trailing
}
