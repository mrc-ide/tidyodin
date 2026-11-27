#' Tidy the parameters from a set of monty samples
#'
#' Convert the parameters in a `monty_samples` object into a long-format
#' tibble with one row per parameter, sample and chain. Array parameters
#' (`beta[1]`, `beta[2]`) are split into a parameter name and index
#' columns, and parameters from a grouped packer (`beta<a>`) into a
#' parameter name and a `group` column.
#'
#' For a wide format with one column per parameter, use
#' [posterior::as_draws_df()], which monty supports directly. Its
#' `.chain`, `.iteration` and `.draw` columns hold the same values as
#' `chain`, `iteration` and `sample` here.
#'
#' @param samples A `monty_samples` object, from [monty::monty_sample()]
#'   or [run_parameter_sets()].
#'
#' @param vars Optional character vector of parameters to keep.
#'
#' @inheritParams tidy_simulation
#'
#' @return A tibble with columns `chain`, `iteration`, `sample`,
#'   `parameter` (a factor, with levels in the order of the samples),
#'   index columns for array parameters, `group` for
#'   grouped parameters, and `value`.
#'
#' @export
tidy_pars <- function(samples, vars = NULL, dims = NULL, labels = NULL) {
  samples <- check_samples(samples)
  x <- samples$pars
  n <- dim(x)
  dim(x) <- c(n[[1]], n[[2]] * n[[3]])
  melt_state(x, rownames(samples$pars),
             list(sample = seq_len(n[[2]] * n[[3]])),
             vars = vars, dims = dims, labels = labels,
             name_col = "parameter", n_samples = n[[2]])
}


#' Tidy the log density from a set of monty samples
#'
#' @inheritParams tidy_pars
#'
#' @return A tibble with columns `chain`, `iteration`, `sample` and
#'   `density`, which can be joined to the output of the other `tidy_`
#'   functions by `sample`.
#'
#' @export
tidy_density <- function(samples) {
  samples <- check_samples(samples)
  n <- dim(samples$density)
  tibble::tibble(
    chain = rep(seq_len(n[[2]]), each = n[[1]]),
    iteration = rep(seq_len(n[[1]]), times = n[[2]]),
    sample = seq_len(n[[1]] * n[[2]]),
    density = as.vector(samples$density))
}


#' Tidy trajectories saved alongside monty samples
#'
#' Convert the trajectories, final state or snapshots saved with a set
#' of samples into a long-format tibble. These are saved when the
#' likelihood is created with `save_trajectories`, `save_state` or
#' `save_snapshots` (see [dust2::dust_likelihood_monty()]), and are also
#' present in the output of [run_parameter_sets()].
#'
#' Samples from [monty::monty_sample()] do not record which state
#' variables were saved or at what times, so the filter used to create
#' the likelihood is needed as `obj`. Output from [run_parameter_sets()]
#' carries this information with it and `obj` can be omitted.
#'
#' A filter only knows the names of its state variables once it has
#' been run. If the samples were saved and the filter has been created
#' again in a new session, run it once with
#' [dust2::dust_likelihood_run()], using any valid parameters, before
#' passing it here.
#'
#' Trajectories can be very large (variables x times x draws). Use
#' `vars` and `n_samples` to reduce them before reshaping.
#'
#' @param samples A `monty_samples` object, from [monty::monty_sample()]
#'   or [run_parameter_sets()].
#'
#' @param obj The filter or unfilter (a `dust_likelihood` object) that
#'   the samples were created with. Not needed for output from
#'   [run_parameter_sets()].
#'
#' @param times Times for the trajectories. By default these are taken
#'   from `obj` (the times of the data) or from the samples. For
#'   `tidy_snapshots()` these are the snapshot times and are required.
#'
#' @param save_trajectories If the likelihood was created with
#'   `save_trajectories` as a character vector of variable names, that
#'   same vector.
#'
#' @param n_samples Optional number of samples to keep, selected at random
#'   from all chains.
#'
#' @inheritParams tidy_simulation
#'
#' @return A tibble with columns `chain`, `iteration`, `sample`,
#'   `group` and `particle` (where present), `time` (except for
#'   `tidy_state()`), `variable`, one column per array dimension and
#'   `value`.
#'
#' @export
tidy_trajectories <- function(samples, obj = NULL, times = NULL,
                              save_trajectories = NULL, n_samples = NULL,
                              vars = NULL, dims = NULL, labels = NULL) {
  tidy_observation(samples, "trajectories", obj, times, save_trajectories,
                   n_samples, vars, dims, labels)
}


#' @rdname tidy_trajectories
#' @export
tidy_state <- function(samples, obj = NULL, n_samples = NULL, vars = NULL,
                       dims = NULL, labels = NULL) {
  tidy_observation(samples, "state", obj, NULL, NULL, n_samples, vars, dims,
                   labels)
}


#' @rdname tidy_trajectories
#' @export
tidy_snapshots <- function(samples, obj = NULL, times, n_samples = NULL,
                           vars = NULL, dims = NULL, labels = NULL) {
  tidy_observation(samples, "snapshots", obj, times, NULL, n_samples, vars,
                   dims, labels)
}


tidy_observation <- function(samples, what, obj, times, save_trajectories,
                             n_samples, vars, dims, labels,
                             call = rlang::caller_env()) {
  meta <- attr(samples, "tidyodin")
  samples <- check_samples(samples, call)
  x <- samples$observations[[what]]
  if (is.null(x)) {
    cli::cli_abort(
      c("These samples do not contain any {what}",
        i = paste("Create the likelihood with",
                  "{.code {observation_arg[[what]]}} in",
                  "{.fn dust2::dust_likelihood_monty}")),
      call = call)
  }
  n <- dim(samples$density)
  d <- dim(x)
  if (!is.numeric(x) || length(d) < 3 || !identical(utils::tail(d, 2), n)) {
    cli::cli_abort(
      paste("Can't tidy {.code observations${what}}: expected an array",
            "whose last two dimensions are samples and chains"),
      call = call)
  }

  if (!is.null(obj)) {
    if (inherits(obj, "dust_likelihood")) {
      names <- likelihood_state_names(obj, call)
      trailing <- likelihood_trailing(obj,
                                      particle = !inherits(obj, "dust_filter"))
      times <- times %||% obj$inputs$time
      packer <- obj$packer_state
    } else if (inherits(obj, "dust_system")) {
      names <- obj$packer_state$names()
      trailing <- list()
      packer <- obj$packer_state
    } else {
      cli::cli_abort(
        "{.arg obj} must be a {.cls dust_likelihood} or a {.cls dust_system}",
        call = call)
    }
    if (what == "trajectories" && is.character(save_trajectories)) {
      names <- names[packer$subset(save_trajectories)$index]
    }
  } else if (!is.null(meta)) {
    names <- if (what == "trajectories") meta$names_trajectories else
      meta$names_state
    trailing <- list()
    if (meta$n_particles > 1) {
      trailing$particle <- seq_len(meta$n_particles)
    }
    times <- times %||% meta$times
  } else {
    cli::cli_abort(
      c("{.arg obj} is needed to tidy these samples",
        i = "Pass the filter that was used to create the likelihood"),
      call = call)
  }

  if (what != "state") {
    if (is.null(times)) {
      cli::cli_abort("{.arg times} must be provided", call = call)
    }
    trailing$time <- times
  }

  n_total <- n[[1]] * n[[2]]
  lead <- utils::head(d, -2)
  draws <- seq_len(n_total)
  if (!is.null(n_samples) && n_samples < n_total) {
    draws <- sort(sample.int(n_total, n_samples))
    x <- matrix(x, prod(lead), n_total)[, draws, drop = FALSE]
  }
  dim(x) <- c(lead, length(draws))
  trailing$sample <- draws

  melt_state(x, names, trailing, vars = vars, dims = dims, labels = labels,
             n_samples = n[[1]], call = call)
}


observation_arg <- list(trajectories = "save_trajectories = TRUE",
                        state = "save_state = TRUE",
                        snapshots = "save_snapshots")


check_samples <- function(samples, call = rlang::caller_env()) {
  if (!inherits(samples, "monty_samples")) {
    cli::cli_abort("{.arg samples} must be a {.cls monty_samples} object",
                   call = call)
  }
  if (!is.null(attr(samples, "chain"))) {
    samples <- monty::monty_unflatten_chains(samples)
  }
  samples
}


#' Summarise tidy trajectories across draws
#'
#' Collapse a tidy tibble of trajectories to a mean, median and central
#' interval at each time, ready for `ggplot2::geom_ribbon()`. The summary
#' is taken over draws, chains and particles; every other column
#' (variable, time, group, array dimensions and any columns you have
#' added, such as a scenario) is kept as a grouping.
#'
#' @param data A tibble from [tidy_trajectories()], [tidy_simulation()]
#'   or [tidy_filter()].
#'
#' @param interval Width of the central interval, by default 0.95.
#'
#' @return A tibble with the grouping columns and `mean`, `median`,
#'   `lower` and `upper`.
#'
#' @export
summarise_trajectories <- function(data, interval = 0.95) {
  if (!is.data.frame(data) || !("value" %in% names(data))) {
    cli::cli_abort(
      "{.arg data} must be a data frame with a {.val value} column")
  }
  over <- c("value", "chain", "iteration", "sample", "particle")
  by <- setdiff(names(data), over)
  p <- (1 - interval) / 2
  dplyr::summarise(
    data,
    mean = mean(.data$value),
    median = stats::median(.data$value),
    lower = stats::quantile(.data$value, p, names = FALSE),
    upper = stats::quantile(.data$value, 1 - p, names = FALSE),
    .by = dplyr::all_of(by))
}
