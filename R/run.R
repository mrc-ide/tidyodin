#' Run a model over many parameter sets
#'
#' Simulate an odin2 model from its initial conditions, once for each of
#' a collection of parameter sets such as posterior draws or scenarios.
#' The result has the same structure as the output of
#' [monty::monty_sample()], so it can be passed to [tidy_trajectories()],
#' [tidy_pars()] and to monty's own tools such as
#' [monty::monty_samples_thin()]. Use [continue_parameter_sets()] to
#' carry on from the end of a fit or of a previous run.
#'
#' dust2 provides [dust2::dust_system_simulate()], which runs a system
#' that has already been created. This function does the work around
#' that: turning the parameter sets into the form the model needs,
#' creating a system with one group per parameter set, running it, and
#' arranging the output so that each trajectory stays with the
#' parameters that produced it.
#'
#' All parameter sets are run together in a single dust2 system, so
#' every set must produce a model with the same state (for example, the
#' same array sizes).
#'
#' To compare scenarios, call this once per scenario and combine the
#' tidied results with `dplyr::bind_rows(..., .id = "scenario")`.
#'
#' @param generator An odin2 model generator, as returned by
#'   [odin2::odin()].
#'
#' @param pars The parameter sets. One of:
#'
#'   * a `monty_samples` object (requires `packer`); chains are kept
#'   * a data frame with one row per parameter set and one numeric
#'     column per parameter
#'   * a list of parameter lists; without a `packer`, elements that are
#'     the same in every set are treated as fixed
#'   * a matrix with one column per parameter set (requires `packer`)
#'
#' @param times Times to save output at, increasing and no earlier than
#'   `time_start`.
#'
#' @param packer A [monty::monty_packer()] used to convert between
#'   parameter vectors and the list of parameters the model needs.
#'   Required if `pars` is a `monty_samples` object or a matrix.
#'
#' @param fixed Optional named list of parameters held constant across
#'   all parameter sets, for example `list(N = 1000)`. These are added
#'   to every set, and replace any fixed values in `packer`. They may
#'   not name a parameter that varies in `pars`.
#'
#' @param modify Optional function applied to each parameter list before
#'   the model is run, for example to scale a parameter for a
#'   counterfactual: `function(p) { p$beta <- p$beta * 0.8; p }`.
#'
#' @param n_particles Number of particles (stochastic realisations) to
#'   run for each parameter set.
#'
#' @param deterministic Logical, indicating that a stochastic model
#'   should be run deterministically, replacing random draws with their
#'   expectations. Has no effect on ODE models, which are always
#'   deterministic.
#'
#' @param n_draws Optional number of parameter sets to select at random
#'   from `pars`. If given, the result has a single chain.
#'
#' @param time_start The time at which the model starts.
#'
#' @param save_trajectories Either `TRUE`, to save all state variables,
#'   or a character vector of variable names to save.
#'
#' @param dt,ode_control,seed,n_threads Passed to
#'   [dust2::dust_system_create()].
#'
#' @return An object of class `monty_samples`, with elements:
#'
#'   * `pars`: parameters x samples x chains
#'   * `density`: samples x chains, all `NA` as no density is computed
#'   * `observations$trajectories`: state x time x samples x chains, or
#'     state x particle x time x samples x chains if `n_particles > 1`
#'   * `observations$state`: the final state, laid out as trajectories
#'     without the time dimension
#'
#'   Everything needed to tidy the result or continue the run is stored
#'   as an attribute, so [tidy_trajectories()] and
#'   [continue_parameter_sets()] need no further information.
#'
#' @export
run_parameter_sets <- function(generator, pars, times, packer = NULL,
                               fixed = NULL, modify = NULL, n_particles = 1,
                               deterministic = FALSE, n_draws = NULL,
                               time_start = 0, save_trajectories = TRUE,
                               dt = NULL, ode_control = NULL, seed = NULL,
                               n_threads = 1) {
  p <- prepare_pars(pars, packer)
  spec <- list(generator = generator,
               packer = p$packer,
               fixed = c(p$fixed, check_fixed(fixed, p$matrix)),
               modify = check_modify(modify),
               dt = dt,
               ode_control = ode_control,
               deterministic = deterministic)
  simulate_sets(spec, p$matrix, c(p$n_samples, p$n_chains), times,
                state = NULL, time_start = time_start,
                n_particles = n_particles, n_draws = n_draws,
                save_trajectories = save_trajectories, seed = seed,
                n_threads = n_threads)
}


#' Continue a simulation from the end of a fit or a previous run
#'
#' Run a model forward from the final state and time of a set of
#' samples, with each sample continuing from its own state using its own
#' parameters. The samples can be a fit from [monty::monty_sample()], for
#' projecting forward from the end of the data, or the result of
#' [run_parameter_sets()] or of this function.
#'
#' For a fit, the final state must have been saved by creating the
#' likelihood with `save_state = TRUE` (see
#' [dust2::dust_likelihood_monty()]); it cannot be recovered afterwards.
#' The filter is needed too, as it holds the model, its time step and
#' the time of the last data point.
#'
#' To start from somewhere other than the end, pass `state` and `time`
#' explicitly, for example a snapshot saved with `save_snapshots`.
#'
#' @param samples A `monty_samples` object, from [monty::monty_sample()],
#'   [run_parameter_sets()] or `continue_parameter_sets()`.
#'
#' @param times Times to save output at, increasing and no earlier than
#'   the time being continued from.
#'
#' @param obj For a fit, the filter or unfilter the likelihood was
#'   created from. Not needed for output from [run_parameter_sets()].
#'
#' @param packer For a fit, the packer the likelihood was created with.
#'   Not needed for output from [run_parameter_sets()].
#'
#' @param fixed Optional named list of parameters to hold constant from
#'   this point on, replacing fixed values used so far.
#'
#' @param modify Optional function applied to each parameter list from
#'   this point on, in addition to any applied so far.
#'
#' @param n_particles Number of particles to run for each parameter
#'   set, each starting from that set's final state. If the samples
#'   already hold several particles per set, each particle continues
#'   from its own state and this must not be changed.
#'
#' @param deterministic Logical, indicating that a stochastic model
#'   should be run deterministically. By default this follows the fit
#'   (`TRUE` for an unfilter, `FALSE` for a filter) or the previous run.
#'
#' @param state,time Optional state and time to continue from, instead
#'   of the final state. `state` has one column per sample, so a saved
#'   snapshot can be passed directly:
#'   `samples$observations$snapshots[, 1, , ]`.
#'
#' @param save_trajectories Either `TRUE`, to save all state variables,
#'   or a character vector of variable names. By default the same
#'   variables as the previous run, or all of them for a fit.
#'
#' @param dt Time step, by default that of the fit or previous run.
#'
#' @param seed,n_threads Passed to [dust2::dust_system_create()].
#'
#' @inheritParams run_parameter_sets
#'
#' @return An object of class `monty_samples`, as for
#'   [run_parameter_sets()].
#'
#' @export
continue_parameter_sets <- function(samples, times, obj = NULL,
                                    packer = NULL, fixed = NULL,
                                    modify = NULL, n_particles = NULL,
                                    deterministic = NULL, n_draws = NULL,
                                    state = NULL, time = NULL,
                                    save_trajectories = NULL, dt = NULL,
                                    seed = NULL, n_threads = 1) {
  spec <- attr(samples, "tidyodin")$spec
  meta <- attr(samples, "tidyodin")
  samples <- check_samples(samples)
  shape <- dim(samples$density)
  m <- samples$pars
  dim(m) <- c(nrow(m), prod(shape))
  rownames(m) <- rownames(samples$pars)

  if (is.null(spec)) {
    if (!inherits(obj, "dust_likelihood")) {
      cli::cli_abort(
        c("{.arg obj} is needed to continue from these samples",
          i = "Pass the filter that was used to create the likelihood"))
    }
    if (isTRUE(obj$preserve_group_dimension)) {
      cli::cli_abort("Continuing from a grouped fit is not supported")
    }
    if (is.null(packer)) {
      cli::cli_abort(
        "{.arg packer} is needed to continue from a fit")
    }
    check_packer(packer, m)
    spec <- list(generator = obj$generator,
                 packer = packer,
                 fixed = NULL,
                 modify = list(),
                 dt = obj$time_control$dt,
                 ode_control = obj$time_control$ode_control,
                 deterministic = obj$deterministic)
    time <- time %||% max(obj$inputs$time)
    n_particles_prev <- 1
  } else {
    time <- time %||% utils::tail(meta$times, 1)
    save_trajectories <- save_trajectories %||%
      unique(parse_names(meta$names_trajectories)$variable)
    n_particles_prev <- if (is.null(state)) meta$n_particles else 1
  }

  spec$fixed <- utils::modifyList(as.list(spec$fixed),
                                  as.list(check_fixed(fixed, m)))
  spec$modify <- c(spec$modify, check_modify(modify))
  spec$dt <- dt %||% spec$dt
  spec$deterministic <- deterministic %||% spec$deterministic

  if (is.null(state)) {
    state <- samples$observations$state
    if (is.null(state)) {
      cli::cli_abort(
        c("These samples do not contain a final state to continue from",
          i = "Create the likelihood with {.code save_state = TRUE} in {.fn dust2::dust_likelihood_monty}, or pass {.arg state} and {.arg time}"))
    }
  }
  if (n_particles_prev > 1) {
    if (!is.null(n_particles) && n_particles != n_particles_prev) {
      cli::cli_abort(
        "These samples hold {n_particles_prev} particles per parameter set, so {.arg n_particles} can't be changed to {n_particles}")
    }
    n_particles <- n_particles_prev
  }

  simulate_sets(spec, m, shape, times, state = state, time_start = time,
                n_particles = n_particles %||% 1, n_draws = n_draws,
                save_trajectories = save_trajectories %||% TRUE, seed = seed,
                n_threads = n_threads,
                state_has_particles = n_particles_prev > 1)
}


## The shared engine. `m` is a matrix of parameters with one column per
## run and `shape` is c(n_samples, n_chains). `spec` holds everything
## needed to turn a column of `m` into the parameters for the model; it
## is stored on the result so that a run can be continued, and it stays
## valid if the result is thinned because parameters are always rebuilt
## from `pars`.
simulate_sets <- function(spec, m, shape, times, state, time_start,
                          n_particles, n_draws, save_trajectories, seed,
                          n_threads, state_has_particles = FALSE,
                          call = rlang::caller_env()) {
  n_total <- ncol(m)
  draws <- seq_len(n_total)
  if (!is.null(n_draws) && n_draws < n_total) {
    draws <- sort(sample.int(n_total, n_draws))
    m <- m[, draws, drop = FALSE]
    shape <- c(length(draws), 1L)
  }
  n_groups <- ncol(m)

  pars_list <- lapply(seq_len(n_groups), function(i) {
    p <- if (is.null(spec$packer)) list() else spec$packer$unpack(unname(m[, i]))
    p <- utils::modifyList(p, as.list(spec$fixed))
    for (f in spec$modify) {
      p <- f(p)
    }
    p
  })

  sys <- dust2::dust_system_create(
    spec$generator, pars_list, n_particles = n_particles,
    n_groups = n_groups, time = time_start, dt = spec$dt,
    ode_control = spec$ode_control, seed = seed,
    deterministic = spec$deterministic, n_threads = n_threads,
    preserve_particle_dimension = TRUE, preserve_group_dimension = TRUE)
  names_state <- sys$packer_state$names()
  n_state <- length(names_state)

  if (is.null(state)) {
    dust2::dust_system_set_state_initial(sys)
  } else {
    n_each <- if (state_has_particles) n_particles else 1
    expected <- n_state * n_each * n_total
    if (length(state) != expected || NROW(state) != n_state) {
      cli::cli_abort(
        c("Unexpected size for the state to continue from",
          x = "Expected {n_state} state variable{?s} for each of {n_total} parameter set{?s}",
          x = "Found {.val {paste(dim(state) %||% length(state), collapse = ' x ')}}"),
        call = call)
    }
    ## state x particle x group, or state x group repeated over particles
    state <- array(state, c(n_state * n_each, n_total))[, draws, drop = FALSE]
    if (!state_has_particles) {
      state <- state[, rep(seq_len(n_groups), each = n_particles),
                     drop = FALSE]
    }
    dust2::dust_system_set_state(
      sys, array(state, c(n_state, n_particles, n_groups)))
  }

  index_state <- NULL
  names_trajectories <- names_state
  if (is.character(save_trajectories)) {
    index_state <- sys$packer_state$subset(save_trajectories)$index
    names_trajectories <- names_state[index_state]
  } else if (!isTRUE(save_trajectories)) {
    cli::cli_abort(
      "{.arg save_trajectories} must be {.code TRUE} or a character vector",
      call = call)
  }

  ## state x particle x group x time
  y <- dust2::dust_system_simulate(sys, times, index_state = index_state)
  y <- aperm(y, c(1, 2, 4, 3))
  particle <- if (n_particles > 1) n_particles else NULL
  dim(y) <- c(length(names_trajectories), particle, length(times), shape)

  state_end <- dust2::dust_system_state(sys)
  dim(state_end) <- c(n_state, particle, shape)

  pars_array <- array(m, c(nrow(m), shape),
                      dimnames = list(rownames(m), NULL, NULL))
  initial <- array(pars_array[, 1, ], c(nrow(m), shape[[2]]),
                   dimnames = list(rownames(m), NULL))

  ret <- list(pars = pars_array,
              density = matrix(NA_real_, shape[[1]], shape[[2]]),
              initial = initial,
              details = NULL,
              state = NULL,
              observations = list(trajectories = y, state = state_end),
              full_chains = NULL)
  class(ret) <- "monty_samples"
  attr(ret, "tidyodin") <- list(names_trajectories = names_trajectories,
                                names_state = names_state,
                                times = times,
                                n_particles = n_particles,
                                spec = spec)
  ret
}


## Convert any supported form of `pars` into a matrix with one column
## per run, and a packer that turns a column into a parameter list.
prepare_pars <- function(pars, packer, call = rlang::caller_env()) {
  need_packer <- function(what) {
    if (is.null(packer)) {
      cli::cli_abort("{.arg packer} is required when {.arg pars} is {what}",
                     call = call)
    }
  }
  if (inherits(packer, "monty_packer_grouped")) {
    cli::cli_abort("Grouped packers are not supported", call = call)
  }
  n_chains <- 1L
  fixed <- NULL

  if (inherits(pars, "monty_samples")) {
    need_packer("a {.cls monty_samples} object")
    pars <- check_samples(pars, call)
    n_chains <- dim(pars$pars)[[3]]
    m <- pars$pars
    dim(m) <- c(nrow(m), prod(dim(m)[-1]))
    rownames(m) <- rownames(pars$pars)
    check_packer(packer, m, call)
  } else if (is.data.frame(pars)) {
    if (nrow(pars) == 0) {
      cli::cli_abort("{.arg pars} has no rows", call = call)
    }
    if (is.null(packer)) {
      ok <- vapply(pars, is.numeric, TRUE)
      if (!all(ok)) {
        cli::cli_abort(
          c("All columns of {.arg pars} must be numeric",
            x = "Not numeric: {.val {names(pars)[!ok]}}",
            i = "Use {.arg fixed} for parameters that do not vary"),
          call = call)
      }
      packer <- monty::monty_packer(names(pars))
    }
    missing <- setdiff(packer$names(), names(pars))
    if (length(missing) > 0) {
      cli::cli_abort("Columns missing from {.arg pars}: {.val {missing}}",
                     call = call)
    }
    m <- t(as.matrix(pars[packer$names()]))
  } else if (is.matrix(pars)) {
    need_packer("a matrix")
    m <- pars
    check_packer(packer, m, call)
    rownames(m) <- packer$names()
  } else if (is.list(pars) && length(pars) > 0 &&
               all(vapply(pars, is.list, TRUE))) {
    pars <- unname(pars)
    if (is.null(packer)) {
      auto <- packer_from_list(pars, call)
      packer <- auto$packer
      fixed <- auto$fixed
      pars <- lapply(pars, function(el) el[setdiff(names(el), names(fixed))])
    }
    if (is.null(packer)) {
      m <- matrix(numeric(), 0, length(pars))
    } else {
      m <- vapply(pars, packer$pack, numeric(length(packer$names())))
      m <- matrix(m, length(packer$names()), length(pars),
                  dimnames = list(packer$names(), NULL))
    }
  } else {
    cli::cli_abort(
      "{.arg pars} must be a {.cls monty_samples} object, a data frame, a list of parameter lists or a matrix",
      call = call)
  }

  dimnames(m) <- list(rownames(m), NULL)
  list(matrix = m, packer = packer, fixed = fixed,
       n_samples = ncol(m) %/% n_chains, n_chains = n_chains)
}


## Work out a packer for a list of parameter lists: elements that
## differ between sets are the parameters (if nothing differs, every
## numeric element is), and everything else is fixed.
packer_from_list <- function(pars, call) {
  nms <- names(pars[[1]]) %||% character()
  same_names <- vapply(pars, function(p) {
    setequal(names(p) %||% character(), nms)
  }, TRUE)
  if (!all(same_names) || anyDuplicated(nms) || any(nms == "")) {
    cli::cli_abort(
      c("Every parameter set in {.arg pars} must have the same named elements",
        i = "Use {.arg fixed} for parameters that do not vary"),
      call = call)
  }
  first <- pars[[1]]
  varies <- vapply(nms, function(nm) {
    !all(vapply(pars, function(p) identical(p[[nm]], first[[nm]]), TRUE))
  }, TRUE)
  numeric <- vapply(first[nms], is.numeric, TRUE)
  if (!any(varies)) {
    varies <- numeric
  }
  if (any(varies & !numeric)) {
    cli::cli_abort(
      c("Parameters that vary between sets must be numeric",
        x = "Not numeric: {.val {nms[varies & !numeric]}}",
        i = "Provide a {.arg packer} to handle these"),
      call = call)
  }
  fixed <- first[nms[!varies]]
  if (!any(varies)) {
    return(list(packer = NULL, fixed = fixed))
  }
  shape <- lapply(first[nms[varies]], function(el) {
    if (length(el) == 1 && is.null(dim(el))) integer() else
      dim(el) %||% length(el)
  })
  list(packer = monty::monty_packer(array = shape), fixed = fixed)
}


check_packer <- function(packer, m, call = rlang::caller_env()) {
  if (inherits(packer, "monty_packer_grouped")) {
    cli::cli_abort("Grouped packers are not supported", call = call)
  }
  if (!inherits(packer, "monty_packer")) {
    cli::cli_abort("{.arg packer} must be a {.cls monty_packer}",
                   call = call)
  }
  nms <- packer$names()
  if (nrow(m) != length(nms)) {
    cli::cli_abort(
      "Expected {length(nms)} parameter{?s} to match {.arg packer}, but found {nrow(m)}",
      call = call)
  }
  if (!is.null(rownames(m)) && !identical(rownames(m), nms)) {
    cli::cli_abort(
      c("Parameter names do not match {.arg packer}",
        x = "Expected: {.val {nms}}",
        x = "Found: {.val {rownames(m)}}"),
      call = call)
  }
}


check_fixed <- function(fixed, m, call = rlang::caller_env()) {
  if (is.null(fixed)) {
    return(NULL)
  }
  nms <- names(fixed)
  if (!is.list(fixed) || is.null(nms) || any(nms == "") ||
        anyDuplicated(nms)) {
    cli::cli_abort("{.arg fixed} must be a named list", call = call)
  }
  clash <- intersect(nms, parse_names(rownames(m) %||% character())$variable)
  if (length(clash) > 0) {
    cli::cli_abort(
      c("{.arg fixed} names parameters that vary in {.arg pars}: {.val {clash}}",
        i = "Remove them from one or the other, or use {.arg modify}"),
      call = call)
  }
  fixed
}


check_modify <- function(modify, call = rlang::caller_env()) {
  if (is.null(modify)) {
    return(list())
  }
  if (!is.function(modify)) {
    cli::cli_abort("{.arg modify} must be a function", call = call)
  }
  list(modify)
}
