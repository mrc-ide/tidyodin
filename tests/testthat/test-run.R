test_that("run over a list of parameter sets", {
  pars <- lapply(c(0.2, 0.4, 0.8), function(b) list(beta = b, gamma = 0.1))
  t <- seq(0, 20, by = 2)
  res <- run_parameter_sets(model("ode"), pars, t)

  expect_s3_class(res, "monty_samples")
  expect_named(res, c("pars", "density", "initial", "details", "state",
                      "observations", "full_chains"))
  ## gamma is the same in every set, so only beta is a parameter
  expect_equal(dim(res$pars), c(1, 3, 1))
  expect_equal(rownames(res$pars), "beta")
  expect_equal(res$pars["beta", , 1], c(0.2, 0.4, 0.8))
  expect_equal(dim(res$density), c(3, 1))
  expect_equal(dim(res$initial), c(1, 1))
  expect_equal(dim(res$observations$trajectories), c(2, 11, 3, 1))
  expect_equal(dim(res$observations$state), c(2, 3, 1))

  sys <- dust2::dust_system_create(model("ode"), pars[[2]])
  dust2::dust_system_set_state_initial(sys)
  y <- dust2::dust_system_simulate(sys, t)
  expect_equal(res$observations$trajectories[, , 2, 1], y)

  tidy <- tidy_trajectories(res)
  expect_named(tidy, c("chain", "iteration", "sample", "time", "variable",
                       "value"))
  expect_equal(tidy$value[tidy$variable == "I" & tidy$sample == 2], y[2, ])
  expect_equal(nrow(tidy_pars(res)), 3)
  expect_equal(nrow(tidy_state(res)), 6)
})

test_that("fixed parameters are used but not recorded as parameters", {
  t <- 0:10
  grid <- data.frame(beta = c(0.2, 0.4, 0.8))
  res <- run_parameter_sets(model("ode"), grid, t,
                            fixed = list(gamma = 0.3, N = 500))
  expect_equal(rownames(res$pars), "beta")

  sys <- dust2::dust_system_create(
    model("ode"), list(beta = 0.8, gamma = 0.3, N = 500))
  dust2::dust_system_set_state_initial(sys)
  expect_equal(res$observations$trajectories[, , 3, 1],
               dust2::dust_system_simulate(sys, t))

  ## fixed replaces values fixed in a packer
  packer <- monty::monty_packer("beta", fixed = list(gamma = 0.1, N = 1000))
  res2 <- run_parameter_sets(model("ode"), grid, t, packer = packer,
                             fixed = list(gamma = 0.3, N = 500))
  expect_equal(res2$observations, res$observations)

  expect_error(
    run_parameter_sets(model("ode"), grid, t, fixed = list(beta = 1)),
    "names parameters that vary")
  expect_error(run_parameter_sets(model("ode"), grid, t, fixed = list(1)),
               "named list")
  expect_error(
    run_parameter_sets(model("ode"), data.frame(beta = 1, x = "a"), t),
    "must be numeric")
})

test_that("lists of parameter sets can hold arrays and fixed elements", {
  pars <- list(list(N = 1000, beta = 0.2, n_age = 2, n_vax = 3),
               list(N = 1000, beta = 0.6, n_age = 2, n_vax = 3))
  res <- run_parameter_sets(model("age"), pars, 0:4, seed = 1)
  expect_equal(rownames(res$pars), "beta")
  expect_equal(dim(res$observations$trajectories), c(15, 5, 2, 1))
  tidy <- tidy_trajectories(res, dims = c("age", "vax"))
  expect_true(all(c("age", "vax") %in% names(tidy)))

  packer <- monty::monty_packer(array = list(x = 2))
  p <- prepare_pars(
    list(list(x = c(1, 2), k = "a"), list(x = c(3, 4), k = "a")), NULL)
  expect_equal(p$packer$names(), packer$names())
  expect_equal(p$matrix,
               matrix(1:4, 2, dimnames = list(c("x[1]", "x[2]"), NULL)))
  expect_equal(p$fixed, list(k = "a"))

  expect_error(prepare_pars(list(list(a = 1), list(b = 1)), NULL),
               "same named elements")
  expect_error(prepare_pars(list(list(a = "x"), list(a = "y")), NULL),
               "must be numeric")
})

test_that("can run a single set, or a model with no parameters given", {
  res <- run_parameter_sets(model("sir"), list(list()), 0:5, n_particles = 6,
                            dt = 0.25)
  expect_equal(dim(res$pars), c(0, 1, 1))
  expect_equal(dim(res$observations$trajectories), c(4, 6, 6, 1, 1))
  expect_equal(nrow(tidy_pars(res)), 0)
  expect_equal(nrow(tidy_trajectories(res)), 4 * 6 * 6)
  expect_message(print(res), "0 parameters x 1 sample")

  res <- run_parameter_sets(model("ode"), list(list(beta = 0.3, gamma = 0.2)),
                            0:5)
  expect_equal(rownames(res$pars), c("beta", "gamma"))
})

test_that("output works with monty's own tools", {
  pars <- data.frame(beta = seq(0.2, 0.9, length.out = 8), gamma = 0.1)
  res <- run_parameter_sets(model("ode"), pars, 0:5)
  expect_message(print(res), "2 parameters x 8 samples x 1 chain")
  thin <- monty::monty_samples_thin(res, thinning_factor = 2)
  expect_equal(dim(thin$observations$trajectories), c(2, 6, 4, 1))
  expect_equal(nrow(tidy_trajectories(thin)), 2 * 6 * 4)
  flat <- monty::monty_flatten_chains(res)
  expect_equal(tidy_trajectories(flat), tidy_trajectories(res))
  skip_if_not_installed("posterior")
  expect_equal(posterior::as_draws_df(res)$beta, pars$beta)
})

test_that("run from samples keeps chains and matches the fit layout", {
  fit <- sir_fit(save_trajectories = TRUE, save_state = TRUE,
                 save_snapshots = c(4, 8))
  s <- fit$samples
  res <- run_parameter_sets(model("sir"), s, 0:10, packer = sir_packer(),
                            dt = 0.25, seed = 1)
  expect_equal(res$pars, s$pars)
  expect_equal(dim(res$observations$trajectories), c(4, 11, 30, 3))
  expect_equal(names(tidy_trajectories(res)),
               names(tidy_trajectories(s, fit$filter)))

  ## each run uses the right parameters: check against a direct run
  p <- lapply(seq_len(90), function(i) {
    sir_packer()$unpack(matrix(s$pars, 2)[, i])
  })
  sys <- dust2::dust_system_create(model("sir"), p, n_groups = 90, dt = 0.25,
                                   seed = 1)
  dust2::dust_system_set_state_initial(sys)
  y <- dust2::dust_system_simulate(sys, 0:10)
  expect_equal(res$observations$trajectories[, , 13, 2], y[, 43, ])

  expect_error(run_parameter_sets(model("sir"), s, 0:10), "packer")
  expect_error(
    run_parameter_sets(model("sir"), s, 0:10,
                       packer = monty::monty_packer(c("gamma", "beta"))),
    "names do not match")
})

test_that("stochastic models can be run deterministically", {
  grid <- data.frame(beta = c(0.2, 0.5))
  a <- run_parameter_sets(model("sir"), grid, 0:10, deterministic = TRUE,
                          dt = 0.25)
  b <- run_parameter_sets(model("sir"), grid, 0:10, deterministic = TRUE,
                          dt = 0.25)
  expect_equal(a$observations, b$observations)
  ## expectations rather than draws, so not whole numbers
  expect_false(all(a$observations$trajectories %% 1 == 0))
  s <- run_parameter_sets(model("sir"), grid, 0:10, dt = 0.25)
  expect_true(all(s$observations$trajectories %% 1 == 0))
})

test_that("can rerun a fit with changed parameters", {
  skip_if_not_installed("posterior")
  fit <- sir_fit(save_trajectories = TRUE, save_state = TRUE,
                 save_snapshots = c(4, 8))
  s <- fit$samples
  pars <- posterior::as_draws_df(s)
  pars$beta <- 0
  none <- run_parameter_sets(fit$filter, pars, 0:5, packer = sir_packer())
  expect_true(all(none$observations$trajectories[4, , , ] == 0))
  ## chains are kept, and parameters are as supplied
  expect_equal(dim(none$pars), dim(s$pars))
  expect_equal(none$pars["gamma", , ], s$pars["gamma", , ])
  expect_true(all(none$pars["beta", , ] == 0))

  ## without a packer, other parameters come from fixed
  none2 <- run_parameter_sets(fit$filter, pars, 0:5,
                              fixed = list(N = 1000, I0 = 10), seed = 1)
  expect_equal(rownames(none2$pars), c("beta", "gamma"))
  expect_equal(dim(none2$pars), dim(s$pars))

  ## draws that are not whole chains become a single chain
  sub <- run_parameter_sets(fit$filter, pars[c(1, 40, 41), ], 0:5,
                            packer = sir_packer())
  expect_equal(dim(sub$pars), c(2, 3, 1))
})

test_that("model settings are taken from a filter", {
  fit <- sir_fit(save_trajectories = TRUE, save_state = TRUE,
                 save_snapshots = c(4, 8))
  s <- fit$samples
  a <- run_parameter_sets(fit$filter, s, 0:10, packer = sir_packer(),
                          seed = 1)
  b <- run_parameter_sets(model("sir"), s, 0:10, packer = sir_packer(),
                          dt = 0.25, seed = 1)
  expect_equal(a$observations, b$observations)
  expect_equal(attr(a, "tidyodin")$spec$dt, 0.25)
  ## and can be overridden
  d <- run_parameter_sets(fit$filter, s, 0:10, packer = sir_packer(),
                          dt = 1, deterministic = TRUE)
  expect_equal(attr(d, "tidyodin")$spec$dt, 1)
  expect_false(all(d$observations$trajectories %% 1 == 0))

  grouped <- sir_fit_grouped(save_trajectories = TRUE, save_state = TRUE,
                             save_snapshots = c(4, 8))
  expect_error(run_parameter_sets(grouped$filter, list(list()), 0:5),
               "grouped fit is not supported")
})

test_that("can subset draws and saved variables", {
  fit <- sir_fit(save_trajectories = TRUE, save_state = TRUE,
                 save_snapshots = c(4, 8))
  s <- fit$samples
  res <- run_parameter_sets(fit$filter, s, 0:5, packer = sir_packer(),
                            n_samples = 5, save_trajectories = c("S", "R"))
  expect_equal(dim(res$pars), c(2, 5, 1))
  expect_equal(dim(res$observations$trajectories), c(2, 6, 5, 1))
  expect_setequal(tidy_trajectories(res)$variable, c("S", "R"))
  expect_setequal(tidy_state(res)$variable, c("S", "I", "R", "incidence"))
})

test_that("can run several particles per parameter set", {
  pars <- matrix(c(0.2, 0.1, 0.5, 0.1), 2)
  res <- run_parameter_sets(model("sir"), pars, 0:6, packer = sir_packer(),
                            n_particles = 5, dt = 0.25, seed = 1)
  expect_equal(dim(res$observations$trajectories), c(4, 5, 7, 2, 1))
  expect_equal(dim(res$observations$state), c(4, 5, 2, 1))
  tidy <- tidy_trajectories(res)
  expect_named(tidy, c("chain", "iteration", "sample", "particle", "time",
                       "variable", "value"))
  r <- tidy[tidy$variable == "I" & tidy$particle == 4 & tidy$sample == 2, ]
  expect_equal(r$value, res$observations$trajectories[2, 4, , 2, 1])
  expect_equal(res$observations$trajectories[1, , 1, , ], matrix(990, 5, 2))
  s <- summarise_trajectories(tidy)
  expect_equal(nrow(s), 4 * 7)

  expect_error(run_parameter_sets(model("sir"), pars, 0:6), "packer")
  expect_error(run_parameter_sets(model("sir"), 1, 0:6), "must be a")
})

test_that("sample columns from this package are not treated as parameters", {
  pars <- data.frame(chain = rep(1:2, each = 2), iteration = rep(1:2, 2),
                     sample = 1:4, beta = c(0.2, 0.3, 0.4, 0.5))
  res <- run_parameter_sets(model("ode"), pars, 0:3)
  expect_equal(rownames(res$pars), "beta")
  expect_equal(dim(res$pars), c(1, 2, 2))
  expect_named(tidy_trajectories(res),
               c("chain", "iteration", "sample", "time", "variable", "value"))
})

test_that("invalid parameter sets give helpful errors", {
  packer <- sir_packer()
  expect_error(
    run_parameter_sets(model("sir"), data.frame(beta = numeric()), 0:3),
    "has no rows")
  expect_error(
    run_parameter_sets(model("sir"), data.frame(beta = 0.2), 0:3,
                       packer = packer),
    "Columns missing from `pars`: \"gamma\"")
  expect_error(
    run_parameter_sets(model("sir"), matrix(0.2, 3, 2), 0:3, packer = packer),
    "Expected 2 parameters to match `packer`, but found 3")
  expect_error(
    run_parameter_sets(model("sir"), matrix(0.2, 2, 2), 0:3,
                       packer = list()),
    "must be a <monty_packer>")
  grouped <- monty::monty_packer_grouped(c("a", "b"), c("beta", "gamma"))
  expect_error(
    run_parameter_sets(model("sir"), matrix(0.2, 4, 2), 0:3,
                       packer = grouped),
    "Grouped packers are not supported")
  expect_error(
    run_parameter_sets(model("sir"), data.frame(beta = 0.2), 0:3,
                       save_trajectories = FALSE),
    "must be `TRUE` or a character vector")
})

test_that("chains are only kept when samples divide evenly into them", {
  expect_equal(chains_from_draws(NULL), 1)
  expect_equal(chains_from_draws(c(1, 1, 1)), 1)
  expect_equal(chains_from_draws(c(1, 1, 2, 2, 3, 3)), 3)
  ## uneven, or not in order, so there is no chain structure to keep
  expect_equal(chains_from_draws(c(1, 2, 2)), 1)
  expect_equal(chains_from_draws(c(1, 2, 1, 2)), 1)
  expect_equal(chains_from_draws(c(2, 2, 3, 3)), 1)

  pars <- data.frame(chain = c(1, 2, 2), beta = c(0.2, 0.3, 0.4))
  res <- run_parameter_sets(model("ode"), pars, 0:2)
  expect_equal(dim(res$pars), c(1, 3, 1))
})
