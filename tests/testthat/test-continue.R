test_that("continuing a deterministic run matches running straight through", {
  grid <- data.frame(beta = c(0.2, 0.4, 0.8))
  fixed <- list(gamma = 0.15)
  full <- run_parameter_sets(model("ode"), grid, 0:20, fixed = fixed)
  first <- run_parameter_sets(model("ode"), grid, 0:10, fixed = fixed)
  second <- continue_parameter_sets(first, 10:20)

  expect_s3_class(second, "monty_samples")
  expect_equal(second$pars, first$pars)
  expect_equal(second$observations$trajectories,
               full$observations$trajectories[, 11:21, , , drop = FALSE],
               tolerance = 1e-6)
  expect_equal(range(tidy_trajectories(second)$time), c(10, 20))

  ## and again, from the result of continuing
  third <- continue_parameter_sets(second, 20:25)
  expect_equal(third$observations$trajectories[, 1, , 1],
               second$observations$state[, , 1], tolerance = 1e-6)
})

test_that("continuing remembers fixed parameters and modifications", {
  grid <- data.frame(beta = c(0.2, 0.5))
  half <- function(p) {
    p$beta <- p$beta / 2
    p
  }
  first <- run_parameter_sets(model("sir"), grid, 0:5, deterministic = TRUE,
                              fixed = list(gamma = 0.3), modify = half,
                              dt = 0.25)
  full <- run_parameter_sets(model("sir"), grid, 0:10, deterministic = TRUE,
                             fixed = list(gamma = 0.3), modify = half,
                             dt = 0.25)
  second <- continue_parameter_sets(first, 5:10, save_trajectories = "I")
  expect_equal(second$observations$trajectories[1, , , ],
               full$observations$trajectories[2, 6:11, , ])

  ## changing things from this point on
  stop <- continue_parameter_sets(first, 5:10, fixed = list(gamma = 0),
                                  modify = function(p) {
                                    p$beta <- 0
                                    p
                                  })
  traj <- stop$observations$trajectories
  expect_equal(traj[2, 6, , 1], traj[2, 1, , 1])
})

test_that("continuing survives thinning", {
  grid <- data.frame(beta = seq(0.2, 0.9, length.out = 8))
  first <- run_parameter_sets(model("ode"), grid, 0:10)
  full <- run_parameter_sets(model("ode"), grid, 0:20)
  thin <- monty::monty_samples_thin(first, thinning_factor = 2)
  second <- continue_parameter_sets(thin, 10:20)
  expect_equal(second$pars["beta", , 1], grid$beta[c(2, 4, 6, 8)])
  expect_equal(second$observations$trajectories[, , , 1],
               full$observations$trajectories[, 11:21, c(2, 4, 6, 8), 1],
               tolerance = 1e-6)
})

test_that("continuing with particles", {
  grid <- data.frame(beta = c(0.2, 0.5))
  one <- run_parameter_sets(model("sir"), grid, 0:5, dt = 0.25, seed = 1)
  many <- continue_parameter_sets(one, 5:8, n_particles = 4, seed = 1)
  expect_equal(dim(many$observations$trajectories), c(4, 4, 4, 2, 1))
  ## every particle starts from its parameter set's final state (S, I, R)
  for (i in 1:4) {
    expect_equal(many$observations$trajectories[1:3, i, 1, , 1],
                 one$observations$state[1:3, , 1])
  }

  ## each particle continues from its own state
  again <- continue_parameter_sets(many, 8:10)
  expect_equal(dim(again$observations$trajectories), c(4, 4, 3, 2, 1))
  expect_equal(again$observations$trajectories[1:3, , 1, , 1],
               many$observations$state[1:3, , , 1])
  expect_named(tidy_trajectories(again),
               c(".chain", ".iteration", ".draw", "particle", "time",
                 "variable", "value"))
  expect_error(continue_parameter_sets(many, 8:10, n_particles = 2),
               "can't be changed")
})

test_that("continue from the end of a fit", {
  fit <- sir_fit(save_trajectories = TRUE, save_state = TRUE,
                 save_snapshots = c(4, 8))
  s <- fit$samples
  res <- continue_parameter_sets(s, 10:15, obj = fit$filter,
                                 packer = sir_packer(),
                                 save_trajectories = c("S", "R"))
  expect_equal(res$pars, s$pars)
  expect_equal(dim(res$observations$trajectories), c(2, 6, 30, 3))
  ## S and R are not reset, so the first time is the saved state
  expect_equal(res$observations$trajectories[, 1, , ],
               s$observations$state[c(1, 3), , ])
  tidy <- tidy_trajectories(res)
  expect_setequal(tidy$variable, c("S", "R"))
  expect_equal(range(tidy$time), c(10, 15))

  ## matches doing it by hand
  p <- lapply(seq_len(90), function(i) {
    sir_packer()$unpack(matrix(s$pars, 2)[, i])
  })
  sys <- dust2::dust_system_create(model("sir"), p, n_groups = 90, time = 10,
                                   dt = 0.25, seed = 1)
  dust2::dust_system_set_state(sys, matrix(s$observations$state, 4))
  y <- dust2::dust_system_simulate(sys, 10:15)
  res <- continue_parameter_sets(s, 10:15, obj = fit$filter,
                                 packer = sir_packer(), seed = 1)
  expect_equal(res$observations$trajectories[, , 7, 3], y[, 67, ])

  many <- continue_parameter_sets(s, 10:12, obj = fit$filter,
                                  packer = sir_packer(), n_particles = 3)
  expect_equal(dim(many$observations$trajectories), c(4, 3, 3, 30, 3))
  expect_equal(many$observations$trajectories[1, 2, 1, , ],
               s$observations$state[1, , ])
})

test_that("continue from a snapshot, and with a subset of draws", {
  fit <- sir_fit(save_trajectories = TRUE, save_state = TRUE,
                 save_snapshots = c(4, 8))
  s <- fit$samples
  res <- continue_parameter_sets(s, 8:10, obj = fit$filter,
                                 packer = sir_packer(),
                                 state = s$observations$snapshots[, 2, , ],
                                 time = 8)
  expect_equal(res$observations$trajectories[c(1, 3), 1, , ],
               s$observations$snapshots[c(1, 3), 2, , ])

  set.seed(1)
  sub <- continue_parameter_sets(s, 10:12, obj = fit$filter,
                                 packer = sir_packer(), n_draws = 5)
  expect_equal(dim(sub$pars), c(2, 5, 1))
  i <- match(sub$pars["beta", 2, 1], matrix(s$pars, 2)[1, ])
  expect_equal(sub$observations$trajectories[1, 1, 2, 1],
               matrix(s$observations$state, 4)[1, i])
})

test_that("continue gives helpful errors", {
  fit <- sir_fit(save_trajectories = TRUE, save_state = TRUE,
                 save_snapshots = c(4, 8))
  s <- fit$samples
  expect_error(continue_parameter_sets(s, 10:12), "is needed to continue")
  expect_error(continue_parameter_sets(s, 10:12, obj = fit$filter),
               "packer")
  expect_error(
    continue_parameter_sets(s, 10:12, obj = fit$filter,
                            packer = sir_packer(), state = matrix(0, 4, 7)),
    "Unexpected size")

  none <- sir_fit(save_trajectories = c("I", "incidence"))
  expect_error(
    continue_parameter_sets(none$samples, 10:12, obj = none$filter,
                            packer = sir_packer()),
    "do not contain a final state")

  grouped <- sir_fit_grouped(save_trajectories = TRUE, save_state = TRUE,
                             save_snapshots = c(4, 8))
  expect_error(
    continue_parameter_sets(grouped$samples, 10:12, obj = grouped$filter,
                            packer = grouped$packer),
    "grouped fit is not supported")
})
