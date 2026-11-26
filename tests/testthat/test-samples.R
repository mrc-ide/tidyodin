test_that("tidy parameters and density", {
  fit <- sir_fit(save_trajectories = TRUE, save_state = TRUE,
                 save_snapshots = c(4, 8))
  s <- fit$samples
  res <- tidy_pars(s)
  expect_named(res, c("chain", "iteration", "sample", "parameter", "value"))
  expect_equal(nrow(res), 2 * 30 * 3)
  r <- res[res$parameter == "gamma" & res$chain == 2, ]
  expect_equal(r$value, s$pars["gamma", , 2])
  expect_equal(r$iteration, 1:30)
  expect_equal(r$sample, 31:60)

  dens <- tidy_density(s)
  expect_equal(dens$density[dens$chain == 3], s$density[, 3])
  expect_equal(dens$sample, 1:90)

  skip_if_not_installed("posterior")
  df <- posterior::as_draws_df(s)
  expect_equal(res$value[res$parameter == "beta"], df$beta)
  expect_equal(res$sample[res$parameter == "beta"], df$.draw)
  expect_equal(res$chain[res$parameter == "beta"], df$.chain)
})

test_that("tidy array and grouped parameters", {
  nms <- c("gamma", "beta<a>", "x[1]<a>", "x[2]<a>", "beta<b>", "x[1]<b>",
           "x[2]<b>")
  pars <- array(runif(7 * 5 * 2), c(7, 5, 2), dimnames = list(nms, NULL, NULL))
  s <- structure(list(pars = pars, density = matrix(0, 5, 2)),
                 class = "monty_samples")
  res <- tidy_pars(s, dims = "strain")
  expect_named(res, c("chain", "iteration", "sample", "parameter", "strain",
                      "group", "value"))
  r <- res[res$parameter == "x" & res$strain %in% 2 & res$group %in% "b", ]
  expect_equal(r$value, c(pars["x[2]<b>", , ]))
  expect_true(all(is.na(res$group[res$parameter == "gamma"])))
  expect_error(tidy_pars(list()), "monty_samples")
})

test_that("tidy trajectories, state and snapshots from a fit", {
  fit <- sir_fit(save_trajectories = TRUE, save_state = TRUE,
                 save_snapshots = c(4, 8))
  s <- fit$samples
  traj <- dust2::dust_unpack_state(fit$filter, s$observations$trajectories)

  res <- tidy_trajectories(s, fit$filter)
  expect_named(res, c("chain", "iteration", "sample", "time", "variable",
                      "value"))
  expect_equal(nrow(res), 4 * 10 * 30 * 3)
  r <- res[res$variable == "incidence" & res$chain == 2 &
             res$iteration == 17, ]
  expect_equal(r$value, traj$incidence[, 17, 2])
  expect_equal(r$time, 1:10)
  expect_equal(unique(r$sample), 47)

  res <- tidy_state(s, fit$filter)
  expect_named(res, c("chain", "iteration", "sample", "variable", "value"))
  r <- res[res$variable == "R" & res$chain == 3, ]
  expect_equal(r$value, s$observations$state[3, , 3])

  res <- tidy_snapshots(s, fit$filter, times = c(4, 8))
  r <- res[res$variable == "I" & res$time == 8 & res$chain == 1, ]
  expect_equal(r$value, s$observations$snapshots[2, 2, , 1])
  expect_error(tidy_snapshots(s, fit$filter), "times")

  expect_error(tidy_trajectories(s), "is needed")
  expect_error(tidy_trajectories(s, "filter"), "must be a")
})

test_that("tidy trajectories from flattened and thinned samples", {
  fit <- sir_fit(save_trajectories = TRUE, save_state = TRUE,
                 save_snapshots = c(4, 8))
  s <- fit$samples
  expect_equal(tidy_trajectories(monty::monty_flatten_chains(s), fit$filter),
               tidy_trajectories(s, fit$filter))
  expect_equal(tidy_pars(monty::monty_flatten_chains(s)), tidy_pars(s))

  thin <- monty::monty_samples_thin(s, thinning_factor = 2, burnin = 10)
  res <- tidy_trajectories(thin, fit$filter, vars = "I")
  expect_equal(nrow(res), 10 * 10 * 3)
  r <- res[res$chain == 2 & res$iteration == 4, ]
  expect_equal(r$value, thin$observations$trajectories[2, , 4, 2])
})

test_that("can subsample draws", {
  fit <- sir_fit(save_trajectories = TRUE, save_state = TRUE,
                 save_snapshots = c(4, 8))
  s <- fit$samples
  set.seed(1)
  res <- tidy_trajectories(s, fit$filter, n_draws = 7, vars = "I")
  expect_equal(nrow(res), 7 * 10)
  expect_length(unique(res$sample), 7)
  r <- res[res$sample == res$sample[[nrow(res)]], ]
  expect_equal(r$value,
               s$observations$trajectories[2, , r$iteration[[1]],
                                           r$chain[[1]]])
})

test_that("tidy trajectories saved for a subset of variables", {
  fit <- sir_fit(save_trajectories = c("I", "incidence"))
  s <- fit$samples
  expect_error(tidy_trajectories(s, fit$filter),
               "Unexpected array dimensions")
  res <- tidy_trajectories(s, fit$filter,
                           save_trajectories = c("I", "incidence"))
  expect_setequal(res$variable, c("I", "incidence"))
  r <- res[res$variable == "incidence" & res$sample == 5, ]
  expect_equal(r$value, s$observations$trajectories[2, , 5, 1])
  expect_error(tidy_state(s, fit$filter), "do not contain any state")
})

test_that("summarise trajectories", {
  data <- tibble::tibble(sample = rep(1:100, 2), time = rep(1:2, each = 100),
                         variable = "x", value = c(1:100, 101:200))
  res <- summarise_trajectories(data, interval = 0.5)
  expect_named(res, c("time", "variable", "mean", "median", "lower", "upper"))
  expect_equal(res$median, c(50.5, 150.5))
  expect_equal(res$lower, unname(c(quantile(1:100, 0.25),
                                   quantile(101:200, 0.25))))
  expect_error(summarise_trajectories(1), "data frame")
})
