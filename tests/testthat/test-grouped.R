test_that("tidy parameters from a grouped fit", {
  fit <- sir_fit_grouped(save_trajectories = TRUE, save_state = TRUE,
                         save_snapshots = c(4, 8))
  s <- fit$samples
  expect_equal(rownames(s$pars), c("gamma", "beta<north>", "beta<south>"))

  res <- tidy_pars(s)
  expect_named(res, c(".chain", ".iteration", ".draw", "parameter", "group",
                      "value"))
  expect_equal(nrow(res), 3 * 20 * 2)
  r <- res[res$parameter == "beta" & res$group %in% "south", ]
  expect_equal(r$value, c(s$pars["beta<south>", , ]))
  r <- res[res$parameter == "beta" & res$group %in% "north" &
             res$.chain == 2, ]
  expect_equal(r$value, s$pars["beta<north>", , 2])
  r <- res[res$parameter == "gamma", ]
  expect_equal(r$value, c(s$pars["gamma", , ]))
  expect_true(all(is.na(r$group)))

  res <- tidy_pars(s, vars = "beta")
  expect_equal(unique(res$parameter), "beta")
  expect_setequal(res$group, c("north", "south"))
})

test_that("tidy trajectories from a grouped fit", {
  fit <- sir_fit_grouped(save_trajectories = TRUE, save_state = TRUE,
                         save_snapshots = c(4, 8))
  s <- fit$samples
  traj <- s$observations$trajectories
  expect_equal(dim(traj), c(4, 2, 10, 20, 2))

  res <- tidy_trajectories(s, fit$filter)
  expect_named(res, c(".chain", ".iteration", ".draw", "group", "time",
                      "variable", "value"))
  expect_equal(nrow(res), length(traj))
  r <- res[res$variable == "incidence" & res$group == "south" &
             res$.chain == 2 & res$.iteration == 11, ]
  expect_equal(r$value, traj[4, 2, , 11, 2])
  expect_equal(r$time, 1:10)
  expect_equal(unique(r$.draw), 31)
  r <- res[res$variable == "I" & res$group == "north" & res$time == 6, ]
  expect_equal(r$value, c(traj[2, 1, 6, , ]))

  ## groups differ: south was fitted to twice the cases
  expect_false(identical(traj[, 1, , , ], traj[, 2, , , ]))

  res <- tidy_trajectories(s, fit$filter, vars = "S",
                           labels = list(group = c("south", "north")))
  expect_equal(levels(res$group), c("south", "north"))
  r <- res[res$group == "south" & res$.draw == 3, ]
  expect_equal(r$value, traj[1, 2, , 3, 1])
})

test_that("tidy state and snapshots from a grouped fit", {
  fit <- sir_fit_grouped(save_trajectories = TRUE, save_state = TRUE,
                         save_snapshots = c(4, 8))
  s <- fit$samples

  res <- tidy_state(s, fit$filter)
  expect_named(res, c(".chain", ".iteration", ".draw", "group", "variable",
                      "value"))
  r <- res[res$variable == "R" & res$group == "south" & res$.chain == 1, ]
  expect_equal(r$value, s$observations$state[3, 2, , 1])

  res <- tidy_snapshots(s, fit$filter, times = c(4, 8))
  expect_named(res, c(".chain", ".iteration", ".draw", "group", "time",
                      "variable", "value"))
  r <- res[res$variable == "I" & res$group == "north" & res$time == 8 &
             res$.chain == 2, ]
  expect_equal(r$value, s$observations$snapshots[2, 1, 2, , 2])
})

test_that("grouped fits work when flattened, thinned or subsampled", {
  fit <- sir_fit_grouped(save_trajectories = TRUE, save_state = TRUE,
                         save_snapshots = c(4, 8))
  s <- fit$samples
  expect_equal(tidy_trajectories(monty::monty_flatten_chains(s), fit$filter),
               tidy_trajectories(s, fit$filter))
  expect_equal(tidy_pars(monty::monty_flatten_chains(s)), tidy_pars(s))

  thin <- monty::monty_samples_thin(s, thinning_factor = 4)
  res <- tidy_trajectories(thin, fit$filter, vars = "I")
  expect_equal(nrow(res), 2 * 10 * 5 * 2)
  r <- res[res$group == "south" & res$.chain == 2 & res$.iteration == 3, ]
  expect_equal(r$value, thin$observations$trajectories[2, 2, , 3, 2])

  set.seed(1)
  res <- tidy_trajectories(s, fit$filter, n_draws = 6, vars = "I")
  expect_equal(nrow(res), 2 * 10 * 6)
  r <- res[res$.draw == res$.draw[[1]] & res$group == "north", ]
  expect_equal(r$value,
               s$observations$trajectories[2, 1, , r$.iteration[[1]],
                                           r$.chain[[1]]])

  summary <- summarise_trajectories(tidy_trajectories(s, fit$filter))
  expect_equal(nrow(summary), 2 * 10 * 4)
  expect_setequal(summary$group, c("north", "south"))
})

test_that("grouped fit saving a subset of trajectories", {
  fit <- sir_fit_grouped(save_trajectories = c("I", "incidence"))
  s <- fit$samples
  expect_equal(dim(s$observations$trajectories), c(2, 2, 10, 20, 2))
  res <- tidy_trajectories(s, fit$filter,
                           save_trajectories = c("I", "incidence"))
  expect_setequal(res$variable, c("I", "incidence"))
  r <- res[res$variable == "incidence" & res$group == "south" &
             res$.draw == 25, ]
  expect_equal(r$value, s$observations$trajectories[2, 2, , 5, 2])
})

test_that("run_parameter_sets rejects grouped packers", {
  fit <- sir_fit_grouped(save_trajectories = TRUE, save_state = TRUE,
                         save_snapshots = c(4, 8))
  expect_error(
    run_parameter_sets(model("sir"), fit$samples, 0:5, packer = fit$packer),
    "Grouped packers are not supported")
})
