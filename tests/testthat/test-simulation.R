test_that("tidy single-particle simulation", {
  sys <- dust2::dust_system_create(model("ode"), list())
  dust2::dust_system_set_state_initial(sys)
  t <- seq(0, 20, by = 2)
  y <- dust2::dust_system_simulate(sys, t)
  res <- tidy_simulation(y, sys, t)
  expect_named(res, c("time", "variable", "value"))
  expect_equal(res$value[res$variable == "I"], y[2, ])
  expect_equal(res$time[res$variable == "I"], t)
})

test_that("tidy simulation with particles, groups and arrays", {
  pars <- list(list(beta = 0.2), list(beta = 0.4), list(beta = 0.6))
  sys <- dust2::dust_system_create(model("age"), pars, n_particles = 4,
                                   n_groups = 3, seed = 1)
  dust2::dust_system_set_state_initial(sys)
  t <- 0:5
  y <- dust2::dust_system_simulate(sys, t)
  u <- dust2::dust_unpack_state(sys, y)
  res <- tidy_simulation(y, sys, t, dims = c("age", "vax"),
                         labels = list(group = c("low", "mid", "high")))
  expect_named(res, c("group", "particle", "time", "variable", "age", "vax",
                      "value"))
  expect_equal(nrow(res), length(y))
  r <- res[res$variable == "S" & res$age == 2 & res$vax == 3 &
             res$particle == 4 & res$group == "mid", ]
  expect_equal(r$value, u$S[2, 3, 4, 2, ])
  r <- res[res$variable == "R" & res$age == 1 & res$particle == 2 &
             res$group == "high", ]
  expect_equal(r$value, u$R[1, 2, 3, ])
  expect_true(all(is.na(r$vax)))
  r <- res[res$variable == "incidence" & res$particle == 1 &
             res$group == "low", ]
  expect_equal(r$value, u$incidence[1, 1, ])
})

test_that("particle and group dimensions are not confused", {
  t <- 0:3
  sys <- dust2::dust_system_create(model("sir"), list(), n_particles = 3)
  dust2::dust_system_set_state_initial(sys)
  res <- tidy_simulation(dust2::dust_system_simulate(sys, t), sys, t)
  expect_named(res, c("particle", "time", "variable", "value"))

  sys <- dust2::dust_system_create(model("sir"), rep(list(list()), 3),
                                   n_groups = 3)
  dust2::dust_system_set_state_initial(sys)
  res <- tidy_simulation(dust2::dust_system_simulate(sys, t), sys, t)
  expect_named(res, c("group", "time", "variable", "value"))
})

test_that("tidy state, index_state and vars", {
  sys <- dust2::dust_system_create(model("sir"), list(), n_particles = 3)
  dust2::dust_system_set_state_initial(sys)
  s <- dust2::dust_system_state(sys)
  res <- tidy_simulation(s, sys, NULL)
  expect_named(res, c("particle", "variable", "value"))
  expect_equal(res$value, c(s))

  t <- 0:4
  idx <- unlist(dust2::dust_unpack_index(sys)[c("I", "incidence")])
  y <- dust2::dust_system_simulate(sys, t, index_state = idx)
  res <- tidy_simulation(y, sys, t, index_state = idx, vars = "incidence")
  expect_equal(unique(res$variable), "incidence")
  expect_equal(res$value, c(y[2, , ]))
  expect_error(tidy_simulation(y, sys, t), "Unexpected array dimensions")
  expect_error(tidy_simulation(y, NULL, t), "dust_system")
})

test_that("tidy filter trajectories", {
  filter <- dust2::dust_filter_create(model("sir"), time_start = 0,
                                      data = sir_data(), n_particles = 20,
                                      dt = 0.25, seed = 1)
  expect_error(tidy_filter(filter), "not yet been run")
  dust2::dust_likelihood_run(filter, list(), save_trajectories = TRUE)
  h <- dust2::dust_likelihood_last_trajectories(filter)
  res <- tidy_filter(filter)
  expect_named(res, c("particle", "time", "variable", "value"))
  r <- res[res$variable == "incidence" & res$particle == 7, ]
  expect_equal(r$value, h[4, 7, ])
  expect_equal(r$time, 1:10)

  res <- tidy_filter(filter, select_random_particle = TRUE)
  expect_named(res, c("time", "variable", "value"))
  expect_equal(nrow(res), 40)
})

test_that("tidy grouped filter and unfilter", {
  d <- rbind(cbind(sir_data(), group = "a"), cbind(sir_data(), group = "b"))
  d <- dust2::dust_filter_data(d, group = "group")
  filter <- dust2::dust_filter_create(model("sir"), time_start = 0,
                                      data = d, n_particles = 5, dt = 0.25)
  dust2::dust_likelihood_run(filter, list(list(), list()),
                             save_trajectories = TRUE)
  h <- dust2::dust_likelihood_last_trajectories(filter)
  res <- tidy_filter(filter)
  expect_named(res, c("group", "particle", "time", "variable", "value"))
  r <- res[res$variable == "I" & res$particle == 3 & res$group == "b", ]
  expect_equal(r$value, h[2, 3, 2, ])

  unfilter <- dust2::dust_unfilter_create(model("sir"), time_start = 0,
                                          data = sir_data(), dt = 0.25)
  dust2::dust_likelihood_run(unfilter, list(), save_trajectories = TRUE)
  res <- tidy_filter(unfilter)
  expect_named(res, c("time", "variable", "value"))
})
