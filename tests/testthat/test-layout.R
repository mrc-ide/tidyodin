## Line coverage shows that the reshaping code runs, not that every
## value ends up on the right row. These tests check that directly.

## Fill an array with its own linear index, so that each value says
## where it came from, then check that every row of the tidy output
## carries the columns that match its value.
check_melt <- function(shape, trailing, vars = NULL) {
  packer <- monty::monty_packer(array = shape)
  names <- packer$names()
  size <- c(length(names), lengths(trailing, use.names = FALSE))
  x <- array(seq_len(prod(size)), size)
  n_samples <- if (is.null(trailing$sample)) NULL else 3
  res <- melt_state(x, names, trailing, vars = vars, n_samples = n_samples)

  ## Position on the first dimension, from the variable and its indices
  index <- packer$index()
  state <- vapply(seq_len(nrow(res)), function(i) {
    idx <- index[[as.character(res$variable[[i]])]]
    cols <- grep("^index", names(res), value = TRUE)
    at <- unlist(res[i, cols], use.names = FALSE)
    at <- at[!is.na(at)]
    if (length(at) == 0) idx[[1]] else idx[rbind(at)]
  }, numeric(1))

  ## Position on every other dimension, from the column of that name
  expected <- state
  stride <- size[[1]]
  for (i in seq_along(trailing)) {
    at <- match(res[[names(trailing)[[i]]]], trailing[[i]])
    expected <- expected + (at - 1) * stride
    stride <- stride * size[[i + 1]]
  }
  expect_equal(res$value, expected)

  if (!is.null(trailing$sample)) {
    expect_equal(res$sample, (res$chain - 1L) * n_samples + res$iteration)
    expect_true(all(res$iteration >= 1 & res$iteration <= n_samples))
  }
  kept <- vars %||% names(shape)
  expect_setequal(as.character(res$variable), kept)
  expect_equal(nrow(res), sum(parse_names(names)$variable %in% kept) *
                 prod(size[-1]))
  invisible(res)
}


test_that("every value lands on the row that describes it", {
  shapes <- list(
    scalars = list(a = NULL, b = NULL),
    vector = list(a = NULL, v = 3),
    matrix = list(m = c(2, 3), a = NULL, v = 2),
    cube = list(a = NULL, k = c(2, 2, 3), v = 4))
  trailings <- list(
    time = list(time = c(0, 0.5, 1, 2.5)),
    state_only = list(particle = 1:3),
    particles = list(particle = 1:3, time = 0:2),
    groups = list(group = c("x", "y"), time = c(1, 2, 4)),
    both = list(particle = 1:2, group = 1:3, time = 0:3),
    samples = list(time = 1:4, sample = 1:6),
    samples_subset = list(time = 1:2, sample = c(2L, 3L, 5L)),
    everything = list(particle = 1:2, group = c("x", "y"), time = 1:3,
                      sample = 1:6))
  for (shape in shapes) {
    for (trailing in trailings) {
      check_melt(shape, trailing)
    }
  }
  check_melt(shapes$cube, trailings$everything, vars = "k")
  check_melt(shapes$matrix, trailings$both, vars = c("v", "m"))
  check_melt(shapes$vector, trailings$samples, vars = "a")
})


test_that("tidy output matches dust2 for every layout of a system", {
  t <- c(0, 1, 2, 5)
  for (name in c("sir", "age")) {
    for (n_particles in c(1, 3)) {
      for (n_groups in c(1, 2)) {
        for (deterministic in c(FALSE, TRUE)) {
          pars <- if (n_groups == 1) list() else
            list(list(beta = 0.2), list(beta = 0.6))
          sys <- dust2::dust_system_create(
            model(name), pars, n_particles = n_particles,
            n_groups = n_groups, deterministic = deterministic, seed = 1)
          dust2::dust_system_set_state_initial(sys)

          state <- dust2::dust_system_state(sys)
          res <- tidy_simulation(state, sys, NULL)
          expect_equal(nrow(res), length(state))
          unpacked <- dust2::dust_unpack_state(sys, state)
          for (v in names(unpacked)) {
            expect_equal(res$value[res$variable == v], c(unpacked[[v]]))
          }

          y <- dust2::dust_system_simulate(sys, t)
          res <- tidy_simulation(y, sys, t)
          expect_equal(nrow(res), length(y))
          expect_equal("particle" %in% names(res), n_particles > 1)
          expect_equal("group" %in% names(res), n_groups > 1)
          expect_setequal(res$time, t)
          ## dust2 unpacks each variable with its own dimensions first
          ## and the system's after, which is the order of the rows here
          unpacked <- dust2::dust_unpack_state(sys, y)
          expect_equal(levels(res$variable), names(unpacked))
          for (v in names(unpacked)) {
            expect_equal(res$value[res$variable == v], c(unpacked[[v]]))
          }
        }
      }
    }
  }
})


test_that("runs over parameter sets match dust2 for every layout", {
  t <- c(0, 2, 4)
  pars <- data.frame(beta = c(0.2, 0.4, 0.6))
  for (name in c("sir", "age")) {
    for (n_particles in c(1, 3)) {
      for (deterministic in c(FALSE, TRUE)) {
        res <- run_parameter_sets(model(name), pars, t,
                                  n_particles = n_particles,
                                  deterministic = deterministic, seed = 1)
        sys <- dust2::dust_system_create(
          model(name), lapply(pars$beta, function(b) list(beta = b)),
          n_particles = n_particles, n_groups = 3,
          deterministic = deterministic, seed = 1,
          preserve_particle_dimension = TRUE)
        dust2::dust_system_set_state_initial(sys)
        ## state x particle x group x time
        y <- dust2::dust_system_simulate(sys, t)
        unpacked <- dust2::dust_unpack_state(sys, y)

        tidy <- tidy_trajectories(res)
        expect_equal(nrow(tidy), length(y))
        expect_equal("particle" %in% names(tidy), n_particles > 1)
        for (set in 1:3) {
          for (particle in seq_len(n_particles)) {
            keep <- tidy$sample == set & tidy$variable == "incidence"
            if (n_particles > 1) {
              keep <- keep & tidy$particle == particle
            }
            expect_equal(tidy$value[keep], unpacked$incidence[particle, set, ])
          }
        }
        end <- tidy_state(res)
        expect_equal(nrow(end), length(y) / length(t))
        expect_equal(
          end$value[end$variable == "incidence" & end$sample == 2],
          unpacked$incidence[, 2, length(t)])
      }
    }
  }
})
