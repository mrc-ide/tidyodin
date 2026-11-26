cache <- new.env()

model <- function(name) {
  testthat::skip_if_not_installed("odin2")
  if (is.null(cache[[name]])) {
    cache[[name]] <- switch(
      name,
      sir = odin2::odin({
        update(S) <- S - n_SI
        update(I) <- I + n_SI - n_IR
        update(R) <- R + n_IR
        update(incidence) <- incidence + n_SI
        p_SI <- 1 - exp(-beta * I / N * dt)
        p_IR <- 1 - exp(-gamma * dt)
        n_SI <- Binomial(S, p_SI)
        n_IR <- Binomial(I, p_IR)
        initial(S) <- N - I0
        initial(I) <- I0
        initial(R) <- 0
        initial(incidence, zero_every = 1) <- 0
        N <- parameter(1000)
        I0 <- parameter(10)
        beta <- parameter(0.2)
        gamma <- parameter(0.1)
        cases <- data()
        cases ~ Poisson(incidence)
      }, quiet = TRUE),
      age = odin2::odin({
        update(S[, ]) <- S[i, j] - n_SI[i, j]
        update(I[, ]) <- I[i, j] + n_SI[i, j] - n_IR[i, j]
        update(R[]) <- R[i] + sum(n_IR[i, ])
        update(incidence) <- incidence + sum(n_SI)
        p_SI <- 1 - exp(-beta * sum(I) / N * dt)
        p_IR <- 1 - exp(-gamma * dt)
        n_SI[, ] <- Binomial(S[i, j], p_SI)
        n_IR[, ] <- Binomial(I[i, j], p_IR)
        initial(S[, ]) <- N / (n_age * n_vax) - 1
        initial(I[, ]) <- 1
        initial(R[]) <- 0
        initial(incidence, zero_every = 1) <- 0
        N <- parameter(1000)
        beta <- parameter(0.2)
        gamma <- parameter(0.1)
        n_age <- parameter(2)
        n_vax <- parameter(3)
        dim(S, I, n_SI, n_IR) <- c(n_age, n_vax)
        dim(R) <- n_age
      }, quiet = TRUE),
      ode = odin2::odin({
        deriv(S) <- -beta * S * I / N
        deriv(I) <- beta * S * I / N - gamma * I
        initial(S) <- N - 1
        initial(I) <- 1
        N <- parameter(1000)
        beta <- parameter(0.4)
        gamma <- parameter(0.1)
      }, quiet = TRUE))
  }
  cache[[name]]
}

sir_data <- function() {
  data.frame(time = 1:10, cases = c(3, 5, 8, 12, 14, 21, 25, 30, 28, 33))
}

sir_packer <- function() {
  monty::monty_packer(c("beta", "gamma"), fixed = list(N = 1000, I0 = 10))
}

## A small fit, cached as it is used by several tests
sir_fit <- function(...) {
  key <- paste("fit", deparse(list(...)), collapse = "")
  if (is.null(cache[[key]])) {
    filter <- dust2::dust_filter_create(model("sir"), time_start = 0,
                                        data = sir_data(), n_particles = 20,
                                        dt = 0.25, seed = 1)
    likelihood <- dust2::dust_likelihood_monty(filter, sir_packer(), ...)
    prior <- monty::monty_dsl({
      beta ~ Exponential(mean = 0.3)
      gamma ~ Exponential(mean = 0.1)
    })
    set.seed(1)
    samples <- suppressMessages(monty::monty_sample(
      likelihood + prior,
      monty::monty_sampler_random_walk(diag(2) * 0.0004),
      30, initial = c(0.2, 0.1), n_chains = 3))
    cache[[key]] <- list(filter = filter, samples = samples)
  }
  cache[[key]]
}

## A small fit of the SIR model to two groups, with gamma shared and
## beta varying by group
sir_fit_grouped <- function(...) {
  key <- paste("fit_grouped", deparse(list(...)), collapse = "")
  if (is.null(cache[[key]])) {
    d <- rbind(cbind(sir_data(), group = "north"),
               cbind(transform(sir_data(), cases = cases * 2),
                     group = "south"))
    d <- dust2::dust_filter_data(d, group = "group")
    filter <- dust2::dust_filter_create(model("sir"), time_start = 0,
                                        data = d, n_particles = 20,
                                        dt = 0.25, seed = 1)
    packer <- monty::monty_packer_grouped(
      c("north", "south"), c("beta", "gamma"), shared = "gamma",
      fixed = list(N = 1000, I0 = 10))
    likelihood <- dust2::dust_likelihood_monty(filter, packer, ...)
    set.seed(1)
    samples <- suppressMessages(monty::monty_sample(
      likelihood,
      monty::monty_sampler_random_walk(diag(3) * 1e-4),
      20, initial = c(0.1, 0.2, 0.3), n_chains = 2))
    cache[[key]] <- list(filter = filter, samples = samples, packer = packer)
  }
  cache[[key]]
}
