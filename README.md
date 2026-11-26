# tidyodin

Tidy output from [odin2](https://mrc-ide.github.io/odin2), [dust2](https://mrc-ide.github.io/dust2) and [monty](https://mrc-ide.github.io/monty), ready for ggplot2.

dust2 and monty return arrays whose layout depends on how the model was run. This package converts them to long-format tibbles with one row per value.

| Function | Converts |
|---|---|
| `tidy_simulation()` | output of `dust_system_simulate()` or `dust_system_state()` |
| `tidy_filter()` | trajectories from a filter or unfilter |
| `tidy_pars()`, `tidy_density()` | parameters and log density from `monty_sample()` |
| `tidy_trajectories()`, `tidy_state()`, `tidy_snapshots()` | trajectories, final state and snapshots saved with samples |
| `summarise_trajectories()` | any of the above to a median and interval per time |

Two functions run a model over many parameter sets and return an object with the same structure as monty samples, so the functions above work on the result unchanged:

| Function | Does |
|---|---|
| `run_parameter_sets()` | runs from the initial state, for posterior draws, a data frame of scenarios or a list, with optional `fixed` parameters, several particles per set, or a deterministic run |
| `continue_parameter_sets()` | carries on from the final state and time of a fit or a previous run |

## Installation

```r
# install.packages("remotes")
remotes::install_github("mrc-ide/tidyodin")
```

odin2, dust2 and monty are installed from <https://mrc-ide.r-universe.dev>.

## Example

```r
library(tidyodin)
library(ggplot2)

# A stochastic SIR model
sir <- odin2::odin({
  update(S) <- S - n_SI
  update(I) <- I + n_SI - n_IR
  update(R) <- R + n_IR
  p_SI <- 1 - exp(-beta * I / N * dt)
  p_IR <- 1 - exp(-gamma * dt)
  n_SI <- Binomial(S, p_SI)
  n_IR <- Binomial(I, p_IR)
  initial(S) <- N - I0
  initial(I) <- I0
  initial(R) <- 0
  N <- parameter(1000)
  I0 <- parameter(10)
  beta <- parameter(0.4)
  gamma <- parameter(0.1)
})

# Run 30 particles and tidy the output
sys <- dust2::dust_system_create(sir, list(), n_particles = 30, dt = 0.25)
dust2::dust_system_set_state_initial(sys)
t <- 0:60
y <- dust2::dust_system_simulate(sys, t)

ggplot(tidy_simulation(y, sys, t), aes(time, value, group = particle)) +
  geom_line(alpha = 0.3) +
  facet_wrap(~variable, scales = "free_y")
```

See the [Get started](https://mrc-ide.github.io/tidyodin/articles/tidyodin.html) article, or `vignette("tidyodin")`, for fits, projections and counterfactuals.
