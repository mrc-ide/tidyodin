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

sys <- dust2::dust_system_create(sir, list(), n_particles = 30, dt = 0.25)
dust2::dust_system_set_state_initial(sys)
t <- 0:60
y <- dust2::dust_system_simulate(sys, t)

ggplot(tidy_simulation(y, sys, t), aes(time, value, group = particle)) +
  geom_line(alpha = 0.3) +
  facet_wrap(~variable, scales = "free_y")
```

See the [Get started](https://mrc-ide.github.io/tidyodin/articles/tidyodin.html) article, or `vignette("tidyodin")`, for fits, projections and counterfactuals.
