# Supplement 1: Monte Carlo stability of the Savage-Dickey Bayes factors
#
# Recomputes the hypothesis() tests reported in Tutorials 1-3 on the default
# fits (4,000 post-warmup draws) and on the long-run fits (80,000 draws).
# hypothesis() estimates both the prior and the posterior density from draws,
# so each test is repeated over 10 seeds to show its Monte Carlo range.
# Requires the default and long-run fits from Tutorials 1-3 in output/.
# Output: output/sd_stability.rds

###############################################################################!
# 0) R Setup ------------------------------------------------------------------
###############################################################################!
pacman::p_load(here, bmm, brms, tidyverse)

seeds <- c(2026, 1:9)  # 2026 is the seed used in the manuscript

###############################################################################!
# 1) Tests reported in the tutorial -------------------------------------------
###############################################################################!
tests <- tribble(
  ~tutorial, ~fit_file,                ~label,                       ~test,
  1, "fit_m3_ss_softmax",       "Set size effect on c",       "(c_expopenset:ss_lin + c_expclosedset:ss_lin) / 2 = 0",
  1, "fit_m3_ss_softmax",       "Set size effect on a",       "(a_expopenset:ss_lin + a_expclosedset:ss_lin) / 2 = 0",
  1, "fit_m3_ss_softmax",       "Set size effect on c by experiment", "c_expopenset:ss_lin = c_expclosedset:ss_lin",
  2, "fit_m3_cs",               "f: pre vs. retro",           "f_conditionpre - f_conditionretro = 0",
  2, "fit_m3_cs",               "a: control vs. distractor",  "a_conditioncontrol - (a_conditionpre + a_conditionretro) / 2 = 0",
  2, "fit_m3_cs",               "c: control vs. distractor",  "c_conditioncontrol - (c_conditionpre + c_conditionretro) / 2 = 0",
  3, "fit_m3_custom_filtering", "ra vs. rc (pre)",            "ra_conditionpre - rc_conditionpre = 0",
  3, "fit_m3_custom_filtering", "ra vs. rc (retro)",          "ra_conditionretro - rc_conditionretro = 0",
  3, "fit_m3_custom_filtering", "ra: pre vs. retro",          "ra_conditionpre - ra_conditionretro = 0",
  3, "fit_m3_custom_filtering", "rc: pre vs. retro",          "rc_conditionpre - rc_conditionretro = 0"
)

###############################################################################!
# 2) Evidence ratios across seeds and sample sizes ----------------------------
###############################################################################!

# Evid.Ratio from hypothesis() is BF01 (posterior / prior density at 0)
bf01_by_seed <- function(fit, test) {
  map_dbl(seeds, \(s) hypothesis(fit, test, seed = s)$hypothesis$Evid.Ratio)
}

run_fit <- function(fit_name, run) {
  file <- if (run == "default") fit_name else paste0(fit_name, "_longrun")
  fit  <- readRDS(here("output", paste0(file, ".rds")))
  out  <- tests |>
    filter(fit_file == fit_name) |>
    mutate(run = run, draws = ndraws(fit))
  out$bf01 <- map(out$test, \(h) bf01_by_seed(fit, h))
  rm(fit)
  gc()
  out
}

sd_draws <- expand_grid(fit_name = unique(tests$fit_file),
                        run = c("default", "longrun")) |>
  pmap(\(fit_name, run) run_fit(fit_name, run)) |>
  list_rbind() |>
  unnest(bf01)

###############################################################################!
# 3) Summary ------------------------------------------------------------------
###############################################################################!
sd_stability <- sd_draws |>
  summarise(
    draws       = first(draws),
    median_bf01 = median(bf01),
    min_bf01    = min(bf01),
    max_bf01    = max(bf01),
    .by = c(tutorial, label, test, run)
  ) |>
  mutate(ratio_max_min = max_bf01 / min_bf01)

print(sd_stability, n = Inf, width = Inf)

saveRDS(sd_stability, here("output", "sd_stability.rds"))
