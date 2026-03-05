#' Tutorial 4: Parameter Recovery for the Custom Filtering M3
#'
#' This script demonstrates parameter recovery as validation for the custom
#' filtering M3 defined in Tutorial 3. We simulate data with known parameter
#' values, fit the model, and check whether the generating values are recovered.
#'
#' The custom filtering M3 replaces the single f parameter from the cs version
#' with two ratio parameters (ra, rc). Recovery focuses on all four parameters:
#' a (item memory), c (context binding), ra (item ratio), rc (binding ratio).
#'
#' In this script, you will see:
#'  1) How to choose generating parameter values on the link scale
#'  2) How to simulate data with rm3()
#'  3) How to evaluate group-level and individual-level recovery
#'  4) How recovery quality informs study design decisions

###############################################################################!
# 0) R Setup: Packages & Settings ---------------------------------------------
###############################################################################!

pacman::p_load(here, bmm, brms, tidyverse, tidybayes, patchwork)
source(here("functions", "clean_plot.R"))

# -- Fitting settings ----------------------------------------------------------
chains <- 4
cores  <- 4
warmup <- 1000
iter   <- 2000

# -- Reproducibility -----------------------------------------------------------
set.seed(2026)

###############################################################################!
# 1) Model Specification -------------------------------------------------------
###############################################################################!

# We use the same custom filtering M3 as in Tutorial 3.
# Five response categories, four free parameters (a, c, ra, rc).
#
# Activation formulas:
#   correct = b + a + c
#   distc   = b + ra*a + rc*c
#   other   = b + a
#   disto   = b + ra*a
#   npl     = b
#
# For simulation, the recovery uses a single condition (no experimental
# manipulation), so ra and rc have a single group-level mean each.

## 1.1) Define the custom M3 object --------------------------------------------

m3_model_custom <- m3(
  resp_cats   = c("corr", "distc", "other", "disto", "npl"),
  num_options = c(
    "n_corr" = 1, "n_distc" = 1, "n_other" = 4,
    "n_disto" = 4, "n_npl" = 5
  ),
  version     = "custom",
  choice_rule = "softmax",
  links = list(
    a  = "identity",
    c  = "identity",
    ra = "logit",
    rc = "logit"
  )
)

## 1.2) Define the activation formulas -----------------------------------------

# For simulation with rm3(), we only need the activation formulas:
m3_act_funs <- bmf(
  corr  ~ a + c + b,
  distc ~ ra * a + rc * c + b,
  other ~ a + b,
  disto ~ ra * a + b,
  npl   ~ b
)

## 1.3) Define the fitting formula ---------------------------------------------

# For model fitting, we also need prediction formulas: intercept + random
# intercept per participant for each parameter.
m3_formula_fit <- bmf(
  corr  ~ a + c + b,
  distc ~ ra * a + rc * c + b,
  other ~ a + b,
  disto ~ ra * a + b,
  npl   ~ b,
  a  ~ 1 + (1 | participant),
  c  ~ 1 + (1 | participant),
  ra ~ 1 + (1 | participant),
  rc ~ 1 + (1 | participant)
)

###############################################################################!
# 2) Data Simulation -----------------------------------------------------------
###############################################################################!

## 2.1) Choose generating parameter values -------------------------------------

# Group-level means and between-person SDs on the link scale.
# a and c: identity link (link = native)
# ra and rc: logit link (native = plogis(link))
#
# Values are chosen to approximate empirical estimates from Li et al. (2026)
# and to produce realistic individual differences.

# Group-level means (on link scale)
mean_a  <- 1.5      # identity link -> native = 1.5
mean_c  <- 2.0      # identity link -> native = 2.0
mean_ra <- 0.85     # logit link -> native = plogis(0.85) ~ 0.70
mean_rc <- 0.25     # logit link -> native = plogis(0.25) ~ 0.56

# Between-person SDs (on link scale)
sd_a  <- 0.50
sd_c  <- 0.50
sd_ra <- 0.40
sd_rc <- 0.40

# Quick check: what do these look like on the native scale?
cat("Parameter ranges (mean +/- 2 SD on native scale):\n")
cat("  a: ", mean_a - 2*sd_a, "to", mean_a + 2*sd_a, "(identity)\n")
cat("  c: ", mean_c - 2*sd_c, "to", mean_c + 2*sd_c, "(identity)\n")
cat("  ra:", round(plogis(mean_ra - 2*sd_ra), 2), "to",
    round(plogis(mean_ra + 2*sd_ra), 2), "(logit -> probability)\n")
cat("  rc:", round(plogis(mean_rc - 2*sd_rc), 2), "to",
    round(plogis(mean_rc + 2*sd_rc), 2), "(logit -> probability)\n")

## 2.2) Simulation settings ----------------------------------------------------

N               <- 50    # number of participants
trials_per_cond <- 60    # trials in the single condition

## 2.3) Draw true individual parameters ----------------------------------------

true_pars <- tibble(
  participant = 1:N,
  a_link  = rnorm(N, mean_a, sd_a),
  c_link  = rnorm(N, mean_c, sd_c),
  ra_link = rnorm(N, mean_ra, sd_ra),
  rc_link = rnorm(N, mean_rc, sd_rc)
)

# Convert to native scale for use in rm3()
true_pars <- true_pars %>%
  mutate(
    a_native  = a_link,          # identity
    c_native  = c_link,          # identity
    ra_native = plogis(ra_link), # logit -> probability
    rc_native = plogis(rc_link)  # logit -> probability
  )

cat("\nTrue parameter summary (native scale):\n")
true_pars %>%
  summarise(
    across(ends_with("_native"),
           list(mean = mean, sd = sd, min = min, max = max),
           .names = "{.col}_{.fn}")
  ) %>%
  pivot_longer(everything(),
               names_to = c("parameter", "stat"),
               names_pattern = "(.+)_native_(.+)") %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  print()

## 2.4) Simulate the full dataset ----------------------------------------------

# Single condition: all participants receive the same number of trials.
# Option counts match the Li et al. design:
#   1 correct, 1 paired distractor, 4 other, 4 other distractor, 5 NPL
data_list <- list()

for (i in seq_len(N)) {

  pars_i <- c(
    a  = true_pars$a_native[i],
    c  = true_pars$c_native[i],
    ra = true_pars$ra_native[i],
    rc = true_pars$rc_native[i]
  )

  sim_i <- rm3(
    n        = 1,
    size     = trials_per_cond,
    pars     = pars_i,
    m3_model = m3_model_custom,
    act_funs = m3_act_funs
  )

  row_i <- data.frame(sim_i)
  row_i$participant <- i
  row_i$n_corr  <- 1L
  row_i$n_distc <- 1L
  row_i$n_other <- 4L
  row_i$n_disto <- 4L
  row_i$n_npl   <- 5L

  data_list[[i]] <- row_i

  if (i %% 10 == 0) cat("  Simulated participant", i, "of", N, "\n")
}

sim_data <- bind_rows(data_list)

cat("\nSimulated dataset:\n")
cat("  Rows:", nrow(sim_data), "(expected:", N, ")\n")
str(sim_data)

## 2.5) Inspect the simulated data ---------------------------------------------

cat("\nMean response proportions:\n")
sim_data %>%
  mutate(total = corr + other + distc + disto + npl) %>%
  summarise(
    corr  = mean(corr / total),
    distc = mean(distc / total),
    other = mean(other / total),
    disto = mean(disto / total),
    npl   = mean(npl / total)
  ) %>%
  print()

###############################################################################!
# 3) Model Fitting -------------------------------------------------------------
###############################################################################!

## 3.1) Define the model object for fitting ------------------------------------

m3_model_fit <- m3(
  resp_cats   = c("corr", "distc", "other", "disto", "npl"),
  num_options = c("n_corr", "n_distc", "n_other", "n_disto", "n_npl"),
  version     = "custom",
  choice_rule = "softmax",
  links = list(
    a  = "identity",
    c  = "identity",
    ra = "logit",
    rc = "logit"
  )
)

## 3.2) Fit the model ----------------------------------------------------------

fit <- bmm(
  formula = m3_formula_fit,
  model   = m3_model_fit,
  data    = sim_data,
  chains  = chains,
  cores   = cores,
  warmup  = warmup,
  iter    = iter,
  backend = "cmdstanr",
  file    = here("output", "fit_m3_recovery_filtering")
)

## 3.3) Convergence checks -----------------------------------------------------
cat("\nModel summary:\n")
summary(fit)

###############################################################################!
# 4) Recovery Evaluation -------------------------------------------------------
###############################################################################!

## 4.1) Group-level recovery ---------------------------------------------------

fe <- fixef(fit)

par_names  <- c("a", "c", "ra", "rc")
true_means <- c(mean_a, mean_c, mean_ra, mean_rc)

group_recovery <- tibble(
  parameter       = par_names,
  true_mean       = true_means,
  recovered_mean  = fe[paste0(par_names, "_Intercept"), "Estimate"],
  recovered_lower = fe[paste0(par_names, "_Intercept"), "Q2.5"],
  recovered_upper = fe[paste0(par_names, "_Intercept"), "Q97.5"]
) %>%
  mutate(
    bias     = recovered_mean - true_mean,
    coverage = true_mean >= recovered_lower & true_mean <= recovered_upper
  )

cat("\nGroup-level recovery:\n")
print(group_recovery)

## 4.2) Group-level recovery plot ----------------------------------------------

group_plot <- ggplot(group_recovery,
                     aes(x = true_mean, y = recovered_mean)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  geom_pointrange(aes(ymin = recovered_lower, ymax = recovered_upper),
                  size = 0.5) +
  geom_text(aes(label = parameter), nudge_x = 0.15, size = 3.5) +
  labs(x = "True generating value (link scale)",
       y = "Recovered estimate (link scale)") +
  clean_plot()

group_plot

ggsave(here("figures", "tutorial4_group_recovery.pdf"),
       group_plot, width = 5, height = 5)

## 4.3) SD recovery ------------------------------------------------------------

vc <- VarCorr(fit)

true_sds <- c(sd_a, sd_c, sd_ra, sd_rc)

sd_recovery <- tibble(
  parameter    = par_names,
  true_sd      = true_sds,
  recovered_sd = NA_real_,
  sd_lower     = NA_real_,
  sd_upper     = NA_real_
)

for (k in seq_along(par_names)) {
  p <- par_names[k]
  sd_est <- vc$participant$sd[paste0(p, "_Intercept"), ]
  sd_recovery$recovered_sd[k] <- sd_est["Estimate"]
  sd_recovery$sd_lower[k]     <- sd_est["Q2.5"]
  sd_recovery$sd_upper[k]     <- sd_est["Q97.5"]
}

sd_recovery <- sd_recovery %>%
  mutate(
    sd_bias     = recovered_sd - true_sd,
    sd_coverage = true_sd >= sd_lower & true_sd <= sd_upper
  )

cat("\nSD recovery (between-person variability):\n")
print(sd_recovery)

## 4.4) Individual-level recovery ----------------------------------------------

co <- coef(fit)

indiv_recovery <- tibble()

for (p in par_names) {
  est <- co$participant[, , paste0(p, "_Intercept")]

  indiv_p <- tibble(
    participant    = 1:N,
    parameter      = p,
    recovered_link = est[1:N, "Estimate"],
    indiv_lower    = est[1:N, "Q2.5"],
    indiv_upper    = est[1:N, "Q97.5"]
  )

  indiv_recovery <- bind_rows(indiv_recovery, indiv_p)
}

# Add true link-scale values
true_long <- true_pars %>%
  select(participant,
         a = a_link, c = c_link, ra = ra_link, rc = rc_link) %>%
  pivot_longer(-participant, names_to = "parameter",
               values_to = "true_link")

indiv_recovery <- indiv_recovery %>%
  left_join(true_long, by = c("participant", "parameter"))

# Individual-level metrics
cat("\nIndividual-level recovery metrics:\n")
indiv_summary <- indiv_recovery %>%
  group_by(parameter) %>%
  summarise(
    correlation = cor(true_link, recovered_link),
    mean_bias   = mean(recovered_link - true_link),
    rmse        = sqrt(mean((recovered_link - true_link)^2)),
    coverage    = mean(true_link >= indiv_lower & true_link <= indiv_upper),
    .groups     = "drop"
  )

print(indiv_summary)

## 4.5) Individual-level scatter plots -----------------------------------------

corr_labels <- indiv_summary %>%
  mutate(label = paste0("r = ", sprintf("%.2f", correlation)))

scatter_plot <- ggplot(indiv_recovery,
                       aes(x = true_link, y = recovered_link)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_text(data = corr_labels,
            aes(x = Inf, y = -Inf, label = label),
            hjust = 1.1, vjust = -0.5, size = 3.5, inherit.aes = FALSE) +
  facet_wrap(~ parameter, scales = "free") +
  labs(x = "True value (link scale)",
       y = "Recovered value (link scale)") +
  clean_plot()

scatter_plot

ggsave(here("figures", "tutorial4_indiv_scatter.pdf"),
       scatter_plot, width = 6.5, height = 5)
