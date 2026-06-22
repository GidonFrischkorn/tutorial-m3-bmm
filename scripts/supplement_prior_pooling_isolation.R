#' >>> RETIRED (2026-06-20): superseded by the recalibrated bmm M3 default priors
#' >>> (venpopov/bmm#366). This analytic prior-predictive simulation isolated the
#' >>> pooling effect behind the OLD default-prior mismatch (softmax ~89% vs simple
#' >>> ~55% correct). Under the recalibrated defaults both rules imply ~90%, so the
#' >>> mismatch this script dissected no longer arises. Supplement 1, Section 1 now
#' >>> keeps only the qualitative pooling mechanism. Kept for archival reference.
#'
#' Supplement 1, Section 1: Isolating the pooling effect from the prior magnitude
#'
#' The simple and softmax choice rules carry DIFFERENT bmm default priors, so the
#' default-prior predictive divergence (softmax ~89% vs simple ~55% correct) confounds
#' two things: (i) the different prior magnitudes and (ii) how each rule pools multiple
#' activation sources (sum-of-exponentials vs exponential-of-a-sum).
#'
#' This script isolates (ii) by placing GENUINELY IDENTICAL priors on a and c for both
#' rules and recomputing the implied correct-response proportion. No model fitting: this
#' is a direct prior-predictive simulation of the representative recognition condition
#' used in Tutorial 1 (n_other = 2, n_npl = 3).

###############################################################################!
# 0) R Setup ------------------------------------------------------------------
###############################################################################!

pacman::p_load(here, tidyverse)

# Large sample so Monte Carlo error is negligible at the reported precision.
n_draws     <- 2e6
n_other_rep <- 2
n_npl_rep   <- 3

#' Implied P(correct) for the representative condition under each choice rule.
#' Softmax uses the identity link with fixed baseline b = 0; the simple rule uses
#' the log link with fixed baseline b = 0.1 (the bmm defaults for each rule).
prior_pred_pcorr <- function(a, c) {
  # Softmax: P(corr) = exp(a + c) / [exp(a + c) + n_other * exp(a) + n_npl * exp(0)]
  denom_soft <- exp(a + c) + n_other_rep * exp(a) + n_npl_rep
  p_soft     <- exp(a + c) / denom_soft

  # Simple: activations act_k = b + exp(.), b = 0.1
  act_corr   <- 0.1 + exp(a) + exp(c)
  act_other  <- 0.1 + exp(a)
  act_npl    <- 0.1
  denom_simp <- act_corr + n_other_rep * act_other + n_npl_rep * act_npl
  p_simp     <- act_corr / denom_simp

  tibble(softmax = p_soft, simple = p_simp)
}

summarise_pcorr <- function(df, label) {
  df |>
    pivot_longer(everything(), names_to = "rule", values_to = "p_corr") |>
    group_by(rule) |>
    summarise(
      median = median(p_corr),
      q05    = quantile(p_corr, 0.05),
      q95    = quantile(p_corr, 0.95),
      .groups = "drop"
    ) |>
    mutate(scenario = label, .before = 1)
}

###############################################################################!
# 1) Identical modest priors: N(1, 0.5^2) on both a and c --------------------
###############################################################################!

# With identical, modest priors the two rules predict nearly the same correct-
# response proportion: pooling alone contributes little at these activation values.
draws_modest <- prior_pred_pcorr(
  a = rnorm(n_draws, mean = 1, sd = 0.5),
  c = rnorm(n_draws, mean = 1, sd = 0.5)
)
res_modest <- summarise_pcorr(draws_modest, "Identical N(1, 0.5^2)")

###############################################################################!
# 2) Identical large priors: softmax-scale N(2,1) / N(3,1) on both rules -----
###############################################################################!

# Giving BOTH rules the larger softmax-default intercept priors shows where pooling
# bites: at large activations the exponential-of-a-sum grows multiplicatively, so the
# softmax correct-response proportion outruns the simple rule even with identical priors.
draws_large <- prior_pred_pcorr(
  a = rnorm(n_draws, mean = 2, sd = 1),
  c = rnorm(n_draws, mean = 3, sd = 1)
)
res_large <- summarise_pcorr(draws_large, "Identical N(2,1)/N(3,1)")

###############################################################################!
# 3) Report -------------------------------------------------------------------
###############################################################################!

results <- bind_rows(res_modest, res_large) |>
  mutate(across(c(median, q05, q95), ~ round(100 * .x, 1)))

cat("Implied P(correct) for the representative condition (n_other = 2, n_npl = 3)\n")
cat("by choice rule, under IDENTICAL priors on a and c (percent):\n\n")
print(results)
