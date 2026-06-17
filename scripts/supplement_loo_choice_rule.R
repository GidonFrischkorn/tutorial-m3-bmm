#' Supplement: Predictive comparison (LOO/WAIC) of the choice-rule models
#'
#' Reviewer 1 questioned whether the Bayes factor comparing the non-nested
#' softmax and simple choice rules is meaningful, given that "fair" priors for
#' non-nested models are hard to specify. The Bayes factor answers a
#' prior-predictive question (which model's priors anticipated the data, with an
#' Occam penalty) and is therefore prior-sensitive. Leave-one-out
#' cross-validation (LOO) and WAIC answer a different, posterior-predictive
#' question (which model predicts new observations), conditioning on the
#' posterior rather than integrating over the prior -- so they are far less
#' sensitive to the non-nested prior that is hard to calibrate.
#'
#' This script computes LOO and WAIC for the two long-run choice-rule fits and
#' compares them with loo_compare(). Both models are fit to the SAME
#' observations, so the comparison is valid without any prior matching -- which
#' is precisely the point.
#'
#' KEY RESULT (see Supplement 1): the two rules are predictively
#' indistinguishable (elpd difference ~ -0.7, SE ~ 4.1), even though the Bayes
#' factor overwhelmingly favors softmax. BF and LOO genuinely diverge here
#' because they ask different questions -- the cleanest illustration of why both
#' should be reported for non-nested comparisons.
#'
#' No refits are required: brms derives the pointwise log-likelihood from the
#' saved 80,000-sample fit objects on demand.

###############################################################################!
# 0) R Setup ------------------------------------------------------------------
###############################################################################!
pacman::p_load(here, bmm, brms, loo, tidyverse)

###############################################################################!
# 1) Load the long-run choice-rule fits ---------------------------------------
###############################################################################!
# 80,000 post-warmup samples each; produced by tutorial1_simple_span.R and the
# bridge-sampling supplement scripts.
fit_softmax <- readRDS(here("output", "fit_m3_ss_softmax_longrun.rds"))
fit_simple  <- readRDS(here("output", "fit_m3_ss_simple_longrun.rds"))

###############################################################################!
# 2) Predictive comparison ----------------------------------------------------
###############################################################################!

## 2.1) LOO (PSIS-LOO) ---------------------------------------------------------
loo_softmax <- loo(fit_softmax)
loo_simple  <- loo(fit_simple)

cat("\n=== LOO: softmax ===\n");  print(loo_softmax)
cat("\n=== LOO: simple ===\n");   print(loo_simple)

# Pareto-k diagnostics: PSIS-LOO is reliable only when all k < 0.7. For the
# Tutorial 1 multinomial fits, all k are good, so no moment-matching or reloo()
# fallback is needed. If a future configuration shows k > 0.7, recompute with
# loo(fit, moment_match = TRUE) or reloo(), or fall back to WAIC with a caveat.
cat("\nPareto-k (softmax) -- any > 0.7? ",
    any(pareto_k_values(loo_softmax) > 0.7), "\n")
cat("Pareto-k (simple)  -- any > 0.7? ",
    any(pareto_k_values(loo_simple) > 0.7), "\n")

## 2.2) loo_compare ------------------------------------------------------------
loo_cmp <- loo_compare(loo_softmax, loo_simple)
cat("\n=== loo_compare (softmax vs simple) ===\n")
print(loo_cmp)
cat("\nInterpretation: an elpd difference smaller than ~2x its SE indicates the",
    "\ntwo specifications predict new data about equally well.\n")

## 2.3) WAIC (corroboration) ---------------------------------------------------
waic_softmax <- waic(fit_softmax)
waic_simple  <- waic(fit_simple)
cat("\n=== WAIC: softmax ===\n"); print(waic_softmax)
cat("\n=== WAIC: simple ===\n");  print(waic_simple)

###############################################################################!
# 3) Save results -------------------------------------------------------------
###############################################################################!
loo_results <- list(
  loo_softmax  = loo_softmax,
  loo_simple   = loo_simple,
  loo_compare  = loo_cmp,
  waic_softmax = waic_softmax,
  waic_simple  = waic_simple
)
saveRDS(loo_results, here("output", "loo_choice_rule.rds"))

cat("\nSaved: output/loo_choice_rule.rds\n")
