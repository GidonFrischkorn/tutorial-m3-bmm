# Supplement 1, Section 2: predictive comparison of the softmax and simple
# choice rules with PSIS-LOO and WAIC.
# Input: the long-run choice-rule fits saved by tutorial1_simple_span.R.
# Output: output/loo_choice_rule.rds

###############################################################################!
# 0) R Setup -------------------------------------------------------------------
###############################################################################!
pacman::p_load(here, bmm, brms, loo)

###############################################################################!
# 1) Load the long-run choice-rule fits ----------------------------------------
###############################################################################!
fit_softmax <- readRDS(here("output", "fit_m3_ss_softmax_longrun.rds"))
fit_simple  <- readRDS(here("output", "fit_m3_ss_simple_longrun.rds"))

###############################################################################!
# 2) Predictive comparison -----------------------------------------------------
###############################################################################!

## 2.1) LOO (PSIS-LOO) ---------------------------------------------------------
loo_softmax <- loo(fit_softmax)
loo_simple  <- loo(fit_simple)

print(loo_softmax)
print(loo_simple)

# PSIS-LOO is reliable only if all Pareto k < 0.7 (TRUE flags a problem); if
# so, use loo(fit, moment_match = TRUE) or reloo()
print(any(pareto_k_values(loo_softmax) > 0.7))
print(any(pareto_k_values(loo_simple) > 0.7))

## 2.2) loo_compare ------------------------------------------------------------
# An elpd difference within about two standard errors means both rules predict
# new data about equally well
loo_cmp <- loo_compare(loo_softmax, loo_simple)
print(loo_cmp)

## 2.3) WAIC -------------------------------------------------------------------
# WAIC as a check on the LOO estimates
waic_softmax <- waic(fit_softmax)
waic_simple  <- waic(fit_simple)

print(waic_softmax)
print(waic_simple)

###############################################################################!
# 3) Save results --------------------------------------------------------------
###############################################################################!
# LOO, WAIC, and loo_compare() results for both choice rules
loo_results <- list(
  loo_softmax  = loo_softmax,
  loo_simple   = loo_simple,
  loo_compare  = loo_cmp,
  waic_softmax = waic_softmax,
  waic_simple  = waic_simple
)
saveRDS(loo_results, here("output", "loo_choice_rule.rds"))
