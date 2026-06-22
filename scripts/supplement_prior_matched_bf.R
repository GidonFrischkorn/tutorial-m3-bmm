#' >>> RETIRED (2026-06-20): superseded by the recalibrated bmm M3 default priors
#' >>> (venpopov/bmm#366). This script documented the OLD default-prior mismatch
#' >>> (softmax ~89% vs simple ~55% prior-predictive correct). Under the
#' >>> recalibrated defaults both rules imply ~90% and the choice-rule Bayes factor
#' >>> (~0.37) agrees with LOO that the rules are indistinguishable, so the
#' >>> prior-sensitivity / matched-prior analysis is no longer reported in
#' >>> Supplement 1. Kept for archival reference only; not sourced by any document.
#'
#' Supplement: Prior-matched Bayes factor for the choice-rule comparison
#'
#' Reviewer 1 asked whether softmax is favored over the simple rule because it
#' fits the data better, or because the two rules imply different prior
#' predictive distributions under the same nominal priors (softmax implies ~89%
#' correct, simple ~55%; see Tutorial 1 and Supplement Section 1).
#'
#' IMPORTANT: the imbalance is NOT driven by the fixed baseline b. In softmax,
#' b is added to every category's activation and therefore cancels (softmax is
#' shift-invariant): changing b leaves the predictions unchanged. We confirmed
#' this empirically -- fixing the softmax baseline at log(0.1) left the prior
#' predictive at ~89%. The imbalance instead arises from the a/c intercept
#' priors interacting with the link-function difference (softmax forms
#' exp(a + c); the simple rule forms a sum of exponentials plus an additive
#' baseline). The correct fix is therefore to retune the a/c intercept priors.
#'
#' We calibrated the softmax intercept priors so that the implied prior
#' predictive P(correct) matches the simple rule's (median ~.55, 90% CI ~.42-.73)
#' for the representative condition (2 other-list options, 3 lures):
#'   a ~ normal(2.10, 0.64)   (default: normal(2, 1))
#'   c ~ normal(1.17, 0.38)   (default: normal(3, 1))
#' These reproduce the simple rule's P(correct) quantiles (.41/.57/.72). Note
#' the full category distribution cannot be matched exactly: softmax gives the
#' lures a fixed weight exp(0)=1 while the simple rule gives them the additive
#' baseline 0.1, so P(NPL) differs structurally regardless of the a/c priors.
#'
#' We then recompute the Bayes factor (prior-matched softmax vs. simple). If
#' softmax is still favored, the BF advantage reflects likelihood fit rather
#' than the prior-predictive imbalance on the correct-response proportion.
#'
#' NOTE: refits with 80,000 post-warmup samples for stable bridge sampling.
#' Fits are cached via `file`; rerunning is cheap once cached.

###############################################################################!
# 0) R Setup ------------------------------------------------------------------
###############################################################################!
pacman::p_load(here, bmm, brms, tidyverse, tidybayes, bridgesampling)

chains <- 4
cores  <- 4
warmup <- 1000
iter   <- 21000  # 80,000 post-warmup samples total (bridge-stable)

data_agg <- read_csv(
  here("data", "Oberauer_2019_SimpleSpan_agg.csv"),
  show_col_types = FALSE
) |>
  mutate(exp = factor(exp), id = factor(id))
if (!"ss_lin" %in% names(data_agg)) {
  data_agg <- mutate(data_agg, ss_lin = setsize - 5)
}

###############################################################################!
# 1) Model Specification ------------------------------------------------------
###############################################################################!

m3_formula <- bmf(
  c ~ 0 + exp + exp:ss_lin + (1 + ss_lin | gr(id, by = exp)),
  a ~ 0 + exp + exp:ss_lin + (1 + ss_lin | gr(id, by = exp))
)

m3_softmax <- m3(
  resp_cats   = c("corr", "other", "npl"),
  num_options = c("n_corr", "n_other", "n_npl"),
  version     = "ss",
  choice_rule = "softmax"
)

# NOTE (bmm#366): the recalibrated bmm M3 defaults now place equal, broad priors
# on softmax a and c (normal(3, 1)) so the two choice rules already imply a
# comparable, broad prior-predictive range. We therefore use the recalibrated
# DEFAULT softmax priors here rather than the previous hand-calibrated
# low-accuracy match. Under these defaults this script reproduces the default
# choice-rule comparison from Tutorial 1 (i.e., the separate prior-matching step
# is no longer needed). Effect priors are the shared normal(0, 0.5).
priors_matched <- c(
  prior(normal(3, 1),   class = "b", nlpar = "a", coef = "expopenset"),
  prior(normal(3, 1),   class = "b", nlpar = "a", coef = "expclosedset"),
  prior(normal(0, 0.5), class = "b", nlpar = "a", coef = "expopenset:ss_lin"),
  prior(normal(0, 0.5), class = "b", nlpar = "a", coef = "expclosedset:ss_lin"),
  prior(normal(3, 1),   class = "b", nlpar = "c", coef = "expopenset"),
  prior(normal(3, 1),   class = "b", nlpar = "c", coef = "expclosedset"),
  prior(normal(0, 0.5), class = "b", nlpar = "c", coef = "expopenset:ss_lin"),
  prior(normal(0, 0.5), class = "b", nlpar = "c", coef = "expclosedset:ss_lin")
)

###############################################################################!
# 2) Model Fitting ------------------------------------------------------------
###############################################################################!

# Prior-matched softmax, long-run (for bridge sampling)
fit_softmax_matched_lr <- bmm(
  formula      = m3_formula,
  model        = m3_softmax,
  data         = data_agg,
  prior        = priors_matched,
  chains       = chains,
  cores        = cores,
  warmup       = warmup,
  iter         = iter,
  sample_prior = "yes",
  save_pars    = save_pars(all = TRUE),
  backend      = "cmdstanr",
  file         = here("output", "fit_m3_ss_softmax_priormatched_longrun")
)

# Simple long-run fit (produced by tutorial1_simple_span.R)
fit_simple_lr <- readRDS(here("output", "fit_m3_ss_simple_longrun.rds"))

###############################################################################!
# 3) Evaluation ---------------------------------------------------------------
###############################################################################!

## 3.1) Confirm the prior-matched softmax predicts P(correct) ~ .55 ------------
prior_draws <- as_draws_df(fit_softmax_matched_lr) |>
  as_tibble() |>
  transmute(a = prior_b_a_expopenset, c = prior_b_c_expopenset)
n_other <- 2; n_npl <- 3
pp <- prior_draws |>
  mutate(den = exp(a + c) + n_other * exp(a) + n_npl,
         p_corr = exp(a + c) / den) |>
  pull(p_corr)
cat("Prior-matched softmax P(correct) [5/50/95]:\n")
print(round(quantile(pp, c(.05, .5, .95), na.rm = TRUE), 3))
cat("(target -- simple rule: 0.42 / 0.55 / 0.73)\n\n")

## 3.2) Bridge sampling: prior-matched softmax vs. simple ----------------------
bridge_softmax_matched <- bridge_sampler(fit_softmax_matched_lr, repetitions = 10)
bridge_simple_lr       <- bridge_sampler(fit_simple_lr, repetitions = 10)

bf_matched <- bf(bridge_softmax_matched, bridge_simple_lr)$bf
cat("Prior-matched BF (softmax vs simple), 10 reps:\n")
cat("  median:", median(bf_matched),
    " range: [", min(bf_matched), ",", max(bf_matched), "]\n")
cat("Compare to the default-prior long-run BF: ~6.2e7\n")

saveRDS(bridge_softmax_matched,
        here("output", "bridge_m3_ss_softmax_priormatched_longrun.rds"))

# Interpretation:
#  - If the matched-prior softmax still yields BF >> 1 in favor of softmax, the
#    advantage reflects likelihood fit, not the prior-predictive imbalance on
#    the correct-response proportion -- the strongest answer to Reviewer 1.
#  - If matching the priors substantially shrinks the BF, the original
#    comparison was partly prior-driven, and the manuscript's
#    "illustration, not adjudication" framing is the correct one.
