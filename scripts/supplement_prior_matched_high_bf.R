#' Supplement: Choice-rule Bayes factor under well-calibrated, matched priors
#'
#' Companion to scripts/supplement_prior_matched_bf.R. That script matched both
#' rules to a LOW prior predictive (~55% correct, the simple rule's default) and
#' found the simple rule favored (BF01 ~ 300). But the data show ~91% correct
#' (median), so 55% is an implausible prior-predictive target that penalizes the
#' (otherwise well-calibrated) softmax rule.
#'
#' Here we match BOTH rules to a data-PLAUSIBLE high-accuracy prior predictive
#' (median P(correct) ~ 0.85, 90% CI ~ 0.62-0.95) so the comparison is on equal,
#' well-calibrated footing. Intercept priors calibrated by simulation to hit
#' that common target for each rule's own link structure:
#'   SOFTMAX:  a ~ normal(1.96, 0.38)   c ~ normal(2.62, 0.74)
#'   SIMPLE:   a ~ normal(0.70, 0.52)   c ~ normal(3.13, 0.75)
#'
#' Comparison summary (correct-response prior predictive in brackets):
#'   default priors        softmax[89%] vs simple[55%]  -> BF10 ~ 6.2e7  (softmax)
#'   matched-low  (this=no) both ~55%                   -> BF10 ~ 0.0034 (simple)
#'   matched-high (this)    both ~85%                   -> ?
#' If softmax wins here too, its advantage is robust to a fair, data-plausible
#' prior. If the result flips with the target, the choice-rule BF is prior-driven
#' and should be reported as illustrative rather than decisive.
#'
#' NOTE: fits TWO models at 80,000 post-warmup samples (~25 min total). Cached.

###############################################################################!
# 0) R Setup ------------------------------------------------------------------
###############################################################################!
pacman::p_load(here, bmm, brms, tidyverse, tidybayes, bridgesampling)

chains <- 4
cores  <- 4
warmup <- 1000
iter   <- 21000

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

m3_softmax <- m3(resp_cats = c("corr", "other", "npl"),
                 num_options = c("n_corr", "n_other", "n_npl"),
                 version = "ss", choice_rule = "softmax")
m3_simple  <- m3(resp_cats = c("corr", "other", "npl"),
                 num_options = c("n_corr", "n_other", "n_npl"),
                 version = "ss", choice_rule = "simple")

# Calibrated priors -> common prior predictive P(correct) ~ 0.85 for each rule
priors_softmax_hi <-
  prior(normal(1.96, 0.38), class = "b", coef = "expclosedset", nlpar = "a") +
  prior(normal(1.96, 0.38), class = "b", coef = "expopenset",   nlpar = "a") +
  prior(normal(2.62, 0.74), class = "b", coef = "expclosedset", nlpar = "c") +
  prior(normal(2.62, 0.74), class = "b", coef = "expopenset",   nlpar = "c")

priors_simple_hi <-
  prior(normal(0.70, 0.52), class = "b", coef = "expclosedset", nlpar = "a") +
  prior(normal(0.70, 0.52), class = "b", coef = "expopenset",   nlpar = "a") +
  prior(normal(3.13, 0.75), class = "b", coef = "expclosedset", nlpar = "c") +
  prior(normal(3.13, 0.75), class = "b", coef = "expopenset",   nlpar = "c")

###############################################################################!
# 2) Model Fitting (both rules, calibrated priors) ----------------------------
###############################################################################!

fit_softmax_hi <- bmm(
  formula = m3_formula, model = m3_softmax, data = data_agg,
  prior = priors_softmax_hi, chains = chains, cores = cores,
  warmup = warmup, iter = iter, sample_prior = "yes",
  save_pars = save_pars(all = TRUE), backend = "cmdstanr",
  file = here("output", "fit_m3_ss_softmax_matchedhi_longrun")
)

fit_simple_hi <- bmm(
  formula = m3_formula, model = m3_simple, data = data_agg,
  prior = priors_simple_hi, chains = chains, cores = cores,
  warmup = warmup, iter = iter, sample_prior = "yes",
  save_pars = save_pars(all = TRUE), backend = "cmdstanr",
  file = here("output", "fit_m3_ss_simple_matchedhi_longrun")
)

###############################################################################!
# 3) Evaluation ---------------------------------------------------------------
###############################################################################!

## 3.1) Confirm both prior predictives sit at ~85% correct ---------------------
n_other <- 2; n_npl <- 3
pp_soft <- as_draws_df(fit_softmax_hi) |> as_tibble() |>
  transmute(a = prior_b_a_expopenset, c = prior_b_c_expopenset) |>
  mutate(d = exp(a + c) + n_other * exp(a) + n_npl, p = exp(a + c) / d) |> pull(p)
pp_simp <- as_draws_df(fit_simple_hi) |> as_tibble() |>
  transmute(a = prior_b_a_expopenset, c = prior_b_c_expopenset) |>
  mutate(ac = 0.1 + exp(a) + exp(c), ao = 0.1 + exp(a),
         d = ac + n_other * ao + n_npl * 0.1, p = ac / d) |> pull(p)
cat("Prior-pred P(correct), softmax [5/50/95]:",
    paste(round(quantile(pp_soft, c(.05, .5, .95), na.rm = TRUE), 3), collapse = "/"), "\n")
cat("Prior-pred P(correct), simple  [5/50/95]:",
    paste(round(quantile(pp_simp, c(.05, .5, .95), na.rm = TRUE), 3), collapse = "/"), "\n\n")

## 3.2) Bridge sampling BF: matched-high softmax vs matched-high simple --------
br_soft <- bridge_sampler(fit_softmax_hi, repetitions = 10)
br_simp <- bridge_sampler(fit_simple_hi, repetitions = 10)
bf_hi <- bf(br_soft, br_simp)$bf
cat("Matched-HIGH (both ~85%) BF10 (softmax vs simple), 10 reps:\n")
cat("  median:", median(bf_hi),
    " range: [", min(bf_hi), ",", max(bf_hi), "]\n")
cat("Reference: default-prior BF10 ~ 6.2e7 (softmax);",
    "matched-low (~55%) BF10 ~ 0.0034 (simple)\n")

saveRDS(br_soft, here("output", "bridge_m3_ss_softmax_matchedhi_longrun.rds"))
saveRDS(br_simp, here("output", "bridge_m3_ss_simple_matchedhi_longrun.rds"))
