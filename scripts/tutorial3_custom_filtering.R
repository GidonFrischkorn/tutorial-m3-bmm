#' Tutorial 3: Custom M3 — Filtering with Separate Ratio Parameters
#'
#' This script demonstrates how to define and fit a custom M3 using
#' the same data as Tutorial 2 (Li, Frischkorn, & Oberauer, 2026).
#'
#' The standard cs version uses a single f parameter that attenuates both item
#' memory and context binding equally for distractors. The custom version here
#' replaces f with two separate ratio parameters:
#'   - ra: proportion of item memory activation transferred to distractors
#'   - rc: proportion of context binding activation transferred to distractors
#'
#' This allows testing whether filtering operates differently on item memory
#' vs. context binding — the key research question in Li et al. (2026).
#'
#' In this script, you will see:
#'  1) How to define a custom M3 with user-specified activation formulas
#'  2) How to specify link functions for bounded parameters
#'  3) How to handle non-identified parameters in the control condition
#'  4) How to compare the custom version with the standard cs version

###############################################################################!
# 0) R Setup: Packages & Data -------------------------------------------------
###############################################################################!

pacman::p_load(here, bmm, brms, tidyverse, tidybayes, patchwork)
source(here("functions", "clean_plot.R"))

# -- Fitting settings ----------------------------------------------------------
chains <- 4
cores  <- 4
warmup <- 1000
iter   <- 2000

## 0.1) Load pre-aggregated data -----------------------------------------------

# Same data as Tutorial 2 (see scripts/prepare_Li2026_data.R for details).
data_agg <- read_csv(here("data", "Li_2026_ComplexSpan_Exp1_agg.csv"),
                     show_col_types = FALSE) %>%
  mutate(condition = factor(condition, levels = c("control", "pre", "retro")))

# 135 rows: 45 participants × 3 conditions (within-subjects)
str(data_agg)

###############################################################################!
# 1) Model Specification -------------------------------------------------------
###############################################################################!

## 1.1) Define the custom M3 object --------------------------------------------

# The custom version uses the same five response categories as cs,
# but replaces the single f parameter with two ratio parameters (ra, rc).
#
# Activation formulas:
#   correct = b + a + c          (full item memory + context binding)
#   distc   = b + ra*a + rc*c    (distractor: separate attenuation for a and c)
#   other   = b + a              (other memoranda: item memory only)
#   disto   = b + ra*a           (other distractors: attenuated item memory)
#   npl     = b                  (lures: baseline only)
#
# Link functions:
#   a, c  -> identity (activations are unbounded)
#   ra, rc -> logit   (ratios bounded between 0 and 1)
m3_model_custom <- m3(
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

## 1.2) Define activation formulas and predictor formulas ----------------------

# For custom versions, bmf() specifies both the activation formulas (how
# parameters combine into category activations) and the regression formulas
# (how experimental conditions predict each parameter).
m3_formula_custom <- bmf(
  # Activation formulas (one per response category)
  corr  ~ a + c + b,
  distc ~ ra * a + rc * c + b,
  other ~ a + b,
  disto ~ ra * a + b,
  npl   ~ b,
  # Predictor formulas (condition effects + random effects)
  a  ~ 0 + condition + (0 + condition || participant),
  c  ~ 0 + condition + (0 + condition || participant),
  ra ~ 0 + condition + (0 + condition || participant),
  rc ~ 0 + condition + (0 + condition || participant)
)

## 1.3) Set priors -------------------------------------------------------------

# In the control condition, no distractors are present (n_distc = 0,
# n_disto = 0), so ra and rc do not influence the likelihood and are not
# identified. We fix both the fixed effects and their random effect SDs
# at constant values.
priors_custom <- c(
  prior(constant(0), nlpar = "ra", coef = "conditioncontrol"),
  prior(constant(0), class = "sd", nlpar = "ra",
        coef = "conditioncontrol", group = "participant"),
  prior(constant(0), nlpar = "rc", coef = "conditioncontrol"),
  prior(constant(0), class = "sd", nlpar = "rc",
        coef = "conditioncontrol", group = "participant")
)

# Verify priors
default_prior(m3_formula_custom, m3_model_custom,
              data = data_agg, prior = priors_custom)

###############################################################################!
# 2) Model Fitting -------------------------------------------------------------
###############################################################################!

fit_m3_custom <- bmm(
  formula      = m3_formula_custom,
  model        = m3_model_custom,
  data         = data_agg,
  prior        = priors_custom,
  chains       = chains,
  cores        = cores,
  warmup       = warmup,
  iter         = iter,
  sample_prior = "yes",
  save_pars    = save_pars(all = TRUE),
  backend      = "cmdstanr",
  file         = here("output", "fit_m3_custom_filtering")
)

## 2.1) Convergence checks -----------------------------------------------------
summary(fit_m3_custom)

###############################################################################!
# 3) Model Evaluation ----------------------------------------------------------
###############################################################################!

## 3.1) Posterior Predictive Checks --------------------------------------------

pp_pred <- posterior_epred(fit_m3_custom)

pred_means <- apply(pp_pred, c(2, 3), mean)
colnames(pred_means) <- c("pred_corr", "pred_distc", "pred_other",
                           "pred_disto", "pred_npl")

pp_data <- bind_cols(data_agg, as_tibble(pred_means)) %>%
  mutate(
    obs_total  = corr + other + distc + disto + npl,
    pred_total = pred_corr + pred_other + pred_distc + pred_disto + pred_npl
  )

pp_long <- pp_data %>%
  transmute(
    participant, condition,
    Observed_Correct  = corr     / obs_total,
    Observed_Other    = other    / obs_total,
    Observed_DistC    = distc   / obs_total,
    Observed_DistO    = disto   / obs_total,
    Observed_NPL      = npl      / obs_total,
    Predicted_Correct = pred_corr   / pred_total,
    Predicted_Other   = pred_other  / pred_total,
    Predicted_DistC   = pred_distc / pred_total,
    Predicted_DistO   = pred_disto / pred_total,
    Predicted_NPL     = pred_npl    / pred_total
  ) %>%
  pivot_longer(
    cols      = Observed_Correct:Predicted_NPL,
    names_to  = c("source", "category"),
    names_sep = "_",
    values_to = "proportion"
  ) %>%
  mutate(
    source   = factor(source, levels = c("Observed", "Predicted")),
    category = factor(category,
                      levels = c("Correct", "Other", "DistC", "DistO", "NPL"),
                      labels = c("Correct", "Other", "Dist (paired)",
                                 "Dist (other)", "NPL"))
  )

pp_summary <- pp_long %>%
  group_by(condition, source, category) %>%
  summarise(
    mean = mean(proportion),
    se   = sd(proportion) / sqrt(n()),
    .groups = "drop"
  )

pp_plot <- ggplot(pp_summary,
                  aes(x = category, y = mean,
                      color = source, shape = source, group = source)) +
  geom_point(size = 2.5, position = position_dodge(0.4)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = 0.2, position = position_dodge(0.4)) +
  facet_wrap(~ condition) +
  scale_color_m3() +
  labs(x = "Response Category", y = "Response Proportion",
       color = NULL, shape = NULL) +
  clean_plot(axis.text.x = element_text(angle = 30, hjust = 1))

pp_plot

ggsave(here("figures", "tutorial3_pp_check.pdf"),
       pp_plot, width = 6.5, height = 4)

## 3.2) Parameter Estimates ----------------------------------------------------

draws <- as_draws_df(fit_m3_custom) %>% as_tibble()

# a and c: identity link, no transformation needed
param_draws_ac <- tibble(
  .draw = seq_len(nrow(draws)),
  a_control = draws$b_a_conditioncontrol,
  a_pre     = draws$b_a_conditionpre,
  a_retro   = draws$b_a_conditionretro,
  c_control = draws$b_c_conditioncontrol,
  c_pre     = draws$b_c_conditionpre,
  c_retro   = draws$b_c_conditionretro
)

param_long_ac <- param_draws_ac %>%
  pivot_longer(
    cols      = a_control:c_retro,
    names_to  = c("parameter", "condition"),
    names_sep = "_",
    values_to = "activation"
  ) %>%
  mutate(
    parameter = factor(parameter, levels = c("a", "c"),
                       labels = c("a (item memory)", "c (context binding)")),
    condition = factor(condition, levels = c("control", "pre", "retro"))
  )

param_summary_ac <- param_long_ac %>%
  group_by(parameter, condition) %>%
  mean_hdci(activation, .width = 0.95)

param_plot_ac <- ggplot(param_summary_ac,
                        aes(x = condition, y = activation,
                            color = condition)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = .lower, ymax = .upper),
                width = 0.15) +
  facet_wrap(~ parameter, scales = "free_y") +
  scale_color_m3() +
  labs(x = "Condition", y = "Activation") +
  clean_plot(legend.position = "none")

# ra and rc: logit link -> apply plogis() to get probability scale
# Only pre-cue and retro-cue are interpretable; control is fixed.
param_draws_ratio <- tibble(
  .draw = seq_len(nrow(draws)),
  ra_pre   = plogis(draws$b_ra_conditionpre),
  ra_retro = plogis(draws$b_ra_conditionretro),
  rc_pre   = plogis(draws$b_rc_conditionpre),
  rc_retro = plogis(draws$b_rc_conditionretro)
)

param_long_ratio <- param_draws_ratio %>%
  pivot_longer(
    cols      = ra_pre:rc_retro,
    names_to  = c("parameter", "condition"),
    names_sep = "_",
    values_to = "ratio"
  ) %>%
  mutate(
    parameter = factor(parameter, levels = c("ra", "rc"),
                       labels = c("ra (item ratio)", "rc (binding ratio)")),
    condition = factor(condition, levels = c("pre", "retro"))
  )

param_summary_ratio <- param_long_ratio %>%
  group_by(parameter, condition) %>%
  mean_hdci(ratio, .width = 0.95)

param_plot_ratio <- ggplot(param_summary_ratio,
                           aes(x = condition, y = ratio,
                               color = condition)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = .lower, ymax = .upper),
                width = 0.15) +
  facet_wrap(~ parameter) +
  scale_color_m3() +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Condition", y = "Ratio (probability scale)") +
  clean_plot(legend.position = "none")

# Combined parameter plot
param_plot_combined <- param_plot_ac + param_plot_ratio +
  plot_layout(widths = c(2, 2))

param_plot_combined

ggsave(here("figures", "tutorial3_param_estimates.pdf"),
       param_plot_combined, width = 6.5, height = 4)

## 3.3) Model Comparison -------------------------------------------------------

# Compare the standard cs version (Tutorial 2) with the custom filtering
# version using bridge sampling (requires save_pars(all = TRUE) in both fits).
fit_m3_cs <- readRDS(here("output", "fit_m3_cs.rds"))

bridge_cs     <- bridgesampling::bridge_sampler(fit_m3_cs, repetitions = 10)
bridge_custom <- bridgesampling::bridge_sampler(fit_m3_custom, repetitions = 10)

bf_result <- bridgesampling::bf(bridge_custom, bridge_cs)
bf_result

# Save bridge sampling results for inline reporting
saveRDS(bridge_cs, here("output", "bridge_m3_cs.rds"))
saveRDS(bridge_custom, here("output", "bridge_m3_custom_filtering.rds"))

## 3.4) Hypothesis Tests -------------------------------------------------------

# Does ra differ from rc within each condition?
# If filtering is symmetric, ra = rc (equivalent to the cs assumption).
h_ra_rc_pre   <- hypothesis(fit_m3_custom,
  "ra_conditionpre - rc_conditionpre = 0")
h_ra_rc_retro <- hypothesis(fit_m3_custom,
  "ra_conditionretro - rc_conditionretro = 0")

h_ra_rc_pre
h_ra_rc_retro

# Does ra differ between pre-cue and retro-cue?
h_ra_pre_retro <- hypothesis(fit_m3_custom,
  "ra_conditionpre - ra_conditionretro = 0")

# Does rc differ between pre-cue and retro-cue?
h_rc_pre_retro <- hypothesis(fit_m3_custom,
  "rc_conditionpre - rc_conditionretro = 0")

h_ra_pre_retro
h_rc_pre_retro
