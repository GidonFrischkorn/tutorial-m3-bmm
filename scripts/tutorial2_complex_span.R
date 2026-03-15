#' Tutorial 2: Complex Span M3 Model with bmm

###############################################################################!
# 0) R Setup: Packages & Data -------------------------------------------------
###############################################################################!

# -- Download pre-fitted model objects from OSF (optional) --------------------
# Set to TRUE to download all fitted model objects before running this script.
# This avoids re-fitting models from scratch (which can take several hours).
download_from_osf <- FALSE
if (download_from_osf) source(here::here("scripts", "00_download_osf.R"))

pacman::p_load(here, bmm, brms, tidyverse, tidybayes, patchwork, gghalves)
source(here("functions", "clean_plot.R"))

# -- Fitting settings ----------------------------------------------------------
chains <- 4
cores  <- 4
warmup <- 1000
iter   <- 2000

## 0.1) Load pre-aggregated data -----------------------------------------------

# Data preparation is handled by the companion script
# (scripts/prepare_Li2026_data.R). See that script for exclusion
# criteria and response classification details.
data_agg <- read_csv(here("data", "Li_2026_ComplexSpan_Exp1_agg.csv"),
                     show_col_types = FALSE)

# Convert condition to a factor with control as reference
data_agg <- data_agg |>
  mutate(condition = factor(condition, levels = c("control", "pre", "retro")))

# 135 rows: 45 participants × 3 conditions (within-subjects)
str(data_agg)
data_agg |> count(condition)

###############################################################################!
# 1) Model Specification -------------------------------------------------------
###############################################################################!

## 1.1) Define the M3 model object ---------------------------------------------

# The complex span M3 model (version = "cs") has five response categories.
# IMPORTANT: The cs version assigns activation formulas by position:
#   1. corr  = b + a + c         (target: full activation)
#   2. distc = b + f*a + f*c     (paired distractor: attenuated by f)
#   3. other = b + a             (other memoranda: item memory only)
#   4. disto = b + f*a           (other distractors: attenuated item memory)
#   5. npl   = b                 (not-presented lures: baseline only)
#
# The f parameter (0–1) captures distractor filtering efficiency:
#   f = 0 → perfect filtering (distractors have baseline activation only)
#   f = 1 → no filtering (distractors are as strong as memoranda)
m3_model_cs <- m3(
  resp_cats   = c("corr", "distc", "other", "disto", "npl"),
  num_options = c("n_corr", "n_distc", "n_other", "n_disto", "n_npl"),
  version     = "cs",
  choice_rule = "softmax"
)

## 1.2) Define the predictor formulas ------------------------------------------

# All parameters predicted by condition using cell means coding (0 + condition).
# We use condition-specific uncorrelated random effects
# (0 + condition || participant) so that each condition gets its own
# random effect variance.
m3_formula_cs <- bmf(
  c ~ 0 + condition + (0 + condition || participant),
  a ~ 0 + condition + (0 + condition || participant),
  f ~ 0 + condition + (0 + condition || participant)
)

## 1.3) Set priors -------------------------------------------------------------

# Check default priors
default_prior(m3_formula_cs, m3_model_cs, data = data_agg)

# Fix f in the control condition: there are no distractors (n_distc = 0,
# n_disto = 0), so f does not influence the likelihood and is not identified.
# We fix both the fixed effect and the random effect SD at constant(0) to
# fully remove f in the control condition from the model. Without fixing the
# SD, it samples from its wide default prior and can become very large,
# potentially affecting other parameters.
priors_cs <- c(
  prior(constant(-10), nlpar = "f", coef = "conditioncontrol"),
  prior(constant(0), class = "sd", nlpar = "f",
        coef = "conditioncontrol", group = "participant")
)

# Verify that custom priors are correctly applied
default_prior(m3_formula_cs, m3_model_cs, data = data_agg, prior = priors_cs)

###############################################################################!
# 2) Model Fitting -------------------------------------------------------------
###############################################################################!

fit_m3_cs <- bmm(
  formula      = m3_formula_cs,
  model        = m3_model_cs,
  data         = data_agg,
  prior        = priors_cs,
  chains       = chains,
  cores        = cores,
  warmup       = warmup,
  iter         = iter,
  sample_prior = "yes",
  save_pars    = save_pars(all = TRUE),
  backend      = "cmdstanr",
  file         = here("output", "fit_m3_cs")
)

## 2.1) Convergence checks -----------------------------------------------------
summary(fit_m3_cs)

###############################################################################!
# 3) Model Evaluation ----------------------------------------------------------
###############################################################################!

## 3.1) Posterior Predictive Checks --------------------------------------------

# posterior_epred() returns expected counts per category per observation.
# We average across draws, then compute proportions per participant × condition.
pp_pred <- posterior_epred(fit_m3_cs)

pred_means <- apply(pp_pred, c(2, 3), mean)
colnames(pred_means) <- c("pred_corr", "pred_distc", "pred_other",
                          "pred_disto", "pred_npl")

pp_data <- bind_cols(data_agg, as_tibble(pred_means)) |>
  mutate(
    obs_total  = corr + other + distc + disto + npl,
    pred_total = pred_corr + pred_other + pred_distc + pred_disto + pred_npl
  )

# Reshape to long format
pp_long <- pp_data |>
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
  ) |>
  pivot_longer(
    cols      = Observed_Correct:Predicted_NPL,
    names_to  = c("source", "category"),
    names_sep = "_",
    values_to = "proportion"
  ) |>
  mutate(
    source   = factor(source, levels = c("Observed", "Predicted")),
    category = factor(category,
                      levels = c("Correct", "Other", "DistC", "DistO", "NPL"),
                      labels = c("Correct", "Other", "Dist (paired)",
                                 "Dist (other)", "NPL"))
  )

pp_summary <- pp_long |>
  group_by(condition, source, category) |>
  summarise(
    mean = mean(proportion),
    se   = sd(proportion) / sqrt(n()),
    .groups = "drop"
  )

pp_plot <- ggplot(pp_summary,
                  aes(x = condition, y = mean,
                      color = source, shape = source, group = source)) +
  geom_point(size = 2.5, position = position_dodge(0.4)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = 0.2, position = position_dodge(0.4)) +
  facet_wrap(~ category, nrow = 1, scales = "free_y") +
  scale_color_m3() +
  labs(x = "Condition", y = "Response Proportion",
       color = NULL, shape = NULL) +
  clean_plot(axis.text.x = element_text(angle = 30, hjust = 1)) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(1.0, 0.85),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
    legend.key.size = unit(0.4, "lines"),
    legend.key.spacing.y = unit(0.3, "lines"),
    legend.text = element_text(size = 8)
  )

# Apply per-facet y-axis limits to highlight error categories
pp_plot <- scale_individual_facet_y_axes(
  pp_plot,
  ylims = list(
    c(0.5, 1),    # Correct
    c(0, 0.20),   # Other
    c(0, 0.20),   # Dist (paired)
    c(0, 0.20),   # Dist (other)
    c(0, 0.20)    # NPL
  )
)

pp_plot

ggsave(here("figures", "tutorial2_pp_check.pdf"),
       pp_plot, width = 6.5, height = 3)

## 3.2) Parameter Estimates ----------------------------------------------------

# Softmax cs model uses identity link for a and c, logit link for f.
# a and c are on native activation scale; f must be transformed with plogis().
fixef(fit_m3_cs)

draws <- as_draws_df(fit_m3_cs) |> as_tibble()

# Extract condition-level estimates for a and c
# (identity link, no transformation)
param_draws_ac <- tibble(
  .draw = seq_len(nrow(draws)),
  a_control = draws$b_a_conditioncontrol,
  a_pre     = draws$b_a_conditionpre,
  a_retro   = draws$b_a_conditionretro,
  c_control = draws$b_c_conditioncontrol,
  c_pre     = draws$b_c_conditionpre,
  c_retro   = draws$b_c_conditionretro
)

param_long_ac <- param_draws_ac |>
  pivot_longer(
    cols      = a_control:c_retro,
    names_to  = c("parameter", "condition"),
    names_sep = "_",
    values_to = "activation"
  ) |>
  mutate(
    parameter = factor(parameter, levels = c("a", "c"),
                       labels = c("a (item memory)", "c (context binding)")),
    condition = factor(condition, levels = c("control", "pre", "retro"))
  )

param_summary_ac <- param_long_ac |>
  group_by(parameter, condition) |>
  mean_hdci(activation, .width = 0.95)

param_plot_ac <- ggplot(param_summary_ac,
                        aes(x = condition, y = activation,
                            color = condition)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = .lower, ymax = .upper),
                width = 0.15) +
  facet_wrap(~ parameter) +
  scale_color_m3() +
  labs(x = "Condition", y = "Activation") +
  clean_plot(legend.position = "none")

param_plot_ac

# Extract f parameter (logit link → apply plogis() to get probability scale).
# Only pre-cue and retro-cue are interpretable; control is fixed.
param_draws_f <- tibble(
  .draw = seq_len(nrow(draws)),
  f_pre   = plogis(draws$b_f_conditionpre),
  f_retro = plogis(draws$b_f_conditionretro)
)

param_long_f <- param_draws_f |>
  pivot_longer(
    cols      = f_pre:f_retro,
    names_to  = "condition",
    names_prefix = "f_",
    values_to = "f"
  ) |>
  mutate(condition = factor(condition, levels = c("pre", "retro")))

param_summary_f <- param_long_f |>
  group_by(condition) |>
  mean_hdci(f, .width = 0.95)

param_plot_f <- ggplot(param_summary_f,
                       aes(x = condition, y = f, color = condition)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = .lower, ymax = .upper),
                width = 0.15) +
  scale_color_m3() +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Condition", y = "f (distractor filtering)") +
  clean_plot(legend.position = "none")

param_plot_f

# Combined parameter plot
param_plot_combined <- param_plot_ac + param_plot_f +
  plot_layout(widths = c(2, 1))

param_plot_combined

ggsave(here("figures", "tutorial2_param_estimates.pdf"),
       param_plot_combined, width = 6.5, height = 3.2)

## 3.3) Hypothesis Tests -------------------------------------------------------

# hypothesis() with sample_prior = "yes" provides Evidence
# Ratios (Bayes Factors).
# a and c are on the native activation scale (identity link).
# f is on the logit scale in the model parameterization.

# H1: Does f differ between pre-cue and retro-cue conditions?
h_f_pre_retro <- hypothesis(fit_m3_cs, "f_conditionpre - f_conditionretro = 0")
h_f_pre_retro

# H2–H3: Do a and c differ between pre-cue and retro-cue?
h_a_pre_retro <- hypothesis(fit_m3_cs, "a_conditionpre - a_conditionretro = 0")
h_c_pre_retro <- hypothesis(fit_m3_cs, "c_conditionpre - c_conditionretro = 0")

h_a_pre_retro
h_c_pre_retro

# H4–H5: Do a and c differ between control and distractor conditions?
# Average of pre and retro vs. control.
h_a_ctrl_dist <- hypothesis(
  fit_m3_cs,
  "a_conditioncontrol - (a_conditionpre + a_conditionretro) / 2 = 0"
)
h_c_ctrl_dist <- hypothesis(
  fit_m3_cs,
  "c_conditioncontrol - (c_conditionpre + c_conditionretro) / 2 = 0"
)

h_a_ctrl_dist
h_c_ctrl_dist