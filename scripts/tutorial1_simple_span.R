#' Tutorial 1: Simple Span M3 Model with bmm
#'
#' This script demonstrates the complete M3 workflow using data from:
#' Oberauer, K. (2019). Working Memory Capacity Limits Memory for Bindings.
#' Journal of Cognition, 2(1), 40. https://doi.org/10/gf8p8c
#' Data available at: https://osf.io/qy5sd/
#'
#' The paper contains two between-subjects experiments that differ in stimulus
#' set type:
#'  - Experiment 1 ("open set"): new word pairs on each trial
#'  - Experiment 2 ("closed set"): same small set of words reused across trials
#' Set size (2, 4, 6, 8) is a within-subjects factor in both experiments.
#'
#' In this script, you will see:
#'  1) How to prepare trial-level data for the M3 simple span model
#'  2) How to specify and fit the M3 model with bmm (simple and softmax rules)
#'  3) How to evaluate model fit with posterior predictive checks
#'  4) How to extract, interpret, and plot parameter estimates

###############################################################################!
# 0) R Setup: Packages & Data -------------------------------------------------
###############################################################################!

pacman::p_load(here, bmm, brms, bridgesampling, tidyverse, tidybayes, patchwork, gghalves)
source(here("functions", "clean_plot.R"))

# -- Fitting settings ----------------------------------------------------------
chains <- 4
cores  <- 4
warmup <- 1000
iter   <- 2000

## 0.1) Read trial-level data --------------------------------------------------

# Raw data: space-delimited, no header. Column names from the original JAGS
# analysis script (https://osf.io/qy5sd/files/geb6v).
# Key columns: setsize (2/4/6/8), rcat (1=correct, 2=other, 3=NPL),
# rsizeList (list items in recognition set, *includes* the correct item),
# rsizeNPL (not-presented lures in recognition set).
col_names <- c("id", "session", "block", "trial", "setsize",
               "rsizeList", "rsizeNPL", "tested", "response", "rcat", "rt")

# Read both experiments (between-subjects). Participant IDs overlap across
# files, so we create unique IDs by combining ID with experiment label.
read_exp <- function(file, exp_label) {
  read.table(
    here("data", file),
    header    = FALSE,
    col.names = col_names,
    fill      = TRUE      # handle trailing whitespace at end of file
  ) %>%
    drop_na() %>%          # remove any incomplete rows
    mutate(exp = exp_label,
           id  = paste(id, exp_label, sep = "_"))
}

data_raw <- bind_rows(
  read_exp("Oberauer_2019_SimpleSpan_Exp1.dat", "openset"),
  read_exp("Oberauer_2019_SimpleSpan_Exp2.dat", "closedset")
)

# Quick overview of the raw data
str(data_raw)

## 0.2) Filter to recognition trials -------------------------------------------

# Keep only recognition trials (rsizeList + rsizeNPL > 0); drop recall trials.
data_recog <- data_raw %>%
  filter(rsizeList + rsizeNPL > 0)

stopifnot(all(data_recog$rcat %in% c(1, 2, 3)))

## 0.3) Verify response set constraints ----------------------------------------

stopifnot(all(data_recog$rsizeList <= data_recog$setsize))

# 24 unique recognition conditions (setsize × rsizeList × rsizeNPL)
recog_conditions <- data_recog %>%
  distinct(setsize, rsizeList, rsizeNPL) %>%
  arrange(setsize, rsizeList, rsizeNPL)
print(recog_conditions)

## 0.4) Aggregate to response frequencies --------------------------------------

# bmm expects aggregated data: count response frequencies per participant ×
# experiment × recognition condition.
data_agg <- data_recog %>%
  group_by(id, exp, setsize, rsizeList, rsizeNPL) %>%
  summarise(
    corr  = sum(rcat == 1),   # correct responses
    other = sum(rcat == 2),   # other list item intrusions
    npl   = sum(rcat == 3),   # not-presented lure responses
    .groups = "drop"
  )

## 0.5) Compute response option counts -----------------------------------------

# Number of response options per category (required by bmm).
# n_other = rsizeList - 1 because rsizeList includes the correct item.
data_agg <- data_agg %>%
  mutate(
    n_corr  = 1,
    n_other = rsizeList - 1,
    n_npl   = rsizeNPL
  )

stopifnot(all(data_agg$n_other >= 0))
stopifnot(all(data_agg$n_npl >= 0))
stopifnot(all(data_agg$other[data_agg$n_other == 0] == 0))

## 0.6) Prepare condition factors -----------------------------------------------

# Create set size codings and convert experiment to factor.
data_agg <- data_agg %>%
  mutate(
    ss_fac = factor(setsize),
    ss_lin = setsize - mean(setsize),  # centered numeric version for linear effects)
    exp     = factor(exp)
  )

# 960 rows: 40 participants (20 per experiment) × 24 recognition conditions
str(data_agg)
data_agg %>% distinct(id, exp) %>% count(exp) %>% print()

## 0.7) Save prepared data -----------------------------------------------------

# Save the aggregated data for reuse without re-running the full pipeline.
data_agg %>%
  select(
    # ID variables
    id, exp,
    # Predictors (all set size codings)
    setsize, ss_fac, ss_lin,
    # Response frequencies (dependent variables)
    corr, other, npl,
    # Response option counts (required by m3)
    n_corr, n_other, n_npl
  ) %>%
  write_csv(here("data", "Oberauer_2019_SimpleSpan_agg.csv"))

###############################################################################!
# 1) Model Specification -------------------------------------------------------
###############################################################################!

## 1.1) Define the M3 model objects --------------------------------------------

# Simple span M3 activations: correct = b + a + c, other = b + a, npl = b.
# Choice rules convert activations to probabilities:
#   "simple":  P(i) = act(i) / sum(act)       — b fixed at 0.1
#   "softmax": P(i) = exp(act(i)) / sum(exp)   — b fixed at 0 (exp(0) = 1)
# We create both model objects to compare choice rules.
m3_model_softmax <- m3(
  resp_cats   = c("corr", "other", "npl"),
  num_options = c("n_corr", "n_other", "n_npl"),
  version     = "ss",
  choice_rule = "softmax"
)

m3_model_simple <- m3(
  resp_cats   = c("corr", "other", "npl"),
  num_options = c("n_corr", "n_other", "n_npl"),
  version     = "ss",
  choice_rule = "simple"
)

m3_model_softmax$parameters
m3_model_softmax$fixed_parameters
m3_model_softmax$links

## 1.2) Define the predictor formulas ------------------------------------------

# Linear set size effects on a and c, mirroring the original JAGS analysis
# (meanA + dA * (setsize - 5)). For a single experiment:
#   bmf(c ~ 1 + ss_lin + (1 + ss_lin | id),
#       a ~ 1 + ss_lin + (1 + ss_lin | id))
#
# With two experiments (between-subjects), we use cell means coding (0 + exp)
# with separate slopes (exp:ss_lin) and random effect variances (gr(id, by = exp)).
#
# Both choice rules use the same formula but different default links:
#   softmax → identity link (parameters are on native activation scale)
#   simple  → log link (exponentiate to get activation values)
m3_formula <- bmf(
  c ~ 0 + exp + exp:ss_lin + (1 + ss_lin | gr(id, by = exp)),
  a ~ 0 + exp + exp:ss_lin + (1 + ss_lin | gr(id, by = exp))
)

## 1.3) Check default priors ---------------------------------------------------
default_prior(m3_formula, m3_model_softmax, data = data_agg)

###############################################################################!
# 2) Model Fitting -------------------------------------------------------------
###############################################################################!

# Softmax is the primary model; simple is fitted for comparison.
# save_pars(all = TRUE) is required for the Bayes factor comparison (Section 3.4).

## 2.1) Softmax model (primary) ------------------------------------------------
fit_m3_ss_softmax <- bmm(
  formula      = m3_formula,
  model        = m3_model_softmax,
  data         = data_agg,
  chains       = chains,
  cores        = cores,
  warmup       = warmup,
  iter         = iter,
  sample_prior = "yes",
  save_pars    = save_pars(all = TRUE),
  backend      = "cmdstanr",
  file         = here("output", "fit_m3_ss_softmax")
)

## 2.2) Simple model (comparison) ----------------------------------------------
fit_m3_ss_simple <- bmm(
  formula      = m3_formula,
  model        = m3_model_simple,
  data         = data_agg,
  chains       = chains,
  cores        = cores,
  warmup       = warmup,
  iter         = iter,
  sample_prior = "yes",
  save_pars    = save_pars(all = TRUE),
  backend      = "cmdstanr",
  file         = here("output", "fit_m3_ss_simple")
)

## 2.3) Convergence checks -----------------------------------------------------
# All Rhat values should be < 1.01 and ESS should be adequate
summary(fit_m3_ss_softmax)
summary(fit_m3_ss_simple)


###############################################################################!
# 3) Model Evaluation ----------------------------------------------------------
###############################################################################!

## 3.1) Posterior Predictive Checks --------------------------------------------

# posterior_epred() returns expected counts per category per observation.
# We average across draws, then aggregate to participant × setsize × experiment.
pp_pred <- posterior_epred(fit_m3_ss_softmax)

pred_means <- apply(pp_pred, c(2, 3), mean)
colnames(pred_means) <- c("pred_corr", "pred_other", "pred_npl")

pp_data <- bind_cols(data_agg, as_tibble(pred_means)) %>%
  group_by(id, exp, setsize) %>%
  summarise(
    obs_corr  = sum(corr),  obs_other  = sum(other),  obs_npl  = sum(npl),
    pred_corr = sum(pred_corr), pred_other = sum(pred_other), pred_npl = sum(pred_npl),
    .groups = "drop"
  ) %>%
  mutate(
    obs_total  = obs_corr  + obs_other  + obs_npl,
    pred_total = pred_corr + pred_other + pred_npl
  )

pp_long <- pp_data %>%
  transmute(
    id, exp, setsize,
    Observed_Correct  = obs_corr  / obs_total,
    Observed_Other    = obs_other / obs_total,
    Observed_NPL      = obs_npl   / obs_total,
    Predicted_Correct = pred_corr  / pred_total,
    Predicted_Other   = pred_other / pred_total,
    Predicted_NPL     = pred_npl   / pred_total
  ) %>%
  pivot_longer(
    cols      = Observed_Correct:Predicted_NPL,
    names_to  = c("source", "category"),
    names_sep = "_",
    values_to = "proportion"
  ) %>%
  mutate(
    source   = factor(source, levels = c("Observed", "Predicted")),
    category = factor(category, levels = c("Correct", "Other", "NPL")),
    exp      = factor(exp, levels = c("openset", "closedset"),
                      labels = c("Open Set", "Closed Set"))
  )

pp_summary <- pp_long %>%
  group_by(exp, setsize, source, category) %>%
  summarise(
    mean = mean(proportion),
    se   = sd(proportion) / sqrt(n()),
    .groups = "drop"
  )

pp_plot <- ggplot(pp_summary,
                  aes(x = factor(setsize), y = mean,
                      color = source, shape = source, group = source)) +
  geom_point(size = 2.5, position = position_dodge(0.3)) +
  geom_line(position = position_dodge(0.3)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = 0.15, position = position_dodge(0.3)) +
  facet_grid(category ~ exp, scales = "free_y") +
  scale_color_m3() +
  labs(x = "Set Size", y = "Response Proportion",
       color = NULL, shape = NULL) +
  clean_plot()

pp_plot

ggsave(here("figures", "tutorial1_pp_check_softmax.pdf"),
       pp_plot, width = 6.5, height = 6)

# PP check for the simple choice rule model
pp_pred_simple <- posterior_epred(fit_m3_ss_simple)
pred_means_simple <- apply(pp_pred_simple, c(2, 3), mean)
colnames(pred_means_simple) <- c("pred_corr", "pred_other", "pred_npl")

pp_data_simple <- bind_cols(data_agg, as_tibble(pred_means_simple)) %>%
  group_by(id, exp, setsize) %>%
  summarise(
    obs_corr  = sum(corr),  obs_other  = sum(other),  obs_npl  = sum(npl),
    pred_corr = sum(pred_corr), pred_other = sum(pred_other), pred_npl = sum(pred_npl),
    .groups = "drop"
  ) %>%
  mutate(
    obs_total  = obs_corr  + obs_other  + obs_npl,
    pred_total = pred_corr + pred_other + pred_npl
  )

pp_long_simple <- pp_data_simple %>%
  transmute(
    id, exp, setsize,
    Observed_Correct  = obs_corr  / obs_total,
    Observed_Other    = obs_other / obs_total,
    Observed_NPL      = obs_npl   / obs_total,
    Predicted_Correct = pred_corr  / pred_total,
    Predicted_Other   = pred_other / pred_total,
    Predicted_NPL     = pred_npl   / pred_total
  ) %>%
  pivot_longer(
    cols      = Observed_Correct:Predicted_NPL,
    names_to  = c("source", "category"),
    names_sep = "_",
    values_to = "proportion"
  ) %>%
  mutate(
    source   = factor(source, levels = c("Observed", "Predicted")),
    category = factor(category, levels = c("Correct", "Other", "NPL")),
    exp      = factor(exp, levels = c("openset", "closedset"),
                      labels = c("Open Set", "Closed Set"))
  )

pp_summary_simple <- pp_long_simple %>%
  group_by(exp, setsize, source, category) %>%
  summarise(
    mean = mean(proportion),
    se   = sd(proportion) / sqrt(n()),
    .groups = "drop"
  )

pp_plot_simple <- ggplot(pp_summary_simple,
                         aes(x = factor(setsize), y = mean,
                             color = source, shape = source, group = source)) +
  geom_point(size = 2.5, position = position_dodge(0.3)) +
  geom_line(position = position_dodge(0.3)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = 0.15, position = position_dodge(0.3)) +
  facet_grid(category ~ exp, scales = "free_y") +
  scale_color_m3() +
  labs(x = "Set Size", y = "Response Proportion",
       color = NULL, shape = NULL) +
  clean_plot()

pp_plot_simple

ggsave(here("figures", "tutorial1_pp_check_simple.pdf"),
       pp_plot_simple, width = 6.5, height = 6)

## 3.2) Parameter Estimates ----------------------------------------------------

# Compute predicted activation at each set size × experiment.
# Softmax uses identity link: parameters are on native scale, no exp() needed.
fixef(fit_m3_ss_softmax)

draws <- as_draws_df(fit_m3_ss_softmax) %>% as_tibble()
ss_vals <- c(2, 4, 6, 8)
ss_centered <- ss_vals - mean(data_agg$setsize)
param_draws <- map_dfr(seq_along(ss_vals), function(i) {
  tibble(
    .draw       = seq_len(nrow(draws)),
    setsize     = ss_vals[i],
    c_openset   = draws$b_c_expopenset   +
                    draws$`b_c_expopenset:ss_lin`   * ss_centered[i],
    c_closedset = draws$b_c_expclosedset +
                    draws$`b_c_expclosedset:ss_lin` * ss_centered[i],
    a_openset   = draws$b_a_expopenset   +
                    draws$`b_a_expopenset:ss_lin`   * ss_centered[i],
    a_closedset = draws$b_a_expclosedset +
                    draws$`b_a_expclosedset:ss_lin` * ss_centered[i]
  )
})

param_long <- param_draws %>%
  pivot_longer(
    cols      = c_openset:a_closedset,
    names_to  = c("parameter", "experiment"),
    names_sep = "_",
    values_to = "activation"
  ) %>%
  mutate(
    parameter  = factor(parameter, levels = c("c", "a"),
                        labels = c("c (context binding)", "a (item memory)")),
    experiment = factor(experiment, levels = c("openset", "closedset"),
                        labels = c("Open Set", "Closed Set"))
  )

param_summary <- param_long %>%
  group_by(setsize, parameter, experiment) %>%
  mean_hdci(activation, .width = 0.95)

param_plot <- ggplot(param_summary,
                     aes(x = factor(setsize), y = activation,
                         color = experiment, group = experiment)) +
  geom_point(size = 2.5, position = position_dodge(0.3)) +
  geom_line(position = position_dodge(0.3)) +
  geom_errorbar(aes(ymin = .lower, ymax = .upper),
                width = 0.15, position = position_dodge(0.3)) +
  facet_wrap(~ parameter, scales = "free_y") +
  scale_color_m3() +
  labs(x = "Set Size", y = "Activation",
       color = "Experiment") +
  clean_plot()

param_plot

ggsave(here("figures", "tutorial1_param_estimates_softmax.pdf"),
       param_plot, width = 6.5, height = 4)

# Simple model: log link — exponentiate to get activation values.
draws_simple <- as_draws_df(fit_m3_ss_simple) %>% as_tibble()

param_draws_simple <- map_dfr(seq_along(ss_vals), function(i) {
  tibble(
    .draw       = seq_len(nrow(draws_simple)),
    setsize     = ss_vals[i],
    c_openset = exp(draws_simple$b_c_expopenset +
                        draws_simple$`b_c_expopenset:ss_lin` *
                          ss_centered[i]),
    c_closedset = exp(draws_simple$b_c_expclosedset +
                          draws_simple$`b_c_expclosedset:ss_lin` *
                            ss_centered[i]),
    a_openset = exp(draws_simple$b_a_expopenset +
                        draws_simple$`b_a_expopenset:ss_lin` *
                          ss_centered[i]),
    a_closedset = exp(draws_simple$b_a_expclosedset +
                          draws_simple$`b_a_expclosedset:ss_lin` *
                            ss_centered[i])
  )
})

param_long_simple <- param_draws_simple %>%
  pivot_longer(
    cols      = c_openset:a_closedset,
    names_to  = c("parameter", "experiment"),
    names_sep = "_",
    values_to = "activation"
  ) %>%
  mutate(
    parameter  = factor(parameter, levels = c("c", "a"),
                        labels = c("c (context binding)", "a (item memory)")),
    experiment = factor(experiment, levels = c("openset", "closedset"),
                        labels = c("Open Set", "Closed Set"))
  )

param_summary_simple <- param_long_simple %>%
  group_by(setsize, parameter, experiment) %>%
  mean_hdci(activation, .width = 0.95)

param_plot_simple <- ggplot(param_summary_simple,
                            aes(x = factor(setsize), y = activation,
                                color = experiment, group = experiment)) +
  geom_point(size = 2.5, position = position_dodge(0.3)) +
  geom_line(position = position_dodge(0.3)) +
  geom_errorbar(aes(ymin = .lower, ymax = .upper),
                width = 0.15, position = position_dodge(0.3)) +
  facet_wrap(~ parameter, scales = "free_y") +
  scale_color_m3() +
  labs(x = "Set Size", y = "Activation",
       color = "Experiment") +
  clean_plot()

param_plot_simple

ggsave(here("figures", "tutorial1_param_estimates_simple.pdf"),
       param_plot_simple, width = 6.5, height = 4)

## 3.3) Hypothesis Tests -------------------------------------------------------

# hypothesis() with sample_prior = "yes" provides Evidence Ratios (Bayes Factors).
# Softmax parameters are on native activation scale (identity link).

# H1–H2: Does c differ from a within each experiment?
h_ca_openset <- hypothesis(fit_m3_ss_softmax,
                           "c_expopenset - a_expopenset = 0")
h_ca_closedset <- hypothesis(fit_m3_ss_softmax,
                             "c_expclosedset - a_expclosedset = 0")

h_ca_openset
h_ca_closedset

# H3–H4: Overall set size effects on c and a (averaged across experiments)
h_c_ss <- hypothesis(fit_m3_ss_softmax,
  "(c_expopenset:ss_lin + c_expclosedset:ss_lin) / 2 = 0")
h_a_ss <- hypothesis(fit_m3_ss_softmax,
  "(a_expopenset:ss_lin + a_expclosedset:ss_lin) / 2 = 0")

h_c_ss
h_a_ss

# H5–H6: Do set size effects on c and a differ between experiments?
h_c_ss_exp <- hypothesis(fit_m3_ss_softmax,
  "c_expopenset:ss_lin - c_expclosedset:ss_lin = 0")
h_a_ss_exp <- hypothesis(fit_m3_ss_softmax,
  "a_expopenset:ss_lin - a_expclosedset:ss_lin = 0")

h_c_ss_exp
h_a_ss_exp

## 3.4) Choice Rule Comparison -------------------------------------------------

# Bayes factor via bridge sampling (requires save_pars(all = TRUE) from Section 2).
# 10 repetitions quantify variability of the marginal likelihood estimates.
bridge_file_softmax <- here("output", "bridge_m3_ss_softmax.rds")
bridge_file_simple  <- here("output", "bridge_m3_ss_simple.rds")

if (file.exists(bridge_file_softmax) & file.exists(bridge_file_simple)) {
  bridge_softmax <- readRDS(bridge_file_softmax)
  bridge_simple  <- readRDS(bridge_file_simple)
} else {
  bridge_softmax <- bridge_sampler(fit_m3_ss_softmax, repetitions = 10)
  bridge_simple  <- bridge_sampler(fit_m3_ss_simple,  repetitions = 10)
  saveRDS(bridge_softmax, bridge_file_softmax)
  saveRDS(bridge_simple,  bridge_file_simple)
}

# bf() with repetitions > 1 reports median BF and range across repetitions.
bf_choice_rule <- bf(bridge_softmax, bridge_simple)
bf_choice_rule
