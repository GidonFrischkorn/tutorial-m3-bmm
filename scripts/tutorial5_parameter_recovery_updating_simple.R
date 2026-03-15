#' Tutorial 5: Custom M3 — Parameter Recovery (Simplified Version)

###############################################################################!
# 0) R Setup: Packages & Settings ---------------------------------------------
###############################################################################!

# -- Download pre-fitted model objects from OSF (optional) --------------------
# Set to TRUE to download all fitted model objects before running this script.
# This avoids re-fitting models from scratch (which can take several hours).
download_from_osf <- FALSE
if (download_from_osf) source(here::here("scripts", "00_download_osf.R"))

pacman::p_load(here, bmm, brms, tidyverse, tidybayes, patchwork)
source(here("functions", "clean_plot.R"))

# -- Fitting settings ----------------------------------------------------------
chains <- 4
cores  <- 4
warmup <- 1000
iter   <- 2000

# -- Reproducibility -----------------------------------------------------------
set.seed(2025)

###############################################################################!
# 1) Model Specification -------------------------------------------------------
###############################################################################!

## 1.1) Define the custom M3 model object --------------------------------------

# m3() tells bmm about the structure of your model: which response categories
# exist, how many response options each has, and the link functions for each
# parameter.
m3_model_custom <- m3(
  resp_cats = c("correct", "other", "oldinpos", "otherold", "npl"),
  num_options = c(
    "n_correct"  = 1,
    "n_other"    = 4,
    "n_oldinpos" = 1,
    "n_otherold" = 4,
    "n_npl"      = 5
  ),
  version     = "custom",
  choice_rule = "softmax",
  links = list(
    a = "identity",
    c = "identity",
    d = "logit",
    e = "log",
    r = "log"
  )
)

## 1.2) Define the activation formulas ----------------------------------------

# bmf() specifies both the activation formulas (how parameters combine for
# each response category) and the prediction formulas (what predicts each
# parameter). bmm parses exp() and arithmetic in these formulas internally.
#
# For simulation with rm3(), we only need the activation formulas:
m3_act_funs <- bmf(
  correct  ~ (1 + e * te) * (a + c) + b,
  other    ~ (1 + e * te) * a + b,
  oldinpos ~ exp(-r * tr) * d * (1 + e * te) * (a + c) + b,
  otherold ~ exp(-r * tr) * d * (1 + e * te) * a + b,
  npl      ~ b
)

## 1.3) Define the fitting formula --------------------------------------------

# For model fitting, we also need prediction formulas that tell bmm what
# predicts each parameter. Here, each parameter has:
#   - an intercept (1) = the group-level mean
#   - a random intercept (1 | participant) = person-level deviation
#
# te and tr are NOT predicted — they enter through the activation formulas
# as data variables.
m3_formula_fit <- bmf(
  correct  ~ (1 + e * te) * (a + c) + b,
  other    ~ (1 + e * te) * a + b,
  oldinpos ~ exp(-r * tr) * d * (1 + e * te) * (a + c) + b,
  otherold ~ exp(-r * tr) * d * (1 + e * te) * a + b,
  npl      ~ b,
  a ~ 1 + (1 | participant),
  c ~ 1 + (1 | participant),
  d ~ 1 + (1 | participant),
  e ~ 1 + (1 | participant),
  r ~ 1 + (1 | participant)
)

###############################################################################!
# 2) Data Simulation -----------------------------------------------------------
###############################################################################!

## 2.1) Choose generating parameter values ------------------------------------

# Group-level means (on link scale)
mean_c <- 2.0     # identity link → native = 2.0
mean_a <- 1.0     # identity link → native = 1.0
mean_d <- 0.85    # logit link → native = plogis(0.85) ≈ 0.70
mean_e <- 0.00    # log link → native = exp(0) = 1.0
mean_r <- -0.69   # log link → native = exp(-0.69) ≈ 0.50

# Between-person SDs (on link scale)
sd_c <- 0.50
sd_a <- 0.50
sd_d <- 0.65
sd_e <- 0.35
sd_r <- 0.38

## 2.2) Simulation settings ---------------------------------------------------

# We simulate a single recovery cell: one sample size and one trial count.
# The full script (tutorial3_parameter_recovery.R) repeats this across
# multiple sample sizes and trial counts.
N               <- 50    # number of participants
trials_per_cond <- 25    # trials per experimental condition

## 2.3) Draw true individual parameters ---------------------------------------

# Each simulated participant gets their own parameter values, drawn from a
# normal distribution around the group mean (on the link scale). This
# mimics real individual differences.
true_pars <- tibble(
  participant = 1:N,
  c_link = rnorm(N, mean_c, sd_c),
  a_link = rnorm(N, mean_a, sd_a),
  d_link = rnorm(N, mean_d, sd_d),
  e_link = rnorm(N, mean_e, sd_e),
  r_link = rnorm(N, mean_r, sd_r)
)

# Convert to native scale for use in rm3() (which needs native-scale values)
true_pars <- true_pars |>
  mutate(
    c_native = c_link,         # identity: no transformation
    a_native = a_link,         # identity: no transformation
    d_native = plogis(d_link), # logit → probability (0 to 1)
    e_native = exp(e_link),    # log → positive
    r_native = exp(r_link)     # log → positive
  )

# Inspect the generated individual parameters
true_pars |>
  summarise(
    across(ends_with("_native"),
           list(mean = mean, sd = sd, min = min, max = max),
           .names = "{.col}_{.fn}")
  ) |>
  pivot_longer(everything(),
               names_to = c("parameter", "stat"),
               names_pattern = "(.+)_native_(.+)") |>
  pivot_wider(names_from = stat, values_from = value) |>
  print()

## 2.4) Set up the experimental design ----------------------------------------

# 3 × 3 factorial design: encoding time (te) × removal time (tr)
# Each participant completes all 9 conditions.
design_grid <- expand_grid(
  te = c(0.25, 0.5, 1.0),
  tr = c(0.25, 0.5, 1.0)
)

# The number of response options is the same for all conditions:
#   1 correct target, 4 other current items, 1 old item in the tested
#   position, 4 old items from other positions, 5 not-presented lures
design_grid <- design_grid |>
  mutate(
    n_correct  = 1L,
    n_other    = 4L,
    n_oldinpos = 1L,
    n_otherold = 4L,
    n_npl      = 5L
  )

print(design_grid)

n_conds <- nrow(design_grid)  # = 9

## 2.5) Demonstrate rm3() for a single participant and condition ---------------

# Before simulating the full dataset, let's see what rm3() does.
# rm3() generates multinomial response counts from the M3 model.

# Parameters for participant 1 (native scale)
pars_p1 <- c(
  a = true_pars$a_native[1],
  c = true_pars$c_native[1],
  d = true_pars$d_native[1],
  e = true_pars$e_native[1],
  r = true_pars$r_native[1]
)

print(round(pars_p1, 3))

# Simulate trials for condition 1 (te = 0.25, tr = 0.25)
demo_data <- rm3(
  n        = 1,
  size     = trials_per_cond,
  pars     = c(pars_p1, b = 0, te = design_grid$te[1], tr = design_grid$tr[1]),
  m3_model = m3_model_custom,
  act_funs = m3_act_funs
)

print(data.frame(demo_data))

## 2.6) Simulate the full dataset ---------------------------------------------

data_list <- list()
row_idx   <- 1

for (i in seq_len(N)) {

  # This participant's parameters (native scale)
  pars_i <- c(
    a = true_pars$a_native[i],
    c = true_pars$c_native[i],
    d = true_pars$d_native[i],
    e = true_pars$e_native[i],
    r = true_pars$r_native[i]
  )

  for (j in seq_len(n_conds)) {

    # Simulate response counts for this participant × condition
    sim_ij <- rm3(
      n        = 1,
      size     = trials_per_cond,
      pars     = c(pars_i, b = 0,
                   te = design_grid$te[j],
                   tr = design_grid$tr[j]),
      m3_model = m3_model_custom,
      act_funs = m3_act_funs
    )

    # Convert to a plain data frame and add identifying columns
    row_ij <- data.frame(sim_ij)
    row_ij$participant <- i
    row_ij$te          <- design_grid$te[j]
    row_ij$tr          <- design_grid$tr[j]
    row_ij$n_correct   <- design_grid$n_correct[j]
    row_ij$n_other     <- design_grid$n_other[j]
    row_ij$n_oldinpos  <- design_grid$n_oldinpos[j]
    row_ij$n_otherold  <- design_grid$n_otherold[j]
    row_ij$n_npl       <- design_grid$n_npl[j]

    data_list[[row_idx]] <- row_ij
    row_idx <- row_idx + 1
  }

  # Progress message every 10 participants
  if (i %% 10 == 0) cat("  Simulated participant", i, "of", N, "\n")
}

sim_data <- bind_rows(data_list)
str(sim_data)

## 2.7) Inspect the simulated data --------------------------------------------

# Check that response proportions look reasonable.
# We expect: correct > other > oldinpos/otherold > npl
cat("\nMean response proportions across all participants and conditions:\n")
sim_data |>
  mutate(total = correct + other + oldinpos + otherold + npl) |>
  summarise(
    correct  = mean(correct / total),
    other    = mean(other / total),
    oldinpos = mean(oldinpos / total),
    otherold = mean(otherold / total),
    npl      = mean(npl / total)
  ) |>
  print()

# How do proportions change across conditions?
cat("\nMean correct proportion by encoding time (te) and removal time (tr):\n")
sim_data |>
  mutate(total = correct + other + oldinpos + otherold + npl) |>
  group_by(te, tr) |>
  summarise(prop_correct = mean(correct / total), .groups = "drop") |>
  pivot_wider(names_from = tr, values_from = prop_correct,
              names_prefix = "tr=") |>
  print()

###############################################################################!
# 3) Model Fitting -------------------------------------------------------------
###############################################################################!

## 3.1) Define the model object for fitting ------------------------------------

m3_model_custom_fit <- m3(
  resp_cats   = c("correct", "other", "oldinpos", "otherold", "npl"),
  num_options = c("n_correct", "n_other", "n_oldinpos", "n_otherold", "n_npl"),
  version     = "custom",
  choice_rule = "softmax",
  links = list(
    a = "identity",
    c = "identity",
    d = "logit",
    e = "log",
    r = "log"
  )
)

## 3.2) Check default priors ---------------------------------------------------

# Before fitting, inspect what priors bmm will use by default.
# Priors that are too wide can slow convergence; priors that are too narrow
# can bias estimates.
default_prior(m3_formula_fit, m3_model_custom_fit, data = sim_data)

## 3.3) Fit the model ----------------------------------------------------------

fit <- bmm(
  formula = m3_formula_fit,
  model   = m3_model_custom_fit,
  data    = sim_data,
  chains  = chains,
  cores   = cores,
  warmup  = warmup,
  iter    = iter,
  backend = "cmdstanr",
  file    = here("output", "fit_m3_custom_simple_N50_T25")
)

## 3.4) Convergence checks -----------------------------------------------------

summary(fit)

###############################################################################!
# 4) Recovery Evaluation -------------------------------------------------------
###############################################################################!

## 4.1) Group-level recovery ---------------------------------------------------

# fixef() extracts the group-level means (fixed effects) with 95% CIs.
# These are on the LINK scale.
fe <- fixef(fit)
print(fe)

# Build a comparison table: true vs. recovered
par_names <- c("a", "c", "d", "e", "r")
true_means <- c(mean_a, mean_c, mean_d, mean_e, mean_r)

group_recovery <- tibble(
  parameter       = par_names,
  true_mean       = true_means,
  recovered_mean  = fe[paste0(par_names, "_Intercept"), "Estimate"],
  recovered_lower = fe[paste0(par_names, "_Intercept"), "Q2.5"],
  recovered_upper = fe[paste0(par_names, "_Intercept"), "Q97.5"]
)

# Compute recovery metrics:
#   bias: recovered minus true (positive means overestimation)
#   coverage: does the 95% CI contain the true value?
group_recovery <- group_recovery |>
  mutate(
    bias     = recovered_mean - true_mean,
    coverage = true_mean >= recovered_lower & true_mean <= recovered_upper
  )

cat("\nGroup-level recovery:\n")
print(group_recovery)

cat("\nSummary: mean absolute bias =",
    round(mean(abs(group_recovery$bias)), 3),
    ", coverage =", sum(group_recovery$coverage), "/", nrow(group_recovery),
    "parameters covered\n")

## 4.2) Group-level recovery plot ----------------------------------------------

group_plot <- ggplot(group_recovery,
                     aes(x = true_mean, y = recovered_mean)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  geom_pointrange(aes(ymin = recovered_lower, ymax = recovered_upper),
                  size = 0.5) +
  geom_text(aes(label = parameter), nudge_x = 0.15, size = 3.5) +
  labs(x = "True generating value (link scale)",
       y = "Recovered estimate (link scale)",
       title = "Group-level parameter recovery") +
  clean_plot()

group_plot

ggsave(here("figures", "tutorial5_simple_group_recovery.pdf"),
       group_plot, width = 5, height = 5)

## 4.3) SD recovery ------------------------------------------------------------

# VarCorr() extracts the random effect standard deviations. These tell us
# about the between-person variability — compare to the true SDs we used
# to generate individual differences.
vc <- VarCorr(fit)

true_sds <- c(sd_a, sd_c, sd_d, sd_e, sd_r)

sd_recovery <- tibble(
  parameter    = par_names,
  true_sd      = true_sds,
  recovered_sd = NA_real_,
  sd_lower     = NA_real_,
  sd_upper     = NA_real_
)

# Extract each parameter's SD estimate from the VarCorr output
for (k in seq_along(par_names)) {
  p <- par_names[k]
  sd_est <- vc$participant$sd[paste0(p, "_Intercept"), ]
  sd_recovery$recovered_sd[k] <- sd_est["Estimate"]
  sd_recovery$sd_lower[k]     <- sd_est["Q2.5"]
  sd_recovery$sd_upper[k]     <- sd_est["Q97.5"]
}

sd_recovery <- sd_recovery |>
  mutate(
    sd_bias    = recovered_sd - true_sd,
    sd_coverage = true_sd >= sd_lower & true_sd <= sd_upper
  )

print(sd_recovery)

## 4.4) Individual-level recovery ----------------------------------------------

# coef() returns individual-level estimates: the group mean + each person's
# random effect deviation. These are on the link scale.
co <- coef(fit)

# Extract estimates for each parameter and participant
indiv_recovery <- tibble()

for (p in par_names) {
  # coef() returns a 3D array: [participant, statistic, parameter]
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

# Add the true link-scale values for comparison
true_long <- true_pars |>
  select(participant,
         a = a_link, c = c_link, d = d_link, e = e_link, r = r_link) |>
  pivot_longer(-participant, names_to = "parameter",
               values_to = "true_link")

indiv_recovery <- indiv_recovery |>
  left_join(true_long, by = c("participant", "parameter"))

# Compute individual-level recovery metrics for each parameter
indiv_summary <- indiv_recovery |>
  group_by(parameter) |>
  summarise(
    correlation = cor(true_link, recovered_link),
    mean_bias   = mean(recovered_link - true_link),
    rmse        = sqrt(mean((recovered_link - true_link)^2)),
    coverage    = mean(true_link >= indiv_lower & true_link <= indiv_upper),
    .groups     = "drop"
  )

print(indiv_summary)

## 4.5) Individual-level scatter plots -----------------------------------------

# True vs. recovered for each parameter. Points near the diagonal = good
# recovery. Spread around the diagonal = uncertainty in individual estimates.

# Prepare correlation labels for annotation in each panel
corr_labels <- indiv_summary |>
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
       y = "Recovered value (link scale)",
       title = "Individual-level parameter recovery") +
  clean_plot()

scatter_plot

ggsave(here("figures", "tutorial5_simple_indiv_scatter.pdf"),
       scatter_plot, width = 6.5, height = 5)

###############################################################################!
# 5) Next Steps ----------------------------------------------------------------
###############################################################################!

# This simplified script demonstrates the full workflow for a single cell
# (N = 50, trials_per_cond as set in Section 2.2). The full tutorial script
# (tutorial3_parameter_recovery.R) extends this by:
#
#   1) Varying sample size (N = 25, 50, 100) and trials per condition
#      (5, 10, 25) to show how data quantity affects recovery quality.
#
#   2) Simulating true parameters once for N_max = 100, then using subsets
#      for smaller sample sizes to ensure fair comparison.
#
#   3) Summarizing recovery across all 9 cells with heatmaps showing how
#      correlation and coverage improve with more data.
#
# For publication-quality simulation studies, each cell should be replicated
# many times (e.g., 100–500 replications) to account for Monte Carlo noise.
# See:
#   - SimDesign package (Sigal & Chalmers, 2016) for structured simulations
#   - Göttmann, Frischkorn, Oberauer, Schaefer, & Schubert (2025) for a full
#     treatment of M3 parameter recovery methodology
