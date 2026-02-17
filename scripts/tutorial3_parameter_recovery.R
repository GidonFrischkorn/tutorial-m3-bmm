#' Tutorial 3: Custom M3 — Parameter Recovery
#'
#' This script demonstrates how to define a fully custom M3 model using bmm
#' and validates it through parameter recovery simulation. The custom model is
#' the memory updating model from Oberauer & Lewandowsky (2019), with 5
#' response categories and 5 estimated parameters.
#'
#' Instead of fitting real data, this tutorial:
#'   1) Specifies a custom M3 model with user-defined activation formulas
#'   2) Simulates data from the model with known generating parameters
#'   3) Fits the model to recover the generating parameters
#'   4) Evaluates recovery quality at the group and individual level
#'
#' The simulation varies sample size (N = 25, 50, 100) and trials per
#' condition (5, 10, 25) to illustrate how data quantity affects recovery.
#'
#' Reference:
#'   Oberauer, K., & Lewandowsky, S. (2019). Simple measurement models for
#'   complex working memory tasks. Psychological Review, 126(6), 880–932.

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
set.seed(2025)

###############################################################################!
# 1) Model Specification -------------------------------------------------------
###############################################################################!

# The memory updating model uses 5 response categories and 5 estimated
# parameters. The task presents set size 5 items, then updates each position
# once. At test, participants identify the current target from a recognition
# set containing current items, replaced items, and not-presented lures.
#
# Activation formulas (softmax choice rule, b fixed at 0):
#
#   correct  = b + (1 + e*te) * (a + c)
#   other    = b + (1 + e*te) * a
#   oldinpos = b + exp(-r*tr) * d * (1 + e*te) * (a + c)
#   otherold = b + exp(-r*tr) * d * (1 + e*te) * a
#   npl      = b
#
# Where:
#   a = general item memory (identity link)
#   c = context-specific binding (identity link)
#   d = residual activation of replaced items, 0–1 (logit link)
#   e = extended encoding strength per unit time (log link)
#   r = removal rate of outdated items (log link)
#   b = baseline activation (fixed at 0 for softmax)
#   te = time available for extended encoding (data variable, seconds)
#   tr = time available for removal (data variable, seconds)
#
# Cognitive interpretation:
#   (1 + e*te) scales activation by encoding opportunity — longer encoding
#   time boosts all non-baseline activation proportionally via e.
#   exp(-r*tr) * d captures residual activation of replaced items — d is the
#   base residual strength and exp(-r*tr) attenuates it over removal time.

## 1.1) Define the custom M3 model object --------------------------------------

m3_model_custom <- m3(
  resp_cats   = c("correct", "other", "oldinpos", "otherold", "npl"),
  num_options = c("n_correct" = 1, "n_other" = 4, "n_oldinpos" = 1, "n_otherold" = 4, "n_npl" = 5),
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

# These formulas define how total activation accumulates for each response
# category. bmm parses exp() and arithmetic expressions internally,
# converting them to non-linear brms formulas.
#
# Note: in bmf(), `*` is arithmetic multiplication (not formula interaction),
# because these are non-linear model expressions evaluated as R code.
m3_act_funs <- bmf(
  correct  ~ (1 + e * te) * (a + c) + b,
  other    ~ (1 + e * te) * a + b,
  oldinpos ~ exp(-r * tr) * d * (1 + e * te) * (a + c) + b,
  otherold ~ exp(-r * tr) * d * (1 + e * te) * a + b,
  npl      ~ b
)

## 1.3) Define the fitting formula --------------------------------------------

# The fitting formula combines activation formulas (Section 1.2) with
# parameter prediction formulas. Each parameter gets an intercept (group mean)
# and a participant-level random intercept for between-person variability.
# No experimental conditions are predicted — te and tr enter through the
# activation formulas directly.
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

## 2.1) Generating parameter values -------------------------------------------

# Group-level means and between-person SDs on the link scale.
# The link function determines the scale on which parameters are estimated:
#   identity: link = native (no transformation)
#   logit:    native = plogis(link), constrains to (0, 1)
#   log:      native = exp(link), constrains to positive values
#
# Between-person SDs are set so that ±2 SD on the link scale spans a
# cognitively plausible range on the native scale.
gen_pars <- tibble(
  parameter = c("c",    "a",    "d",     "e",     "r"),
  link      = c("identity", "identity", "logit", "log",   "log"),
  mean_link = c( 2.0,   1.0,    0.85,    0.00,   -0.69),
  sd_link   = c( 0.50,  0.50,   0.65,    0.35,    0.38)
)

gen_pars

## 2.2) Draw true individual parameters for N_max = 100 -----------------------

# We generate parameters for the maximum sample size (N = 100) once.
# Smaller sample sizes use the first N participants from this set, ensuring
# that effects of sample size are not confounded with different true values.
N_max <- 100

true_link <- tibble(participant = 1:N_max)
for (i in seq_len(nrow(gen_pars))) {
  p <- gen_pars$parameter[i]
  true_link[[p]] <- rnorm(N_max, gen_pars$mean_link[i], gen_pars$sd_link[i])
}

# Convert to native scale for use in rm3()
true_native <- true_link %>%
  mutate(
    c_native = c,                    # identity
    a_native = a,                    # identity
    d_native = plogis(d),            # logit → probability
    e_native = exp(e),               # log → positive
    r_native = exp(r)                # log → positive
  )

# Verify generating distributions look reasonable
true_native %>%
  select(participant, ends_with("_native")) %>%
  pivot_longer(-participant, names_to = "parameter", values_to = "value") %>%
  group_by(parameter) %>%
  summarise(mean = mean(value), sd = sd(value),
            min = min(value), max = max(value)) %>%
  print()

## 2.3) Experimental design grid ----------------------------------------------

# 3 × 3 factorial: te (encoding time) × tr (removal time)
# Each participant completes all 9 conditions.
design_grid <- expand_grid(
  te = c(0.25, 0.5, 1.0),
  tr = c(0.25, 0.5, 1.0)
)

# Fixed response option counts (same for all conditions)
design_grid <- design_grid %>%
  mutate(
    n_correct  = 1L,
    n_other    = 4L,
    n_oldinpos = 1L,
    n_otherold = 4L,
    n_npl      = 5L
  )

design_grid

## 2.4) Simulation function ----------------------------------------------------

#' Simulate M3 data for one cell of the recovery design
#'
#' @param N Number of participants (uses the first N from true_native)
#' @param trials_per_cond Number of trials per te × tr condition
#' @param true_native Tibble with true native-scale parameters per participant
#' @param design_grid Tibble with te, tr, and n_options columns
#' @param m3_model The custom m3 model object
#' @param act_funs The activation formulas (bmf object)
#' @return Tibble with simulated response frequencies (aggregated)
simulate_cell <- function(N, trials_per_cond, true_native,
                          design_grid, m3_model, act_funs) {
  n_conds <- nrow(design_grid)
  data_list <- vector("list", N)

  iter <- 1
  for (i in seq_len(N)) {
    pars_i <- c(
      a = true_native$a_native[i],
      c = true_native$c_native[i],
      d = true_native$d_native[i],
      e = true_native$e_native[i],
      r = true_native$r_native[i]
    )

    for(j in seq_len(n_conds)) {
      cat(sprintf("Simulating participant %d, condition %d/%d\n",
                  i, j, n_conds))

      sim_ij <- rm3(
        n        = 1,
        size     = trials_per_cond,
        pars     = pars_i,
        m3_model = m3_model,
        act_funs = act_funs,
        te       = design_grid$te[j],
        tr       = design_grid$tr[j]
      )

      sim_ij <- data.frame(sim_ij)

      sim_ij$participant <- i
      sim_ij$te <- design_grid$te[j]
      sim_ij$tr <- design_grid$tr[j]
      sim_ij$n_correct  <- design_grid$n_correct[j]
      sim_ij$n_other    <- design_grid$n_other[j]
      sim_ij$n_oldinpos <- design_grid$n_oldinpos[j]
      sim_ij$n_otherold <- design_grid$n_otherold[j]
      sim_ij$n_npl      <- design_grid$n_npl[j]
      data_list[[iter]] <- sim_ij
      iter <- iter + 1
    }
  }

  bind_rows(data_list)
}

## 2.5) Simulate data for all 9 cells -----------------------------------------

# Recovery design: 3 sample sizes × 3 trial counts = 9 cells
sample_sizes  <- c(25, 50, 100)
trial_counts  <- c(5, 10, 25)

recovery_grid <- expand_grid(
  N               = sample_sizes,
  trials_per_cond = trial_counts
)

# Simulate data for each cell. Each cell uses the first N participants
# from the same set of true parameters, ensuring comparability.
sim_data <- recovery_grid %>%
  mutate(
    cell_id = row_number(),
    data = map2(N, trials_per_cond, ~ {
      cat(sprintf("  Simulating: N = %d, trials/cond = %d\n", .x, .y))
      simulate_cell(.x, .y, true_native, design_grid,
                    m3_model_custom, m3_act_funs)
    })
  )

## 2.6) Verify simulated data -------------------------------------------------

# Quick check: inspect one cell (N = 50, trials = 25)
example_data <- sim_data %>%
  filter(N == 50, trials_per_cond == 25) %>%
  pull(data) %>%
  .[[1]]

# Should have 50 participants × 9 conditions = 450 rows
cat("Example cell dimensions:", nrow(example_data), "rows\n")
str(example_data)

# Response proportions should show correct > other > oldinpos/otherold > npl
example_data %>%
  summarise(across(correct:npl, ~ mean(.x / trials_per_cond))) %>%
  print()

###############################################################################!
# 3) Model Fitting -------------------------------------------------------------
###############################################################################!

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

## 3.1) Check default priors ---------------------------------------------------

# Inspect the default priors bmm assigns to the custom model.
# These may need adjustment if they are too wide or misspecified.
default_prior(m3_formula_fit, m3_model_custom_fit, data = example_data)

## 3.2) Fit all 9 cells -------------------------------------------------------

# Each cell is fitted independently. The file argument caches the fit so
# re-running the script skips completed fits. With N = 100 and 9 conditions,
# the largest cells may take substantial time.
#
# Fitting progress is tracked by cell_id. On a multi-core machine, cells
# could be parallelized across separate R sessions.
fit_results <- sim_data %>%
  mutate(
    fit = pmap(list(data, N, trials_per_cond, cell_id), function(d, n, t, id) {
      cat("Fitting cell ", id, ": N = ", n, ", trials/cond = ", t, "\n", sep = "")
      bmm(
        formula = m3_formula_fit,
        model   = m3_model_custom_fit,
        data    = d,
        chains  = chains,
        cores   = cores,
        warmup  = warmup,
        iter    = iter,
        backend = "cmdstanr",
        file    = here("output", paste0("fit_m3_custom_N", n, "_T", t))
      )
    })
  )

## 3.3) Convergence checks -----------------------------------------------------

# Check Rhat and ESS for each cell. All Rhat should be < 1.01.
convergence_summary <- fit_results %>%
  mutate(
    rhat_max = map_dbl(fit, ~ max(rhat(.x), na.rm = TRUE)),
    ess_min  = map_dbl(fit, ~ min(
      c(neff_ratio(.x) * ((iter - warmup) * chains)), na.rm = TRUE
    ))
  ) %>%
  select(cell_id, N, trials_per_cond, rhat_max, ess_min)

print(convergence_summary)

###############################################################################!
# 4) Recovery Evaluation -------------------------------------------------------
###############################################################################!

## 4.1) Extract group-level estimates ------------------------------------------

# For each cell, extract the fixed effects (group-level means on link scale)
# and compare to the true generating means.
group_recovery <- fit_results %>%
  mutate(
    fixefs = map(fit, ~ {
      fe <- fixef(.x)
      tibble(
        parameter = c("a", "c", "d", "e", "r"),
        recovered_mean  = fe[paste0(c("a", "c", "d", "e", "r"),
                                     "_Intercept"), "Estimate"],
        recovered_lower = fe[paste0(c("a", "c", "d", "e", "r"),
                                     "_Intercept"), "Q2.5"],
        recovered_upper = fe[paste0(c("a", "c", "d", "e", "r"),
                                     "_Intercept"), "Q97.5"]
      )
    })
  ) %>%
  select(cell_id, N, trials_per_cond, fixefs) %>%
  unnest(fixefs)

# Add true generating means (on link scale)
group_recovery <- group_recovery %>%
  left_join(gen_pars %>% select(parameter, true_mean = mean_link),
            by = "parameter")

# Compute recovery metrics
group_recovery <- group_recovery %>%
  mutate(
    bias     = recovered_mean - true_mean,
    coverage = true_mean >= recovered_lower & true_mean <= recovered_upper
  )

## 4.2) Group-level recovery plot ----------------------------------------------

group_plot <- ggplot(group_recovery,
                     aes(x = true_mean, y = recovered_mean,
                         color = factor(N))) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  geom_pointrange(aes(ymin = recovered_lower, ymax = recovered_upper),
                  position = position_dodge(width = 0.15), size = 0.4) +
  facet_grid(trials_per_cond ~ parameter, scales = "free",
             labeller = labeller(
               trials_per_cond = ~ paste(.x, "trials/cond")
             )) +
  scale_color_m3() +
  labs(x = "True value (link scale)",
       y = "Recovered value (link scale)",
       color = "Sample size") +
  clean_plot()

group_plot

ggsave(here("figures", "tutorial3_group_recovery.pdf"),
       group_plot, width = 6.5, height = 6)

# Group-level recovery summary table
group_summary <- group_recovery %>%
  group_by(N, trials_per_cond) %>%
  summarise(
    mean_bias     = mean(bias),
    mean_abs_bias = mean(abs(bias)),
    coverage      = mean(coverage),
    .groups = "drop"
  )

print(group_summary)

## 4.3) Extract group-level SD estimates ---------------------------------------

# Compare recovered random effect SDs to the true generating SDs.
sd_recovery <- fit_results %>%
  mutate(
    sd_ests = map(fit, ~ {
      vc <- VarCorr(.x)
      # Extract SD estimates for each parameter from the participant group
      par_names <- c("a", "c", "d", "e", "r")
      map_dfr(par_names, function(p) {
        sd_est <- vc$participant$sd[paste0(p, "_Intercept"), ]
        tibble(
          parameter      = p,
          recovered_sd   = sd_est["Estimate"],
          sd_lower       = sd_est["Q2.5"],
          sd_upper       = sd_est["Q97.5"]
        )
      })
    })
  ) %>%
  select(cell_id, N, trials_per_cond, sd_ests) %>%
  unnest(sd_ests) %>%
  left_join(gen_pars %>% select(parameter, true_sd = sd_link),
            by = "parameter") %>%
  mutate(
    sd_bias    = recovered_sd - true_sd,
    sd_coverage = true_sd >= sd_lower & true_sd <= sd_upper
  )

sd_summary <- sd_recovery %>%
  group_by(N, trials_per_cond) %>%
  summarise(
    mean_sd_bias  = mean(sd_bias),
    sd_coverage   = mean(sd_coverage),
    .groups = "drop"
  )

print(sd_summary)

## 4.4) Individual-level recovery ----------------------------------------------

# For each cell, extract individual-level estimates (fixed + random effect)
# and compare to the true generating values on the link scale.
#
# coef() returns the total estimate per participant: fixed effect + random
# effect deviation. These are on the link scale.
indiv_recovery <- fit_results %>%
  mutate(
    indiv = map2(fit, N, function(fit_obj, n) {
      par_names <- c("a", "c", "d", "e", "r")
      co <- coef(fit_obj)

      map_dfr(par_names, function(p) {
        # coef() returns a list with one element per group variable
        est <- co$participant[, , paste0(p, "_Intercept")]
        tibble(
          participant    = 1:n,
          parameter      = p,
          recovered_link = est[1:n, "Estimate"],
          indiv_lower    = est[1:n, "Q2.5"],
          indiv_upper    = est[1:n, "Q97.5"]
        )
      })
    })
  ) %>%
  select(cell_id, N, trials_per_cond, indiv) %>%
  unnest(indiv)

# Add true individual link-scale values
true_long <- true_native %>%
  select(participant, a = a, c = c, d = d, e = e, r = r) %>%
  pivot_longer(-participant, names_to = "parameter",
               values_to = "true_link")

indiv_recovery <- indiv_recovery %>%
  left_join(true_long, by = c("participant", "parameter"))

# Compute individual-level recovery metrics per cell × parameter
indiv_summary <- indiv_recovery %>%
  group_by(cell_id, N, trials_per_cond, parameter) %>%
  summarise(
    correlation = cor(true_link, recovered_link),
    mean_bias   = mean(recovered_link - true_link),
    rmse        = sqrt(mean((recovered_link - true_link)^2)),
    coverage    = mean(true_link >= indiv_lower & true_link <= indiv_upper),
    .groups     = "drop"
  )

print(indiv_summary, n = 45)

## 4.5) Individual-level scatter plots -----------------------------------------

# Show true vs. recovered for each parameter in the N=100, trials=25 cell
indiv_best <- indiv_recovery %>%
  filter(N == 100, trials_per_cond == 25)

# Correlation labels for each panel
corr_labels <- indiv_best %>%
  group_by(parameter) %>%
  summarise(correlation = cor(true_link, recovered_link), .groups = "drop") %>%
  mutate(label = paste0("r = ", sprintf("%.2f", correlation)))

scatter_plot <- ggplot(indiv_best,
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

ggsave(here("figures", "tutorial3_indiv_recovery_scatter.pdf"),
       scatter_plot, width = 6.5, height = 5)

## 4.6) Recovery across cells: correlation heatmap -----------------------------

# How does individual-level correlation change with N and trials per parameter?
corr_plot <- ggplot(indiv_summary,
                    aes(x = factor(trials_per_cond), y = factor(N),
                        fill = correlation)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", correlation)), size = 2.8) +
  facet_wrap(~ parameter, nrow = 1) +
  scale_fill_gradient(low = "white", high = "#0072B2",
                      limits = c(0, 1), name = "r") +
  labs(x = "Trials per condition",
       y = "Sample size (N)") +
  clean_plot()

corr_plot

ggsave(here("figures", "tutorial3_recovery_heatmap.pdf"),
       corr_plot, width = 6.5, height = 3.5)

## 4.7) Note on simulation design ---------------------------------------------

# This tutorial demonstrates the parameter recovery workflow with a single
# replication per cell. Because each cell is fitted only once, the results
# are subject to Monte Carlo noise — a different random seed would produce
# somewhat different recovery statistics.
#
# For publication-quality simulation studies, each cell should be replicated
# many times (e.g., 100–500 replications) and results averaged across
# replications. This can be implemented using:
#   - SimDesign package (Sigal & Chalmers, 2016) for structured simulation
#   - furrr package for parallel execution across replications
#   - For a full treatment of M3 parameter recovery methodology, see:
#     Göttmann, Frischkorn, Oberauer, Schaefer, & Schubert (2025)
