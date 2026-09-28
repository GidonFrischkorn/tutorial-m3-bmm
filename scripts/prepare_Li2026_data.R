# Data preparation for Tutorial 2 (also used in Tutorial 3): classifies the
# retrieval responses of Experiment 1 in Li, Frischkorn, & Oberauer (2026),
# JEP:LMC, https://doi.org/10.1037/xlm0001585 (data: https://osf.io/wpcx5/).
# Input:  data/Li_2026_ComplexSpan_Exp1.csv (trial level)
# Output: data/Li_2026_ComplexSpan_Exp1_agg.csv (person-by-condition counts)

###############################################################################!
# 0) R Setup -------------------------------------------------------------------
###############################################################################!

pacman::p_load(here, tidyverse)

###############################################################################!
# 1) Read Raw Data -------------------------------------------------------------
###############################################################################!

data_raw <- read_csv(
  here("data", "Li_2026_ComplexSpan_Exp1.csv"),
  show_col_types = FALSE
)

# Each trial shows 3 image pairs with a size judgment; one image per pair is
# cued as the memory target. Conditions: control (no distractors), pre-cue and
# retro-cue (target cued before or after the judgment). At retrieval,
# participants pick each target from 12 images.

# 18,000 rows: 50 participants x 60 trials (3 blocks of 20, one condition per
# block) x 6 screens (3 memory + 3 retrieval)

###############################################################################!
# 2) Participant Exclusions ----------------------------------------------------
###############################################################################!

# Exclude participants with a size judgment response rate < 50% or accuracy
# < 70% (among trials with a response), following the original analysis.
# The original analysis also used a post-experiment survey that is not part
# of the public OSF data, so only the performance criterion is applied here.

# First, identify trials with missing images (failed to load)
trials_with_na <- data_raw |>
  filter(screenID == "memory") |>
  group_by(participant, blockID, trialID) |>
  summarise(has_missing = any(is.na(leftImage) | is.na(rightImage)),
            .groups = "drop")

# Compute judgment accuracy per participant
# (excluding trials with missing images)
judge_performance <- data_raw |>
  filter(screenID == "memory") |>
  left_join(trials_with_na, by = c("participant", "blockID", "trialID")) |>
  filter(!has_missing) |>
  mutate(responded = response != "NULL" & !is.na(response)) |>
  group_by(participant) |>
  summarise(
    pResp = mean(responded),
    pAcc  = mean(acc[responded], na.rm = TRUE)
  )

excluded_ids <- judge_performance |>
  filter(pResp < 0.5 | pAcc < 0.7) |>
  pull(participant)

cat("Excluded", length(excluded_ids), "participants for low judgment accuracy:",
    paste(excluded_ids, collapse = ", "), "\n")
cat(
  "Remaining:",
  n_distinct(data_raw$participant) - length(excluded_ids),
  "participants\n\n"
)

###############################################################################!
# 3) Process Memory Screens ----------------------------------------------------
###############################################################################!

# Identify target and distractor for each memory screen.
# The cue indicates which image was the memory target.
data_memory <- data_raw |>
  filter(screenID == "memory",
         !participant %in% excluded_ids) |>
  left_join(trials_with_na, by = c("participant", "blockID", "trialID")) |>
  filter(!has_missing) |>
  mutate(
    target     = ifelse(cue == "left", leftImage, rightImage),
    distractor = ifelse(cue == "left", rightImage, leftImage)
  ) |>
  select(participant, blockID, trialID, condition, target, distractor)

# Build a trial-level lookup: all targets and all distractors per trial
trial_lookup <- data_memory |>
  group_by(participant, blockID, trialID, condition) |>
  summarise(
    all_targets     = list(target),
    all_distractors = list(distractor),
    .groups = "drop"
  )

# Build a paired-distractor lookup: which distractor was paired with each target
pair_lookup <- data_memory |>
  select(participant, blockID, trialID, target, distractor)

###############################################################################!
# 4) Classify Retrieval Responses ----------------------------------------------
###############################################################################!

data_retrieval <- data_raw |>
  filter(screenID == "retrieval",
         !participant %in% excluded_ids) |>
  select(participant, blockID, trialID, condition, correctObject, response)

# Remove retrieval tests from trials with missing images
data_retrieval <- data_retrieval |>
  semi_join(trial_lookup, by = c("participant", "blockID", "trialID"))

# Join with pair lookup to get the distractor paired with the tested target
data_retrieval <- data_retrieval |>
  left_join(pair_lookup,
            by = c("participant", "blockID", "trialID",
                   "correctObject" = "target"))

# Join with trial lookup to get all targets and distractors for the trial
data_retrieval <- data_retrieval |>
  left_join(trial_lookup |> select(-condition),
            by = c("participant", "blockID", "trialID"))

# Classify each response: corr (tested target), distc (distractor paired with
# the tested target), other (another target from the trial), disto (another
# distractor from the trial), npl (not-presented lure), or noRes (no response)
data_retrieval <- data_retrieval |>
  mutate(
    rcat = case_when(
      response == "NULL"                                   ~ "noRes",
      response == correctObject                            ~ "corr",
      response == distractor                               ~ "distc",
      map2_lgl(response, all_targets, ~ .x %in% .y)       ~ "other",
      map2_lgl(response, all_distractors, ~ .x %in% .y)   ~ "disto",
      TRUE                                                 ~ "npl"
    )
  )

# Verify: control condition should have no distractor responses
n_ctrl_dist <- data_retrieval |>
  filter(condition == "control", rcat %in% c("distc", "disto")) |>
  nrow()
stopifnot(n_ctrl_dist == 0)

# Report classification summary
cat("Response classification (all participants):\n")
data_retrieval |>
  count(condition, rcat) |>
  pivot_wider(names_from = rcat, values_from = n, values_fill = 0) |>
  print()
cat("\n")

###############################################################################!
# 5) Remove No-Response Trials and Aggregate -----------------------------------
###############################################################################!

# Drop trials where no response was given (response == "NULL")
data_retrieval <- data_retrieval |>
  filter(rcat != "noRes")

# Aggregate: count response frequencies per participant × condition
data_agg <- data_retrieval |>
  count(participant, condition, rcat) |>
  pivot_wider(names_from = rcat, values_from = n, values_fill = 0)

# Ensure all 5 response columns exist (even if all zeros)
for (col in c("corr", "other", "distc", "disto", "npl")) {
  if (!col %in% names(data_agg)) data_agg[[col]] <- 0L
}

###############################################################################!
# 6) Add Response Option Counts ------------------------------------------------
###############################################################################!

# Number of response options per category among the 12 retrieval images.
# Control: 1 tested target, 2 other targets, 9 lures (no distractors).
# Pre/retro: 1 tested target, 2 other targets, 1 paired distractor,
# 2 other distractors, 6 lures.
response_options <- tibble(
  condition = c("control", "pre", "retro"),
  n_corr    = c(1, 1, 1),
  n_other   = c(2, 2, 2),
  n_distc   = c(0, 1, 1),
  n_disto   = c(0, 2, 2),
  n_npl     = c(9, 6, 6)
)

data_agg <- data_agg |>
  left_join(response_options, by = "condition")

# Order condition factor: control, pre, retro
data_agg <- data_agg |>
  mutate(condition = factor(condition, levels = c("control", "pre", "retro")))

# Select and order columns for the final output
data_agg <- data_agg |>
  select(participant, condition,
         corr, other, distc, disto, npl,
         n_corr, n_other, n_distc, n_disto, n_npl)

###############################################################################!
# 7) Summary and Save ----------------------------------------------------------
###############################################################################!

cat("Final aggregated data:\n")
cat("  Dimensions:", nrow(data_agg), "rows ×", ncol(data_agg), "columns\n")
cat("  Participants:", n_distinct(data_agg$participant), "\n")
cat("  Conditions:", paste(levels(data_agg$condition), collapse = ", "), "\n\n")

cat("Response frequency summary:\n")
data_agg |>
  group_by(condition) |>
  summarise(across(corr:npl, ~ sprintf("%.1f (%.1f)", mean(.x), sd(.x)))) |>
  print()

# Save to CSV
write_csv(data_agg, here("data", "Li_2026_ComplexSpan_Exp1_agg.csv"))
cat("\nSaved to: data/Li_2026_ComplexSpan_Exp1_agg.csv\n")
