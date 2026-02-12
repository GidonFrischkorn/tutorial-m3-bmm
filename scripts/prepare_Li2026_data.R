#' Data Preparation: Li, Frischkorn, & Oberauer (accepted, JEP:LMC)
#'
#' This companion script prepares the trial-level data from Experiment 1
#' of Li, Frischkorn, & Oberauer for the M3 complex span model in bmm.
#' The prepared (aggregated) data is saved to data/Li_2026_ComplexSpan_Exp1_agg.csv
#' and used directly in Tutorial 2.
#'
#' Paper: Li, C., Frischkorn, G. T., & Oberauer, K. (accepted). Can We Process
#'        Information Without Encoding It into Working Memory? Journal of
#'        Experimental Psychology: Learning, Memory, and Cognition.
#' Data:  https://osf.io/wpcx5/overview
#'
#' Task overview:
#'   Novel complex-span task. On each trial, participants see 3 image pairs
#'   and judge which object in each pair is larger (size judgment task).
#'   One image per pair is cued as the memory target. At retrieval,
#'   participants identify each target from a set of 12 images.
#'
#' Conditions:
#'   - control: no distractors (non-target images are repeated fillers)
#'   - pre-cue: target/distractor identified BEFORE size judgment
#'   - retro-cue: target/distractor identified AFTER size judgment
#'
#' Response categories for the M3 complex span model:
#'   corr    = correct target selected
#'   other   = a different target from the same trial selected
#'   dist_c  = the distractor paired with the tested target selected
#'   dist_o  = a distractor from a different pair in the same trial selected
#'   npl     = a not-presented lure selected

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

# 50 participants, 3 conditions × 60 trials = 180 trials per participant
# Each trial: 3 memory screens + 3 retrieval screens = 6 rows per trial
# Total: 50 × 180 × 6 = 54,000... but actual is 18,000 rows because
# participants are assigned to ONE condition per block (3 blocks of 20 trials)

###############################################################################!
# 2) Participant Exclusions ----------------------------------------------------
###############################################################################!

# Exclude participants with poor performance on the size judgment task.
# Criteria (following the original analysis):
#   - response rate < 50%, OR
#   - accuracy < 70% (among trials with a response)
#
# Note: The published analysis also excluded participants based on a
# post-experiment survey (self-reported distraction or aid use). That survey
# data is not included in the public OSF dataset, so only the performance-
# based criterion is applied here.

# First, identify trials with missing images (failed to load)
trials_with_NA <- data_raw %>%
  filter(screenID == "memory") %>%
  group_by(participant, blockID, trialID) %>%
  summarise(has_missing = any(is.na(leftImage) | is.na(rightImage)),
            .groups = "drop")

# Compute judgment accuracy per participant (excluding trials with missing images)
judge_performance <- data_raw %>%
  filter(screenID == "memory") %>%
  left_join(trials_with_NA, by = c("participant", "blockID", "trialID")) %>%
  filter(!has_missing) %>%
  mutate(responded = response != "NULL" & !is.na(response)) %>%
  group_by(participant) %>%
  summarise(
    pResp = mean(responded),
    pAcc  = mean(acc[responded], na.rm = TRUE)
  )

excluded_ids <- judge_performance %>%
  filter(pResp < 0.5 | pAcc < 0.7) %>%
  pull(participant)

cat("Excluded", length(excluded_ids), "participants for low judgment accuracy:",
    paste(excluded_ids, collapse = ", "), "\n")
cat("Remaining:", n_distinct(data_raw$participant) - length(excluded_ids), "participants\n\n")

###############################################################################!
# 3) Process Memory Screens ----------------------------------------------------
###############################################################################!

# Identify target and distractor for each memory screen.
# The cue indicates which image was the memory target.
data_memory <- data_raw %>%
  filter(screenID == "memory",
         !participant %in% excluded_ids) %>%
  left_join(trials_with_NA, by = c("participant", "blockID", "trialID")) %>%
  filter(!has_missing) %>%
  mutate(
    target     = ifelse(cue == "left", leftImage, rightImage),
    distractor = ifelse(cue == "left", rightImage, leftImage)
  ) %>%
  select(participant, blockID, trialID, condition, target, distractor)

# Build a trial-level lookup: all targets and all distractors per trial
trial_lookup <- data_memory %>%
  group_by(participant, blockID, trialID, condition) %>%
  summarise(
    all_targets     = list(target),
    all_distractors = list(distractor),
    .groups = "drop"
  )

# Build a paired-distractor lookup: which distractor was paired with each target
pair_lookup <- data_memory %>%
  select(participant, blockID, trialID, target, distractor)

###############################################################################!
# 4) Classify Retrieval Responses ----------------------------------------------
###############################################################################!

data_retrieval <- data_raw %>%
  filter(screenID == "retrieval",
         !participant %in% excluded_ids) %>%
  select(participant, blockID, trialID, condition, correctObject, response)

# Remove retrieval tests from trials with missing images
data_retrieval <- data_retrieval %>%
  semi_join(trial_lookup, by = c("participant", "blockID", "trialID"))

# Join with pair lookup to get the distractor paired with the tested target
data_retrieval <- data_retrieval %>%
  left_join(pair_lookup,
            by = c("participant", "blockID", "trialID",
                    "correctObject" = "target"))

# Join with trial lookup to get all targets and distractors for the trial
data_retrieval <- data_retrieval %>%
  left_join(trial_lookup %>% select(-condition),
            by = c("participant", "blockID", "trialID"))

# Classify each response into one of 5 categories (+ no response)
data_retrieval <- data_retrieval %>%
  mutate(
    rcat = case_when(
      response == "NULL"                                   ~ "noRes",
      response == correctObject                            ~ "corr",
      response == distractor                               ~ "dist_c",
      map2_lgl(response, all_targets, ~ .x %in% .y)       ~ "other",
      map2_lgl(response, all_distractors, ~ .x %in% .y)   ~ "dist_o",
      TRUE                                                 ~ "npl"
    )
  )

# Verify: control condition should have no distractor responses
n_ctrl_dist <- data_retrieval %>%
  filter(condition == "control", rcat %in% c("dist_c", "dist_o")) %>%
  nrow()
stopifnot(n_ctrl_dist == 0)

# Report classification summary
cat("Response classification (all participants):\n")
data_retrieval %>%
  count(condition, rcat) %>%
  pivot_wider(names_from = rcat, values_from = n, values_fill = 0) %>%
  print()
cat("\n")

###############################################################################!
# 5) Remove No-Response Trials and Aggregate -----------------------------------
###############################################################################!

# Drop trials where no response was given (response == "NULL")
data_retrieval <- data_retrieval %>%
  filter(rcat != "noRes")

# Aggregate: count response frequencies per participant × condition
data_agg <- data_retrieval %>%
  count(participant, condition, rcat) %>%
  pivot_wider(names_from = rcat, values_from = n, values_fill = 0)

# Ensure all 5 response columns exist (even if all zeros)
for (col in c("corr", "other", "dist_c", "dist_o", "npl")) {
  if (!col %in% names(data_agg)) data_agg[[col]] <- 0L
}

###############################################################################!
# 6) Add Response Option Counts ------------------------------------------------
###############################################################################!

# The number of response options per category depends on the condition.
# At retrieval, participants choose from 12 images:
#   - control: 3 targets + 9 NPL (no distractors)
#   - pre/retro: 3 targets + 3 distractors + 6 NPL
#
# For the tested item:
#   n_corr   = 1 (the tested target)
#   n_other  = 2 (the other 2 targets from the same trial)
#   n_dist_c = 0 (control) or 1 (pre/retro: distractor paired with tested target)
#   n_dist_o = 0 (control) or 2 (pre/retro: distractors from other pairs)
#   n_npl    = 9 (control) or 6 (pre/retro: remaining images)
response_options <- tibble(
  condition = c("control", "pre", "retro"),
  n_corr   = c(1, 1, 1),
  n_other  = c(2, 2, 2),
  n_dist_c = c(0, 1, 1),
  n_dist_o = c(0, 2, 2),
  n_npl    = c(9, 6, 6)
)

data_agg <- data_agg %>%
  left_join(response_options, by = "condition")

# Order condition factor: control, pre, retro
data_agg <- data_agg %>%
  mutate(condition = factor(condition, levels = c("control", "pre", "retro")))

# Select and order columns for the final output
data_agg <- data_agg %>%
  select(participant, condition,
         corr, other, dist_c, dist_o, npl,
         n_corr, n_other, n_dist_c, n_dist_o, n_npl)

###############################################################################!
# 7) Summary and Save ----------------------------------------------------------
###############################################################################!

cat("Final aggregated data:\n")
cat("  Dimensions:", nrow(data_agg), "rows ×", ncol(data_agg), "columns\n")
cat("  Participants:", n_distinct(data_agg$participant), "\n")
cat("  Conditions:", paste(levels(data_agg$condition), collapse = ", "), "\n\n")

cat("Response frequency summary:\n")
data_agg %>%
  group_by(condition) %>%
  summarise(across(corr:npl, ~ sprintf("%.1f (%.1f)", mean(.x), sd(.x)))) %>%
  print()

# Save to CSV
write_csv(data_agg, here("data", "Li_2026_ComplexSpan_Exp1_agg.csv"))
cat("\nSaved to: data/Li_2026_ComplexSpan_Exp1_agg.csv\n")
