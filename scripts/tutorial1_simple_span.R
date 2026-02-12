#' Tutorial 1: Simple Span M3 Model with bmm
#'
#' This script demonstrates the complete M3 workflow using data from:
#' Oberauer, K. (2019). Working Memory Capacity Limits Memory for Bindings.
#' Journal of Cognition, 2(1), 40. https://doi.org/10/gf8p8c
#' Data available at: https://osf.io/qy5sd/
#'
#' In this script, you will see:
#'  1) How to prepare trial-level data for the M3 simple span model
#'  2) How to specify and fit the M3 model with bmm (simple and softmax rules)
#'  3) How to evaluate model fit with posterior predictive checks
#'  4) How to extract, interpret, and plot parameter estimates

###############################################################################!
# 0) R Setup: Packages & Data -------------------------------------------------
###############################################################################!

pacman::p_load(here, bmm, brms, tidyverse, tidybayes, patchwork, gghalves)
source(here("functions", "clean_plot.R"))

# -- Fitting settings ----------------------------------------------------------
chains <- 4
cores  <- 4
warmup <- 1000
iter   <- 2000

## 0.1) Read trial-level data --------------------------------------------------

# The raw data files are space-delimited with no header row.
# Column names are taken from the original JAGS analysis script
# (Oberauer_2019_MMM_Analyses.R):
#   id        - participant identifier
#   session   - session number
#   block     - block number
#   trial     - trial number within block
#   setsize   - number of items in the memory set (2, 4, 6, or 8)
#   rsizeList - recognition set: number of list items presented (includes the
#               correct item, so the number of *other* list lures is
#               rsizeList - 1)
#   rsizeNPL  - recognition set: number of not-presented lures
#   tested    - serial position of the tested item
#   response  - participant's response
#   rcat      - response category (1 = correct, 2 = other list item, 3 = NPL)
#   rt        - response time (log-transformed; -1 for recall trials)
col_names <- c("id", "session", "block", "trial", "setsize",
               "rsizeList", "rsizeNPL", "tested", "response", "rcat", "rt")

data_raw <- read.table(
  here("data", "Oberauer_2019_SimpleSpan_Exp1.dat"),
  header = FALSE,
  col.names = col_names,
  fill = TRUE       # handle trailing whitespace at end of file
) %>%
  drop_na()          # remove any incomplete rows

# Quick overview of the raw data
str(data_raw)

## 0.2) Filter to recognition trials -------------------------------------------

# The dataset contains both recognition and recall trials. Recognition trials
# are identified by having at least one response option (rsizeList + rsizeNPL > 0).
# Recall trials have rsizeList = 0 and rsizeNPL = 0.
data_recog <- data_raw %>%
  filter(rsizeList + rsizeNPL > 0)

# Verify: recognition trials should have valid response categories (1, 2, or 3)
stopifnot(all(data_recog$rcat %in% c(1, 2, 3)))

## 0.3) Verify response set constraints ----------------------------------------

# Sanity check: the number of list items in the recognition set should never
# exceed the set size
stopifnot(all(data_recog$rsizeList <= data_recog$setsize))

# Check the unique recognition conditions (setsize × rsizeList × rsizeNPL)
recog_conditions <- data_recog %>%
  distinct(setsize, rsizeList, rsizeNPL) %>%
  arrange(setsize, rsizeList, rsizeNPL)

print(recog_conditions)
# 24 unique conditions: set size determines how many recognition conditions
# are possible (more items = more possible recognition set compositions)

## 0.4) Aggregate to response frequencies --------------------------------------

# The M3 model in bmm expects aggregated data: response frequencies per
# participant and condition. We count how often each response category
# occurred for each participant × (setsize, rsizeList, rsizeNPL) combination.
data_agg <- data_recog %>%
  group_by(id, setsize, rsizeList, rsizeNPL) %>%
  summarise(
    corr  = sum(rcat == 1),   # correct responses
    other = sum(rcat == 2),   # other list item intrusions
    npl   = sum(rcat == 3),   # not-presented lure responses
    .groups = "drop"
  )

## 0.5) Compute response option counts -----------------------------------------

# bmm needs to know how many response options exist in each category.
# For the simple span model (version = "ss"):
#   n_corr  = 1 (there is always exactly one correct item)
#   n_other = rsizeList - 1 (other list items in the recognition set;
#             rsizeList includes the correct item, so we subtract 1)
#   n_npl   = rsizeNPL (not-presented lures in the recognition set)
data_agg <- data_agg %>%
  mutate(
    n_corr  = 1,
    n_other = rsizeList - 1,
    n_npl   = rsizeNPL
  )

# Verify: no negative response option counts
stopifnot(all(data_agg$n_other >= 0))
stopifnot(all(data_agg$n_npl >= 0))

# Verify: when n_other = 0, there should be no "other" responses
stopifnot(all(data_agg$other[data_agg$n_other == 0] == 0))

## 0.6) Prepare condition factor ------------------------------------------------

# Set size is the main experimental manipulation. Convert to a factor for
# categorical effects in the model formula (separate estimates per set size).
data_agg <- data_agg %>%
  mutate(setsize = factor(setsize))

# Final data overview
str(data_agg)
head(data_agg)

# Summary: 20 participants × 24 conditions = 480 rows
# Each row contains response frequencies and response option counts
# for one participant in one recognition condition
cat("\nData dimensions:", nrow(data_agg), "rows ×", ncol(data_agg), "columns\n")
cat("Participants:", n_distinct(data_agg$id), "\n")
cat("Conditions per participant:", nrow(data_agg) / n_distinct(data_agg$id), "\n")

###############################################################################!
# 1) Model Specification -------------------------------------------------------
###############################################################################!

# TODO: specify m3 model with m3() and bmf()

###############################################################################!
# 2) Model Fitting -------------------------------------------------------------
###############################################################################!

# TODO: fit model with bmm()

###############################################################################!
# 3) Model Evaluation ----------------------------------------------------------
###############################################################################!

# TODO: posterior predictive checks, parameter estimates, hypothesis testing
