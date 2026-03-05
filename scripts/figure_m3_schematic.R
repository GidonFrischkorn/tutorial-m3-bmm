#' Generate the M3 framework schematic figure
#'
#' This script creates a diagram illustrating the M3 architecture:
#' activation sources -> competitive selection -> response probabilities.
#' Uses the simple span example with three response categories.

pacman::p_load(here, ggplot2, grid, gridExtra)
source(here("functions", "clean_plot.R"))

# Build the figure using grid graphics
create_m3_schematic <- function() {

  # --- Layout parameters ---
  # Three columns: activation sources | choice rule | response probabilities
  # Three rows: correct, other, NPL

  pdf(here("figures", "m3_framework_schematic.pdf"), width = 6.5, height = 4)

  grid.newpage()

  # Use a viewport with some margin
  pushViewport(viewport(x = 0.5, y = 0.5, width = 0.95, height = 0.90))

  # --- Column positions ---
  col_params  <- 0.08   # activation source labels (left)
  col_formula <- 0.38   # activation formula column
  col_rule    <- 0.62   # choice rule
  col_prob    <- 0.85   # response probability

  # --- Row positions (top to bottom) ---
  row_title <- 0.95
  row1 <- 0.75  # correct
  row2 <- 0.50  # other
  row3 <- 0.25  # NPL

  # --- Title labels ---
  grid.text("Activation Sources", x = col_params, y = row_title,
            gp = gpar(fontsize = 10, fontface = "bold"))
  grid.text("Activation Formulas", x = col_formula, y = row_title,
            gp = gpar(fontsize = 10, fontface = "bold"))
  grid.text("Choice Rule", x = col_rule, y = row_title,
            gp = gpar(fontsize = 10, fontface = "bold"))
  grid.text("Response\nProbabilities", x = col_prob, y = row_title,
            gp = gpar(fontsize = 10, fontface = "bold"))

  # --- Parameter boxes (left side) ---
  box_w <- 0.10
  box_h <- 0.08
  param_col <- 0.08

  # b parameter (spans all rows)
  grid.roundrect(x = param_col, y = row3, width = box_w, height = box_h,
                 gp = gpar(fill = "#E8E8E8", col = "grey40"))
  grid.text("b", x = param_col, y = row3,
            gp = gpar(fontsize = 11, fontface = "italic"))

  # a parameter
  grid.roundrect(x = param_col, y = row2, width = box_w, height = box_h,
                 gp = gpar(fill = "#B3CDE3", col = "grey40"))
  grid.text("a", x = param_col, y = row2,
            gp = gpar(fontsize = 11, fontface = "italic"))

  # c parameter
  grid.roundrect(x = param_col + 0.12, y = row1, width = box_w, height = box_h,
                 gp = gpar(fill = "#FBB4AE", col = "grey40"))
  grid.text("c", x = param_col + 0.12, y = row1,
            gp = gpar(fontsize = 11, fontface = "italic"))

  # --- Activation formula boxes ---
  formula_w <- 0.18
  formula_h <- 0.10

  # Correct: b + a + c
  grid.roundrect(x = col_formula, y = row1, width = formula_w, height = formula_h,
                 gp = gpar(fill = "white", col = "grey40"))
  grid.text(expression(A[correct] == b + a + c),
            x = col_formula, y = row1, gp = gpar(fontsize = 9))

  # Other: b + a
  grid.roundrect(x = col_formula, y = row2, width = formula_w, height = formula_h,
                 gp = gpar(fill = "white", col = "grey40"))
  grid.text(expression(A[other] == b + a),
            x = col_formula, y = row2, gp = gpar(fontsize = 9))

  # NPL: b
  grid.roundrect(x = col_formula, y = row3, width = formula_w, height = formula_h,
                 gp = gpar(fill = "white", col = "grey40"))
  grid.text(expression(A[NPL] == b),
            x = col_formula, y = row3, gp = gpar(fontsize = 9))

  # --- Arrows: params -> formulas ---
  arrow_gp <- gpar(col = "grey40", lwd = 1)
  arr <- arrow(length = unit(0.08, "inches"), type = "closed")

  # b -> all formulas
  grid.lines(x = c(param_col + box_w/2, col_formula - formula_w/2),
             y = c(row3, row3), gp = arrow_gp, arrow = arr)
  grid.lines(x = c(param_col + box_w/2, col_formula - formula_w/2),
             y = c(row3, row2), gp = arrow_gp, arrow = arr)
  grid.lines(x = c(param_col + box_w/2, col_formula - formula_w/2),
             y = c(row3, row1), gp = arrow_gp, arrow = arr)

  # a -> correct and other
  grid.lines(x = c(param_col + box_w/2, col_formula - formula_w/2),
             y = c(row2, row2), gp = arrow_gp, arrow = arr)
  grid.lines(x = c(param_col + box_w/2, col_formula - formula_w/2),
             y = c(row2, row1), gp = arrow_gp, arrow = arr)

  # c -> correct only
  grid.lines(x = c(param_col + 0.12 + box_w/2, col_formula - formula_w/2),
             y = c(row1, row1), gp = arrow_gp, arrow = arr)

  # --- Choice rule box ---
  rule_w <- 0.14
  rule_h <- 0.30
  grid.roundrect(x = col_rule, y = row2, width = rule_w, height = rule_h,
                 gp = gpar(fill = "#FFFFCC", col = "grey40"))
  grid.text("Softmax\nor\nSimple", x = col_rule, y = row2,
            gp = gpar(fontsize = 9))

  # --- Arrows: formulas -> choice rule ---
  grid.lines(x = c(col_formula + formula_w/2, col_rule - rule_w/2),
             y = c(row1, row1 - 0.05), gp = arrow_gp, arrow = arr)
  grid.lines(x = c(col_formula + formula_w/2, col_rule - rule_w/2),
             y = c(row2, row2), gp = arrow_gp, arrow = arr)
  grid.lines(x = c(col_formula + formula_w/2, col_rule - rule_w/2),
             y = c(row3, row3 + 0.05), gp = arrow_gp, arrow = arr)

  # --- Response probability labels ---
  prob_w <- 0.16
  prob_h <- 0.10

  grid.roundrect(x = col_prob, y = row1, width = prob_w, height = prob_h,
                 gp = gpar(fill = "white", col = "grey40"))
  grid.text(expression(P(correct)), x = col_prob, y = row1,
            gp = gpar(fontsize = 9))

  grid.roundrect(x = col_prob, y = row2, width = prob_w, height = prob_h,
                 gp = gpar(fill = "white", col = "grey40"))
  grid.text(expression(P(other)), x = col_prob, y = row2,
            gp = gpar(fontsize = 9))

  grid.roundrect(x = col_prob, y = row3, width = prob_w, height = prob_h,
                 gp = gpar(fill = "white", col = "grey40"))
  grid.text(expression(P(NPL)), x = col_prob, y = row3,
            gp = gpar(fontsize = 9))

  # --- Arrows: choice rule -> probabilities ---
  grid.lines(x = c(col_rule + rule_w/2, col_prob - prob_w/2),
             y = c(row1 - 0.05, row1), gp = arrow_gp, arrow = arr)
  grid.lines(x = c(col_rule + rule_w/2, col_prob - prob_w/2),
             y = c(row2, row2), gp = arrow_gp, arrow = arr)
  grid.lines(x = c(col_rule + rule_w/2, col_prob - prob_w/2),
             y = c(row3 + 0.05, row3), gp = arrow_gp, arrow = arr)

  # --- Option counts annotation ---
  grid.text(expression(paste(n[k], " options")),
            x = col_prob, y = 0.08, gp = gpar(fontsize = 8, col = "grey50"))

  popViewport()
  dev.off()
}

create_m3_schematic()
cat("Figure saved to figures/m3_framework_schematic.pdf\n")
