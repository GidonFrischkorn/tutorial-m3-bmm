# Plotting theme and color-blind friendly palette for M3 tutorial paper
#
# Figure guidelines: Helvetica/Arial font, 9pt tick labels, 10pt axis
# titles at final published size. Figures are designed for two-column width
# (~6.875 inches). With base_size = 11, axis text at rel(0.82) ≈ 9pt and
# axis titles at 11pt.
#
# Usage:
#   source(here("functions", "clean_plot.R"))
#
#   # basic usage
#   ggplot(data, aes(x, y)) + geom_point() + clean_plot()
#
#   # with ad-hoc overrides
#   ggplot(data, aes(x, y)) + geom_point() +
#     clean_plot(axis.text.y = element_blank())
#
#   # with color palette
#   ggplot(data, aes(x, y, color = condition)) + geom_point() +
#     scale_color_m3() + clean_plot()

###############################################################################!
# Color Palette ----------------------------------------------------------------
###############################################################################!

# Okabe-Ito color-blind friendly palette
m3_palette <- c(
  "#0072B2",
  "#E69F00",
  "#009E73",
  "#D55E00",
  "#56B4E9",
  "#CC79A7",
  "#F0E442",
  "#000000"
)

# convenience scale functions for consistent color mapping across tutorials
scale_color_m3 <- function(...) {
  scale_color_manual(values = m3_palette, ...)
}

scale_fill_m3 <- function(...) {
  scale_fill_manual(values = m3_palette, ...)
}

###############################################################################!
# Plot Theme -------------------------------------------------------------------
###############################################################################!

# base_size  - overall text size scaling (default 11pt, suitable for figures
#              saved at two-column width ~6.875 inches; use 9 for one-column
#              figures ~3.25 inches)
# base_family - font family (default "" uses ggplot2 default sans-serif;
#               set to "Arial" or "Helvetica" for APS submission)
# ...        - additional theme() arguments passed through for ad-hoc overrides
clean_plot <- function(base_size = 11, base_family = "", ...) {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      # remove grid and panel border for a clean look
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border     = element_blank(),

      # axis lines and ticks
      axis.line  = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.4),
      axis.text  = element_text(color = "black", size = rel(0.82)),

      # legend
      legend.key        = element_rect(fill = "white", color = NA),
      legend.background = element_rect(fill = "white", color = NA),

      # facet strip labels
      strip.background = element_rect(fill = "grey92", color = NA),
      strip.text       = element_text(face = "bold"),

      # plot title and subtitle
      plot.title    = element_text(face = "bold", hjust = 0),
      plot.subtitle = element_text(color = "grey30", hjust = 0),
      plot.margin   = margin(10, 10, 10, 10),

      # pass-through overrides
      ...
    )
}

###############################################################################!
# Per-Facet Y-Axis Scaling -----------------------------------------------------
###############################################################################!

# Set individual y-axis limits for each facet panel in a facet_wrap plot.
# Requires facet_wrap(scales = "free") to work. This function modifies the
# ggproto object directly — call it on a completed ggplot object, not as a
# ggplot layer.
#
# Usage:
#   p <- ggplot(data, aes(x, y)) + geom_point() + facet_wrap(~cat, scales = "free")
#   p <- scale_individual_facet_y_axes(p, ylims = list(c(0.4, 1), c(0, 0.3)))
#
# Adapted from: https://stackoverflow.com/questions/51735481
scale_individual_facet_y_axes <- function(plot, ylims) {
  init_scales_orig <- plot$facet$init_scales

  init_scales_new <- function(...) {
    r <- init_scales_orig(...)
    y <- r$y
    if (is.null(y)) return(r)
    for (i in seq_along(y)) {
      ylim <- ylims[[i]]
      if (!is.null(ylim)) {
        y[[i]]$limits <- ylim
      }
    }
    r$y <- y
    return(r)
  }

  plot$facet$init_scales <- init_scales_new
  return(plot)
}

###############################################################################!
# Figure Dimensions (reference) ------------------------------------------------
###############################################################################!

# Two-column width: 6.875 inches (17.5 cm)
# One-column width: 3.25 inches (8.45 cm)
#
# Recommended ggsave settings for this project:
#
#   Single panel:
#     ggsave(filename, plot, width = 6.5, height = 5, dpi = 300)
#
#   Two-panel (side by side via patchwork):
#     ggsave(filename, plot, width = 6.5, height = 4, dpi = 300)
#
#   Multi-panel (2x2 grid):
#     ggsave(filename, plot, width = 6.5, height = 6, dpi = 300)
#
# Save as PDF for vector graphics (preferred by AMPPS/APS):
#     ggsave(filename, plot, width = 6.5, height = 5, device = "pdf")
