# Plot theme and color-blind friendly palette for the M3 tutorial figures.
# Usage:
#   source(here("functions", "clean_plot.R"))
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

# base_size   - base text size (11 for full-width figures, 9 for ~3.25 in wide)
# base_family - font family ("" uses the ggplot2 default sans-serif)
# ...         - additional theme() arguments for ad-hoc overrides
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

# Set y-axis limits per panel of a facet_wrap() plot with free y scales
# (ylims: one c(min, max) per panel, in panel order). It modifies the facet
# object, so call it on a finished ggplot object, not as a layer.
# Adapted from https://stackoverflow.com/questions/51735481
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
