# Figure 2: stacked bars showing how the activation sources b, a, and c combine
# for each response option in a simple span cued recall task. Uses the same
# letters as the task example (Figure 1).
# Output: figures/m3_activations.pdf and .png

###############################################################################!
# 0) R Setup -------------------------------------------------------------------
###############################################################################!

pacman::p_load(here, ggplot2)
source(here("functions", "clean_plot.R"))

dir.create(here("figures"), showWarnings = FALSE)

## 0.1) Illustrative activation values -----------------------------------------
b_val <- 0.5   # baseline
a_val <- 1.5   # item memory
c_val <- 2.0   # context binding

###############################################################################!
# 1) Build Figure --------------------------------------------------------------
###############################################################################!

## 1.1) Response options by category -------------------------------------------
items      <- c("D", "K", "L", "S", "B", "Y",
                "F", "G", "H", "J", "M")
categories <- c("Correct", rep("Other", 5), rep("NPL", 5))

# Build stacked data in long format
act_data <- do.call(rbind, lapply(seq_along(items), function(i) {
  cat_i <- categories[i]
  vals <- if (cat_i == "Correct") {
    c(b_val, a_val, c_val)
  } else if (cat_i == "Other") {
    c(b_val, a_val, 0)
  } else {
    c(b_val, 0, 0)
  }
  data.frame(
    item       = items[i],
    category   = cat_i,
    source     = c("b (baseline)",
                   "a (item memory)",
                   "c (context binding)"),
    activation = vals
  )
}))

act_data$item <- factor(act_data$item, levels = items)
act_data$source <- factor(
  act_data$source,
  levels = c("c (context binding)",
             "a (item memory)",
             "b (baseline)")
)
act_data$category <- factor(
  act_data$category,
  levels = c("Correct", "Other", "NPL")
)

# Colors
source_colors <- c(
  "b (baseline)"        = "grey75",
  "a (item memory)"     = m3_palette[1],
  "c (context binding)" = m3_palette[2]
)

## 1.2) Formula annotations (one per category) ---------------------------------
formula_data <- data.frame(
  category = factor(c("Correct", "Other", "NPL"),
                    levels = c("Correct", "Other", "NPL")),
  item     = factor(c("D", "S", "H"), levels = items),
  y        = c(b_val + a_val + c_val + 0.3,
               b_val + a_val + 0.3,
               b_val + 0.3),
  label    = c("A == b + a + c",
               "A == b + a",
               "A == b")
)

## 1.3) Plot -------------------------------------------------------------------

# Facets with space = "free_x" keep bar widths equal across categories
act_plot <- ggplot(act_data,
                   aes(x = item, y = activation, fill = source)) +
  geom_col(position = "stack", width = 0.7,
           color = "white", linewidth = 0.3) +
  geom_text(data = formula_data[formula_data$category != "Correct", ],
            aes(x = item, y = y, label = label),
            inherit.aes = FALSE, parse = TRUE, size = 3.2) +
  geom_text(data = formula_data[formula_data$category == "Correct", ],
            aes(x = item, y = y, label = label),
            inherit.aes = FALSE, parse = TRUE, size = 3.2,
            hjust = 0.35) +
  facet_grid(~ category, scales = "free_x", space = "free_x",
             switch = "x") +
  scale_fill_manual(
    values = source_colors,
    name   = "Activation Source",
    breaks = c("c (context binding)",
               "a (item memory)",
               "b (baseline)")
  ) +
  labs(x = NULL, y = "Activation") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  coord_cartesian(clip = "off") +
  clean_plot(
    strip.placement = "outside",
    panel.spacing   = unit(0.8, "lines"),
    legend.position = "bottom"
  ) +
  theme(
    strip.background = element_rect(fill = NA, color = NA),
    strip.text       = element_text(face = "bold", size = 10)
  )

act_plot

###############################################################################!
# 2) Save ----------------------------------------------------------------------
###############################################################################!

# PNG copy is used for the Word version of the manuscript
ggsave(here("figures", "m3_activations.pdf"),
       act_plot, width = 6.5, height = 4)
ggsave(here("figures", "m3_activations.png"),
       act_plot, width = 6.5, height = 4, dpi = 300)
