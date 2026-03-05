#' Generate the simple span cued recall task example figure
#'
#' This script creates a diagram illustrating a simple span cued recall trial:
#' Left panel shows encoding (6 letters at serial positions), right panel shows
#' retrieval (position probe and three response categories).

pacman::p_load(here, ggplot2, patchwork)
source(here("functions", "clean_plot.R"))

# --- Task parameters ---
letters_stim <- c("K", "L", "D", "S", "B", "Y")
positions    <- 1:6
probe_pos    <- 3   # probe position 3 -> correct answer is "D"

# ============================================================================
# Panel A: Encoding phase
# ============================================================================

enc_data <- data.frame(
  x = positions,
  y = 0,
  label = letters_stim
)

p_encoding <- ggplot(enc_data) +
  # boxes
  geom_tile(aes(x = x, y = y),
            width = 0.85, height = 0.85,
            fill = "white", color = "grey30", linewidth = 0.6) +
  # letters
  geom_text(aes(x = x, y = y, label = label),
            size = 7, fontface = "bold") +
  # position labels below
  geom_text(aes(x = x, y = y - 0.7),
            label = paste0("Pos ", positions),
            size = 2.8, color = "grey50") +
  # title
  annotate("text", x = 3.5, y = 0.85, label = "Encoding",
           size = 4.5, fontface = "bold") +
  coord_fixed(ratio = 1, xlim = c(0.2, 6.8), ylim = c(-1.2, 1.3)) +
  theme_void()

# ============================================================================
# Panel B: Retrieval phase
# ============================================================================

ret_boxes <- data.frame(
  x = positions,
  y = 0,
  label = c("", "", "?", "", "", ""),
  fill  = c("white", "white", "#FFF3CD", "white", "white", "white"),
  border = c("grey80", "grey80", "grey30", "grey80", "grey80", "grey80")
)

# Response options with category labels
resp_data <- data.frame(
  x     = c(3, 5, 6.5),
  y     = c(-2.2, -2.2, -2.2),
  label = c("D", "K", "F"),
  cat   = c("Correct", "Intrusion", "NPL"),
  col   = c(m3_palette[3], m3_palette[1], m3_palette[4])
)

p_retrieval <- ggplot() +
  # boxes (greyed out except probe)
  geom_tile(data = ret_boxes, aes(x = x, y = y),
            width = 0.85, height = 0.85,
            fill = ret_boxes$fill, color = ret_boxes$border, linewidth = 0.6) +
  # question mark at probe position
  geom_text(data = ret_boxes[probe_pos, ],
            aes(x = x, y = y), label = "?",
            size = 8, fontface = "bold", color = "grey30") +
  # position labels below boxes
  geom_text(data = data.frame(x = positions, y = rep(-0.7, 6)),
            aes(x = x, y = y),
            label = paste0("Pos ", positions),
            size = 2.8, color = "grey50") +
  # arrows from probe to response options
  geom_segment(data = resp_data,
               aes(x = 3, xend = x, y = -0.55, yend = y + 0.45),
               arrow = arrow(length = unit(0.12, "inches"), type = "closed"),
               color = resp_data$col, linewidth = 0.6) +
  # response option circles
  geom_point(data = resp_data, aes(x = x, y = y),
             shape = 21, size = 10, fill = "white",
             color = resp_data$col, stroke = 1.5) +
  # response letters
  geom_text(data = resp_data, aes(x = x, y = y, label = label),
            size = 5, fontface = "bold", color = resp_data$col) +
  # category labels below

  geom_text(data = resp_data, aes(x = x, y = y - 0.65, label = cat),
            size = 3.2, color = resp_data$col, fontface = "bold") +
  # title
  annotate("text", x = 3.5, y = 0.85, label = "Retrieval",
           size = 4.5, fontface = "bold") +
  coord_fixed(ratio = 1, xlim = c(0.2, 7.8), ylim = c(-3.2, 1.3)) +
  theme_void()

# ============================================================================
# Combine panels
# ============================================================================

task_fig <- p_encoding + p_retrieval +
  plot_layout(widths = c(1, 1.2)) +
  plot_annotation(tag_levels = "A")

task_fig

ggsave(here("figures", "task_example.pdf"),
       task_fig, width = 6.5, height = 3)

cat("Figure saved to figures/task_example.pdf\n")
