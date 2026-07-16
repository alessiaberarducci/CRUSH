################## FULL LEGEND
######### BETA1 E BETA2

library(ggplot2)
library(tidyr)
library(dplyr)
library(grid)
library(scales)

best_gamma_num <- as.numeric(best_gamma)

# --------------------------------------------------
# Colours
# --------------------------------------------------

col_beta1 <- "#D6523C"  # red
col_beta2 <- "#0077B6"  # blue
col_best  <- "black"    # gamma* points

col_true_b1 <- "#D6523C" # true beta1 values
col_true_b2 <- "#0077B6" # true beta2 values

# --------------------------------------------------
# Prepare data
# --------------------------------------------------

df_beta <- data.frame(
  gamma = gamma_grid,
  beta1 = beta_gamma[1, ],
  beta2 = beta_gamma[2, ]
)

df_long <- df_beta |>
  pivot_longer(
    cols = c(beta1, beta2),
    names_to = "component",
    values_to = "value"
  ) |>
  mutate(
    component = factor(
      component,
      levels = c("beta1", "beta2")
    )
  )

# Closest gamma in the grid to best_gamma
best_idx <- which.min(abs(gamma_grid - best_gamma_num))
best_gamma_plot <- gamma_grid[best_idx]

df_best <- df_long |>
  filter(gamma == best_gamma_plot)

# --------------------------------------------------
# Plot range
# --------------------------------------------------

x_max_plot <- 50

df_long_plot <- df_long |>
  filter(gamma <= x_max_plot)

df_best_plot <- df_best |>
  filter(gamma <= x_max_plot)

# --------------------------------------------------
# Horizontal reference lines
# --------------------------------------------------

df_reference <- data.frame(
  reference = factor(
    c(
      "beta1_cp",
      "beta2_cp",
      "beta1_true",
      "beta2_true"
    ),
    levels = c(
      "beta1_cp",
      "beta2_cp",
      "beta1_true",
      "beta2_true"
    )
  ),
  value = c(
    beta_cp[1],
    beta_cp[2],
    1,
    0
  )
)

# --------------------------------------------------
# Dynamic y-axis limits
# --------------------------------------------------

y_min_plot <- min(
  c(df_long_plot$value, beta_cp, 0),
  na.rm = TRUE
) - 0.08

y_max_plot <- max(
  c(df_long_plot$value, beta_cp, 1),
  na.rm = TRUE
) + 0.08

# --------------------------------------------------
# Plot
# --------------------------------------------------

p <- ggplot() +
  
  # Horizontal reference lines: CP and true values
  geom_hline(
    data = df_reference,
    aes(
      yintercept = value,
      colour = reference,
      linetype = reference
    ),
    linewidth = 0.8,
    alpha = 0.85
  ) +
  
  # Regularisation paths
  geom_line(
    data = df_long_plot,
    aes(
      x = gamma,
      y = value,
      colour = component,
      group = component
    ),
    linewidth = 1.2
  ) +
  
  # Points along the paths
  geom_point(
    data = df_long_plot,
    aes(
      x = gamma,
      y = value,
      colour = component
    ),
    size = 2.4,
    alpha = 0.85
  ) +
  
  # Vertical line at gamma*
  geom_vline(
    xintercept = best_gamma_plot,
    linetype = "dotdash",
    colour = "grey25",
    linewidth = 0.75,
    show.legend = FALSE
  ) +
  
  # Highlighted points at gamma*
  geom_point(
    data = df_best_plot,
    aes(
      x = gamma,
      y = value
    ),
    colour = col_best,
    size = 3.8,
    shape = 16,
    show.legend = FALSE
  ) +
  
  # gamma* annotation
  annotate(
    "text",
    x = min(best_gamma_plot + 1, x_max_plot - 13),
    y = y_max_plot +0.04,
    label = paste0("gamma^'*' == ", round(best_gamma_plot, 2)),
    parse = TRUE,
    colour = col_best,
    hjust = 0,
    vjust = 1,
    size = 4,
    family = "serif"
  ) +
  
  # ------------------------------------------------
# Colour scale
# ------------------------------------------------

scale_colour_manual(
  name = NULL,
  breaks = c(
    "beta1",
    "beta1_cp",
    "beta1_true",
    "beta2",
    "beta2_cp",
    "beta2_true"
  ),
  values = c(
    beta1      = col_beta1,
    beta2      = col_beta2,
    beta1_cp   = col_beta1,
    beta2_cp   = col_beta2,
    beta1_true = col_true_b1,
    beta2_true = col_true_b2
  ),
  labels = expression(
    hat(beta)[1](gamma),
    beta[1]^"*",
    beta[1*","*true],
    hat(beta)[2](gamma),
    beta[2]^"*",
    beta[2*","*true]
  )
) +
  
  # ------------------------------------------------
# Linetype scale
# ------------------------------------------------

scale_linetype_manual(
  name = NULL,
  breaks = c(
    "beta1_cp",
    "beta2_cp",
    "beta1_true",
    "beta2_true"
  ),
  values = c(
    beta1_cp   = "dotted",
    beta2_cp   = "dotted",
    beta1_true = "longdash",
    beta2_true = "longdash"
  ),
  labels = expression(
    beta[1*","*CP],
    beta[2*","*CP],
    beta[1*","*true],
    beta[2*","*true]
  )
) +
  
  # ------------------------------------------------
# Axes
# ------------------------------------------------

scale_x_continuous(
  breaks = seq(0, x_max_plot, by = 10),
  minor_breaks = seq(0, x_max_plot, by = 5),
  expand = expansion(mult = c(0.01, 0.03))
) +
  
  scale_y_continuous(
    breaks = pretty_breaks(n = 9),
    labels = label_number(
      accuracy = 0.1,
      trim = TRUE
    ),
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  
  coord_cartesian(
    xlim = c(0, x_max_plot),
    ylim = c(y_min_plot, y_max_plot),
    clip = "off"
  ) +
  
  labs(
    x = expression(gamma),
    y = expression(beta(gamma))
  ) +
  
  # ------------------------------------------------
# Legend
# ------------------------------------------------

guides(
  colour = guide_legend(
    nrow = 1,
    byrow = TRUE,
    order = 1,
    override.aes = list(
      linewidth = c(
        1.2, 0.8, 0.8,
        1.2, 0.8, 0.8
      ),
      linetype = c(
        "solid", "dotted", "longdash",
        "solid", "dotted", "longdash"
      ),
      shape = c(
        16, NA, NA,
        16, NA, NA
      ),
      size = c(
        2.4, NA, NA,
        2.4, NA, NA
      ),
      alpha = 1
    )
  ),
  linetype = "none"
) +
  
  # ------------------------------------------------
# Theme
# ------------------------------------------------

theme_classic(
  base_size = 18,
  base_family = "serif"
) +
  
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    
    axis.line = element_line(
      colour = "black",
      linewidth = 0.7
    ),
    
    axis.ticks = element_line(
      colour = "black",
      linewidth = 0.6
    ),
    
    axis.ticks.length = unit(
      0.14,
      "cm"
    ),
    
    axis.title = element_text(
      size = 15,
      colour = "black"
    ),
    
    axis.title.x = element_text(
      margin = margin(t = 12)
    ),
    
    axis.title.y = element_text(
      margin = margin(r = 12)
    ),
    
    axis.text = element_text(
      size = 11,
      colour = "black"
    ),
    
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.justification = "center",
    
    legend.text = element_text(
      size = 10.5,
      colour = "black"
    ),
    
    legend.key.width = unit(
      1.2,
      "cm"
    ),
    
    legend.key.height = unit(
      0.45,
      "cm"
    ),
    
    legend.spacing.x = unit(
      0.4,
      "cm"
    ),
    
    legend.spacing.y = unit(
      0.15,
      "cm"
    ),
    
    legend.margin = margin(
      t = 6,
      b = 2
    ),
    
    legend.background = element_blank(),
    
    plot.margin = margin(
      15,
      15,
      15,
      15
    )
  )

p

# --------------------------------------------------
# Save PDF
# --------------------------------------------------

out <- file.path(
  "/Users/alessiaberarducci/Library/CloudStorage/OneDrive-USI/File di Lembo Melania - Code CR New/visualization",
  "example2.pdf"
)

ggsave(
  filename = out,
  plot = p,
  device = "pdf",
  width = 7,
  height = 5,
  units = "in",
  useDingbats = FALSE
)

file.exists(out)
out
