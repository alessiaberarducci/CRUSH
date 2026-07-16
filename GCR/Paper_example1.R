###################################
###### Paper example 2. ###########
###################################

source("./cross_valid_fun.R")
set.seed(123)

n<-100
############################
## Environment e1
############################

x1_e1 <- rexp(n, rate = 1)
y_e1  <- rpois(n, lambda = x1_e1)

eps_2_e1 <- runif(n, min = -0.5, max = 0.5)

x2_e1 <- -y_e1 +
  eps_2_e1


############################
## Environment e2
############################

x1_e2 <- rpois(n, lambda = 1)

eps_Y_e2 <- runif(n, min = -1, max = 1)
y_e2 <- log1p(x1_e2) + eps_Y_e2

eps_2_e2 <- runif(n, min = -0.5, max = 0.5)

x2_e2 <- -0.5 * y_e2 +
  eps_2_e2


############################
## Predictor matrices
############################

X_e1 <- cbind(
  X1 = x1_e1,
  X2 = x2_e1
)

X_e2 <- cbind(
  X1 = x1_e2,
  X2 = x2_e2
)
############################
## Combined data object
############################

data <- list(
  Xe = X_e2,
  ye = y_e2,
  Xo = X_e1,
  yo = y_e1
)



############################
## Cross-validation
############################

gamma_grid <- seq(0, 300, by = 1.5)

estimators <- setNames(
  lapply(gamma_grid, make_estimator),
  gamma_grid
)

res <- cross_validation(
  folds      = 5,
  data       = data,
  estimators = estimators
)

best_gamma <- as.numeric(res$best_gamma)

cat("Best gamma:", best_gamma, "\n")



mean_abs_rdelta     <- res$mean_abs_rdelta
foldwise_rdelta_mat <- res$foldwise_rdelta            # rows = γ, cols = fold



m <- moments(data)

beta_cp<-compute_cd(m)
beta_ols<-compute_ols(m)



n <- length(gamma_grid)


dims <- dim(m$Gplus)


beta_gamma <- matrix(NA, nrow = dims[1], ncol = n)
eigen_Gdelta <- eigen(m$Gdelta)

V <- eigen_Gdelta$vectors              
Lambda <- diag(abs(eigen_Gdelta$values)) 

G_delta_plus <- V %*% Lambda %*% t(V)
m$Gdelta_plus <- G_delta_plus

for (i in seq_along(gamma_grid)) {
  Gtemp <- m$Gplus +  gamma_grid[i] * m$Gdelta_plus
  Ztemp <- m$Zplus +  gamma_grid[i] * m$Gdelta_plus %*% beta_cp
  beta_gamma[, i] <- solve(Gtemp, Ztemp)
}





######## CV plot

library(ggplot2)
library(dplyr)
library(grid)
library(scales)

# ==================================================
# Basic objects
# ==================================================

best_gamma_num <- as.numeric(best_gamma)
m <- as.matrix(res$foldwise_rdelta)

# --------------------------------------------------
# Determine matrix orientation
#
# Required final orientation:
# rows    = folds
# columns = gamma values
# --------------------------------------------------

if (ncol(m) == length(gamma_grid)) {
  
  m_plot <- m
  gamma_fold <- as.numeric(gamma_grid)
  
} else if (nrow(m) == length(gamma_grid)) {
  
  m_plot <- t(m)
  gamma_fold <- as.numeric(gamma_grid)
  
} else {
  
  # Try to recover gamma values from matrix column names
  gamma_from_columns <- suppressWarnings(
    as.numeric(colnames(m))
  )
  
  gamma_from_rows <- suppressWarnings(
    as.numeric(rownames(m))
  )
  
  if (
    !is.null(colnames(m)) &&
    length(gamma_from_columns) == ncol(m) &&
    all(!is.na(gamma_from_columns))
  ) {
    
    m_plot <- m
    gamma_fold <- gamma_from_columns
    
  } else if (
    !is.null(rownames(m)) &&
    length(gamma_from_rows) == nrow(m) &&
    all(!is.na(gamma_from_rows))
  ) {
    
    m_plot <- t(m)
    gamma_fold <- gamma_from_rows
    
  } else {
    
    stop(
      paste0(
        "The dimensions of res$foldwise_rdelta are ",
        nrow(m), " x ", ncol(m),
        ", whereas gamma_grid has length ",
        length(gamma_grid),
        ". The number of rows or columns must match ",
        "the length of gamma_grid."
      )
    )
  }
}

# Ensure increasing gamma order
gamma_order <- order(gamma_fold)

gamma_fold <- gamma_fold[gamma_order]
m_plot <- m_plot[, gamma_order, drop = FALSE]

# ==================================================
# Average criterion
# ==================================================

# Use mean_abs_rdelta when it has the correct length.
# Otherwise, compute the foldwise average directly.
if (length(mean_abs_rdelta) == length(gamma_fold)) {
  
  average_risk <- as.numeric(mean_abs_rdelta)[gamma_order]
  
} else {
  
  warning(
    paste0(
      "mean_abs_rdelta has length ",
      length(mean_abs_rdelta),
      ", but the foldwise matrix contains ",
      length(gamma_fold),
      " gamma values. The average is therefore ",
      "computed directly from res$foldwise_rdelta."
    )
  )
  
  average_risk <- colMeans(
    abs(m_plot),
    na.rm = TRUE
  )
}

# ==================================================
# Prepare data frames
# ==================================================

df_avg <- data.frame(
  gamma = gamma_fold,
  abs_risk = average_risk
)

df_fold <- data.frame(
  fold = factor(
    rep(
      seq_len(nrow(m_plot)),
      each = ncol(m_plot)
    )
  ),
  gamma = rep(
    gamma_fold,
    times = nrow(m_plot)
  ),
  abs_risk = as.vector(
    t(abs(m_plot))
  )
)

# ==================================================
# Selected gamma
# ==================================================

best_idx <- which.min(
  abs(gamma_fold - best_gamma_num)
)

best_gamma_plot <- gamma_fold[best_idx]
best_risk_plot <- average_risk[best_idx]

# ==================================================
# Plot range
# ==================================================

x_max_plot <- min(
  50,
  max(gamma_fold, na.rm = TRUE)
)

df_avg_plot <- df_avg |>
  filter(gamma <= x_max_plot)

df_fold_plot <- df_fold |>
  filter(gamma <= x_max_plot)

# Dynamic y-axis range
y_min_plot <- 0

y_max_data <- max(
  c(
    df_fold_plot$abs_risk,
    df_avg_plot$abs_risk
  ),
  na.rm = TRUE
)

y_padding <- max(
  0.02,
  0.08 * y_max_data
)

y_max_plot <- y_max_data + y_padding

# ==================================================
# Colours
# ==================================================

col_average <- "black"
col_folds   <- "grey55"
col_best    <- "black"

# ==================================================
# Plot
# ==================================================

p_risk <- ggplot() +
  
  # Fold-specific curves
  geom_line(
    data = df_fold_plot,
    mapping = aes(
      x = gamma,
      y = abs_risk,
      group = fold,
      linetype = "Fold-specific"
    ),
    colour = col_folds,
    linewidth = 0.7,
    alpha = 0.65
  ) +
  
  # Average curve
  geom_line(
    data = df_avg_plot,
    mapping = aes(
      x = gamma,
      y = abs_risk,
      linetype = "Average"
    ),
    colour = col_average,
    linewidth = 1.2
  ) +
  
  # Points on the average curve
  geom_point(
    data = df_avg_plot,
    mapping = aes(
      x = gamma,
      y = abs_risk
    ),
    inherit.aes = FALSE,
    colour = col_average,
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
  
  # Highlighted point at gamma*
  annotate(
    geom = "point",
    x = best_gamma_plot,
    y = best_risk_plot,
    colour = col_best,
    size = 3.8,
    shape = 16
  ) +
  
  # gamma* annotation
  annotate(
    geom = "text",
    x = min(
      best_gamma_plot + 1,
      x_max_plot - 8
    ),
    y = y_max_plot,
    label = paste0(
      "gamma^'*' == ",
      round(best_gamma_plot, 2)
    ),
    parse = TRUE,
    colour = col_best,
    hjust = 0,
    vjust = 1,
    size = 4,
    family = "serif"
  ) +
  
  # ------------------------------------------------
# Linetype scale
# ------------------------------------------------

scale_linetype_manual(
  name = NULL,
  breaks = c(
    "Average",
    "Fold-specific"
  ),
  values = c(
    "Average" = "solid",
    "Fold-specific" = "dotted"
  )
) +
  
  # ------------------------------------------------
# Axes
# ------------------------------------------------

scale_x_continuous(
  breaks = pretty_breaks(n = 6),
  minor_breaks = NULL,
  expand = expansion(
    mult = c(0.01, 0.03)
  )
) +
  
  scale_y_continuous(
    breaks = pretty_breaks(n = 7),
    labels = label_number(
      accuracy = 0.01,
      trim = TRUE
    ),
    expand = expansion(
      mult = c(0.02, 0.04)
    )
  ) +
  
  coord_cartesian(
    xlim = c(
      0,
      x_max_plot
    ),
    ylim = c(
      y_min_plot,
      y_max_plot
    ),
    clip = "off"
  ) +
  
  labs(
    x = expression(gamma),
    y = expression(
      widehat(R)[Delta]^"+"
    )
  ) +
  
  # ------------------------------------------------
# Legend
# ------------------------------------------------

guides(
  linetype = guide_legend(
    nrow = 1,
    byrow = TRUE,
    order = 1,
    override.aes = list(
      colour = c(
        col_average,
        col_folds
      ),
      linewidth = c(
        1.2,
        0.8
      ),
      shape = c(
        16,
        NA
      ),
      size = c(
        2.4,
        NA
      ),
      alpha = 1
    )
  )
) +
  
  # ------------------------------------------------
# Theme matching the beta-path plot
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

# Display plot
p_risk

# ==================================================
# Save PDF
# ==================================================

out_risk <- file.path(
  Sys.getenv("HOME"),
  "/Users/alessiaberarducci/Library/CloudStorage/OneDrive-USI/File di Lembo Melania - Code CR New/visualization",
  "example1_risk.pdf"
)

ggsave(
  filename = out_risk,
  plot = p_risk,
  device = "pdf",
  width = 7,
  height = 5,
  units = "in",
  useDingbats = FALSE
)

file.exists(out_risk)
out_risk








