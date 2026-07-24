# Visualizations for the paper examples

# Open a larger interactive device when the current Plots pane is too small.
ensure_interactive_plot_device <- function(width, height) {
  if (!interactive()) return(invisible(FALSE))

  current_size <- grDevices::dev.size("in")
  needs_device <- grDevices::dev.cur() == 1L ||
    any(!is.finite(current_size)) ||
    current_size[1] < 0.8 * width ||
    current_size[2] < 0.8 * height

  if (needs_device) {
    grDevices::dev.new(width = width, height = height)
  }
  invisible(needs_device)
}

# Cross-validation criterion across gamma values.
plot_cv_risk <- function(
    cv_result, gamma_grid, file = NULL, x_max = 50, rdelta = NULL,
    ci_lower = NULL, ci_upper = NULL,
    rdelta_lower = NULL, rdelta_upper = NULL) {
  foldwise <- as.matrix(cv_result$foldwise_rdelta)
  if (ncol(foldwise) != length(gamma_grid)) {
    if (nrow(foldwise) == length(gamma_grid)) {
      foldwise <- t(foldwise)
    } else {
      stop("The foldwise results do not match gamma_grid.")
    }
  }

  # Put gamma on an increasing scale before plotting.
  gamma_order <- order(gamma_grid)
  gamma <- gamma_grid[gamma_order]
  foldwise <- foldwise[, gamma_order, drop = FALSE]
  average <- as.numeric(cv_result$mean_abs_rdelta)[gamma_order]
  reorder_interval <- function(values, label) {
    if (is.null(values)) return(NULL)
    if (length(values) != length(gamma_grid)) {
      stop(label, " must have the same length as gamma_grid.")
    }
    as.numeric(values)[gamma_order]
  }
  ci_lower <- reorder_interval(ci_lower, "ci_lower")
  ci_upper <- reorder_interval(ci_upper, "ci_upper")
  rdelta_lower <- reorder_interval(rdelta_lower, "rdelta_lower")
  rdelta_upper <- reorder_interval(rdelta_upper, "rdelta_upper")
  if (!is.null(rdelta)) {
    if (length(rdelta) != length(gamma_grid)) {
      stop("rdelta must have the same length as gamma_grid.")
    }
    rdelta <- as.numeric(rdelta)[gamma_order]
  }
  best <- which.min(abs(gamma - as.numeric(cv_result$best_gamma)))
  x_max <- min(x_max, max(gamma))
  keep <- gamma <= x_max
  plotted_risks <- c(
    foldwise[, keep], average[keep], rdelta[keep],
    ci_lower[keep], ci_upper[keep],
    rdelta_lower[keep], rdelta_upper[keep]
  )
  risk_range <- range(0, plotted_risks, finite = TRUE)
  y_padding <- max(0.02, 0.08 * diff(risk_range))
  y_min <- if (risk_range[1] < 0) risk_range[1] - y_padding else 0
  y_max <- risk_range[2] + y_padding

  # Open a PDF only when an output path is supplied.
  if (!is.null(file)) {
    grDevices::pdf(file, width = 7, height = 5, family = "serif", useDingbats = FALSE)
  } else {
    ensure_interactive_plot_device(width = 7, height = 5)
    layout(1)
    par(mfrow = c(1, 1))
  }
  on.exit(if (!is.null(file)) grDevices::dev.off(), add = TRUE)
  interactive_plot <- is.null(file)
  old_par <- graphics::par(
    family = "serif", lend = "round", ljoin = "round",
    mar = if (interactive_plot) {
      c(2.2, 2.4, 0.4, 0.4)
    } else {
      c(8.2, 5.2, 1.2, 1.2)
    },
    mgp = if (interactive_plot) c(1.25, 0.35, 0) else c(2.7, 0.8, 0),
    cex.axis = if (interactive_plot) 0.72 else 0.9,
    cex.lab = if (interactive_plot) 0.85 else 1.15
  )
  on.exit(graphics::par(old_par), add = TRUE)

  plot(
    gamma[keep], average[keep], type = "n",
    xlim = c(0, x_max), ylim = c(y_min, y_max),
    xlab = expression(gamma),
    ylab = if (is.null(rdelta)) {
      expression(widehat(R)[Delta]^"+")
    } else {
      expression(widehat(R)[Delta]^"+"~"and"~widehat(R)[Delta])
    },
    axes = FALSE, xaxs = "i", yaxs = "i"
  )
  x_ticks <- pretty(c(0, x_max), n = 6)
  axis(1, at = x_ticks, labels = FALSE)
  mtext(x_ticks, side = 1, at = x_ticks, line = 0.6, cex = 0.9)
  axis(2, las = 1)
  box(bty = "l", lwd = 0.7)
  if (!is.null(ci_lower) && !is.null(ci_upper)) {
    polygon(
      c(gamma[keep], rev(gamma[keep])),
      c(ci_lower[keep], rev(ci_upper[keep])),
      border = NA, col = grDevices::adjustcolor("grey45", alpha.f = 0.18)
    )
  }
  if (!is.null(rdelta_lower) && !is.null(rdelta_upper)) {
    polygon(
      c(gamma[keep], rev(gamma[keep])),
      c(rdelta_lower[keep], rev(rdelta_upper[keep])),
      border = NA, col = grDevices::adjustcolor("#D6523C", alpha.f = 0.14)
    )
  }
  # Fold-specific curves stay in the background; the average is highlighted.
  for (i in seq_len(nrow(foldwise))) {
    lines(gamma[keep], foldwise[i, keep], lty = "dotted", lwd = 0.7, col = "grey55")
  }
  lines(gamma[keep], average[keep], lty = "solid", lwd = 1.2, col = "black")
  if (!is.null(rdelta)) {
    abline(h = 0, lty = "dashed", lwd = 0.6, col = "grey70")
    lines(
      gamma[keep], rdelta[keep],
      lty = "longdash", lwd = 1.2, col = "#D6523C"
    )
  }
  abline(v = gamma[best], lty = "dotdash", lwd = 0.75, col = "grey25")
  text(
    min(gamma[best] + 1, x_max - 8), y_max,
    labels = bquote(gamma[CV] == .(round(gamma[best], 2))),
    adj = c(0, 1), cex = 0.9
  )
  legend_labels <- c("Average", "Fold-specific")
  legend_types <- c("solid", "dotted")
  legend_widths <- c(1.2, 0.8)
  legend_colours <- c("black", "grey55")
  if (!is.null(rdelta)) {
    legend_labels <- c(legend_labels, "Rdelta")
    legend_types <- c(legend_types, "longdash")
    legend_widths <- c(legend_widths, 1.2)
    legend_colours <- c(legend_colours, "#D6523C")
  }
  if (interactive_plot) {
    legend(
      "bottomright", inset = 0.02,
      legend = legend_labels,
      lty = legend_types, lwd = legend_widths,
      col = legend_colours, bty = "n", cex = 0.75
    )
  } else {
    plot_limits <- par("usr")
    legend(
      x = mean(plot_limits[1:2]),
      y = plot_limits[3] - 0.24 * diff(plot_limits[3:4]),
      xjust = 0.5, yjust = 1, xpd = NA, horiz = TRUE,
      legend = legend_labels,
      lty = legend_types, lwd = legend_widths,
      col = legend_colours, bty = "n", cex = 0.8
    )
  }
  invisible(file)
}

plot_beta_path <- function(
    beta_gamma, gamma_grid, beta_cp, best_gamma, file = NULL,
    beta_true = NULL, x_max = 50,
    beta_lower = NULL, beta_upper = NULL,
    legend_ncol = NULL, legend_cex = 0.75,
    pdf_width = 7, pdf_height = 5, bottom_margin = 8.2,
    y_limits = NULL, coefficient_names = NULL,
    reference_label = "star", truth_label = "true") {
  # Keep the same colours across Examples 1 and 2.
  colours <- c("#D6523C", "#0077B6")
  if (nrow(beta_gamma) > 2) {
    colours <- grDevices::hcl.colors(nrow(beta_gamma), "Dark 3")
  }
  keep <- gamma_grid <= min(x_max, max(gamma_grid))
  x_max <- min(x_max, max(gamma_grid))
  best <- which.min(abs(gamma_grid - best_gamma))
  if (xor(is.null(beta_lower), is.null(beta_upper))) {
    stop("Both beta_lower and beta_upper must be supplied.")
  }
  if (!is.null(coefficient_names) &&
      length(coefficient_names) != nrow(beta_gamma)) {
    stop("coefficient_names must match the number of coefficient paths.")
  }
  if (!reference_label %in% c("star", "RL", "CP")) {
    stop("reference_label must be 'star', 'RL', or 'CP'.")
  }
  if (!truth_label %in% c("true", "CP")) {
    stop("truth_label must be either 'true' or 'CP'.")
  }
  if (!is.null(beta_lower) &&
      (!all(dim(beta_lower) == dim(beta_gamma)) ||
       !all(dim(beta_upper) == dim(beta_gamma)))) {
    stop("The beta confidence bands must match beta_gamma.")
  }
  data_limits <- range(
    beta_gamma[, keep], beta_cp, beta_true,
    if (!is.null(beta_lower)) beta_lower[, keep],
    if (!is.null(beta_upper)) beta_upper[, keep],
    finite = TRUE
  )
  limits <- if (is.null(y_limits)) {
    data_limits + c(-0.08, 0.16)
  } else {
    if (length(y_limits) != 2L || any(!is.finite(y_limits)) ||
        y_limits[1] >= y_limits[2]) {
      stop("y_limits must contain two increasing finite values.")
    }
    as.numeric(y_limits)
  }

  if (!is.null(file)) {
    grDevices::pdf(
      file, width = pdf_width, height = pdf_height,
      family = "serif", useDingbats = FALSE
    )
  } else {
    ensure_interactive_plot_device(
      width = pdf_width, height = pdf_height
    )
    layout(1)
    par(mfrow = c(1, 1))
  }
  on.exit(if (!is.null(file)) grDevices::dev.off(), add = TRUE)
  interactive_plot <- is.null(file)
  old_par <- graphics::par(
    family = "serif", lend = "round", ljoin = "round",
    mar = if (interactive_plot) {
      c(2.2, 2.4, 0.4, 0.4)
    } else {
      c(bottom_margin, 5.2, 1.2, 1.2)
    },
    mgp = if (interactive_plot) c(1.25, 0.35, 0) else c(2.7, 0.8, 0),
    cex.axis = if (interactive_plot) 0.72 else 0.9,
    cex.lab = if (interactive_plot) 0.85 else 1.15
  )
  on.exit(graphics::par(old_par), add = TRUE)

  matplot(
    gamma_grid[keep], t(beta_gamma[, keep, drop = FALSE]),
    type = "l", lty = "solid", lwd = 1.2, col = colours,
    xlim = c(0, x_max), ylim = limits,
    xlab = expression(gamma), ylab = expression(beta(gamma)),
    axes = FALSE, xaxs = "i"
  )
  if (!is.null(beta_lower) && !is.null(beta_upper)) {
    for (j in seq_len(nrow(beta_gamma))) {
      polygon(
        c(gamma_grid[keep], rev(gamma_grid[keep])),
        c(beta_lower[j, keep], rev(beta_upper[j, keep])),
        border = NA, col = grDevices::adjustcolor(colours[j], alpha.f = 0.16)
      )
    }
    matlines(
      gamma_grid[keep], t(beta_gamma[, keep, drop = FALSE]),
      lty = "solid", lwd = 1.2, col = colours
    )
  }
  axis(1, at = seq(0, x_max, by = 10))
  axis(2, las = 1)
  box(bty = "l", lwd = 0.7)
  abline(h = beta_cp, col = colours, lty = "dotted", lwd = 0.8)
  if (!is.null(beta_true)) {
    abline(h = beta_true, col = colours, lty = "longdash", lwd = 0.8)
  }
  abline(v = gamma_grid[best], lty = "dotdash", lwd = 0.75, col = "grey25")
  text(
    min(gamma_grid[best] + 1, x_max - 13), limits[2],
    labels = bquote(gamma[CV] == .(round(gamma_grid[best], 2))),
    adj = c(0, 1), cex = 0.9
  )

  # Build mathematical labels for the legend.
  path_labels <- lapply(seq_len(nrow(beta_gamma)), function(i) {
    index <- if (is.null(coefficient_names)) i else as.name(coefficient_names[i])
    bquote(hat(beta)[.(index)](gamma))
  })
  cp_labels <- lapply(seq_len(nrow(beta_gamma)), function(i) {
    index <- if (is.null(coefficient_names)) i else as.name(coefficient_names[i])
    if (reference_label %in% c("RL", "CP")) {
      suffix <- as.name(reference_label)
      bquote(hat(beta)[.(index) * "," * .(suffix)])
    } else {
      bquote(hat(beta)[.(index)]^"*")
    }
  })
  legend_labels <- as.expression(c(path_labels, cp_labels))
  legend_colours <- c(colours, colours)
  legend_types <- c(rep("solid", nrow(beta_gamma)), rep("dotted", nrow(beta_gamma)))

  if (!is.null(beta_true)) {
    true_labels <- lapply(seq_len(nrow(beta_gamma)), function(i) {
      index <- if (is.null(coefficient_names)) i else as.name(coefficient_names[i])
      if (truth_label == "CP") {
        bquote(beta[.(index) * ",CP"])
      } else {
        bquote(beta[.(index) * ",true"])
      }
    })
    legend_labels <- as.expression(c(path_labels, cp_labels, true_labels))
    legend_colours <- c(colours, colours, colours)
    legend_types <- c(
      rep("solid", nrow(beta_gamma)), rep("dotted", nrow(beta_gamma)),
      rep("longdash", nrow(beta_gamma))
    )
  }
  legend_arguments <- list(
    legend = legend_labels, col = legend_colours,
    lty = legend_types, lwd = c(rep(1.2, nrow(beta_gamma)), rep(0.8, length(legend_types) - nrow(beta_gamma))),
    bty = "n", cex = legend_cex
  )
  if (interactive_plot) {
    legend_arguments$x <- "bottom"
    legend_arguments$inset <- 0.01
    legend_arguments$xpd <- FALSE
  } else {
    plot_limits <- par("usr")
    legend_arguments$x <- mean(plot_limits[1:2])
    legend_arguments$y <- plot_limits[3] - 0.24 * diff(plot_limits[3:4])
    legend_arguments$xjust <- 0.5
    legend_arguments$yjust <- 1
    legend_arguments$xpd <- NA
  }
  if (is.null(legend_ncol)) {
    legend_arguments$horiz <- TRUE
  } else {
    legend_arguments$ncol <- legend_ncol
  }
  do.call(legend, legend_arguments)
  invisible(file)
}

# Direct-cause graph with edge width proportional to coefficient strength.
plot_direct_cause_graph <- function(
    causes, beta_hat, response_sd, file = NULL) {
  if (length(causes) != length(beta_hat)) {
    stop("causes and beta_hat must have the same length.")
  }
  if (!is.finite(response_sd) || response_sd <= 0) {
    stop("response_sd must be positive.")
  }

  standardized_effect <- as.numeric(beta_hat) / response_sd
  magnitude <- abs(standardized_effect)
  edge_width <- if (max(magnitude) > 0) {
    1 + 7 * magnitude / max(magnitude)
  } else {
    rep(1, length(magnitude))
  }
  edge_colour <- ifelse(
    standardized_effect >= 0, "#D6523C", "#0077B6"
  )

  if (!is.null(file)) {
    grDevices::pdf(
      file, width = 8, height = 5.5,
      family = "serif", useDingbats = FALSE
    )
  }
  on.exit(if (!is.null(file)) grDevices::dev.off(), add = TRUE)
  old_par <- graphics::par(
    family = "serif", mar = c(1.5, 1.5, 1.5, 1.5),
    lend = "round", xpd = NA
  )
  on.exit(graphics::par(old_par), add = TRUE)

  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1), asp = 1)
  cause_y <- seq(0.86, 0.14, length.out = length(causes))
  target_y <- 0.5

  for (i in seq_along(causes)) {
    arrows(
      0.28, cause_y[i], 0.71, target_y,
      length = 0.1, angle = 22,
      lwd = edge_width[i], col = edge_colour[i]
    )
    text(
      0.49, (cause_y[i] + target_y) / 2,
      labels = sprintf("%+.3f", standardized_effect[i]),
      col = edge_colour[i], cex = 0.82,
      pos = if (cause_y[i] >= target_y) 3 else 1,
      offset = 0.25
    )
  }

  rect(0.08, cause_y - 0.045, 0.27, cause_y + 0.045, col = "white")
  cause_labels <- as.expression(lapply(causes, function(cause) {
    if (grepl("^l_[0-9]+$", cause)) {
      index <- sub("^l_", "", cause)
      bquote(L[.(index)])
    } else {
      as.name(cause)
    }
  }))
  text(0.175, cause_y, labels = cause_labels, cex = 0.95)
  rect(0.72, target_y - 0.075, 0.9, target_y + 0.075, col = "white")
  text(0.81, target_y, labels = expression(Y == tilde(I)[2]), cex = 1.1)

  legend(
    "bottom", inset = -0.02, horiz = TRUE, bty = "n",
    legend = c("Positive", "Negative"),
    col = c("#D6523C", "#0077B6"), lwd = 3, cex = 0.85
  )
  invisible(data.frame(
    cause = causes,
    beta_hat = as.numeric(beta_hat),
    standardized_effect = standardized_effect,
    edge_width = edge_width
  ))
}
