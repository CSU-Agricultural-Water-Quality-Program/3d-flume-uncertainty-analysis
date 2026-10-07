# Additive point summaries of the existing fixed-stage model; no fitting.
# Horizontal bars are the recorded gauge accuracy, not posterior stage intervals.

build_flume_point_data <- function(post, dat, prob = 0.95,
                                   flume_names = c("1" = "3D Printed 1",
                                                   "2" = "3D Printed 2",
                                                   "3" = "Fiberglass")) {
  stopifnot(all(c("Q_obs", "H_obs", "F", "sigma_H") %in% names(dat)),
            length(prob) == 1L, is.finite(prob), prob > 0, prob < 1,
            is.matrix(post$a), is.matrix(post$b), nrow(post$a) > 1L,
            identical(dim(post$a), dim(post$b)),
            is.matrix(post$Q_true), nrow(post$Q_true) == nrow(post$a),
            ncol(post$Q_true) == length(dat$H_obs), all(is.finite(post$Q_true)),
            all(dat$F >= 1), all(dat$F <= ncol(post$a)),
            length(dat$Q_obs) == length(dat$H_obs),
            length(dat$F) == length(dat$H_obs),
            all(is.finite(post$a)), all(post$a > 0),
            all(is.finite(post$b)), all(post$b > 0), all(is.finite(dat$H_obs)),
            all(dat$H_obs >= 0), max(dat$H_obs) > 0,
            all(is.finite(dat$Q_obs)),
            all(as.character(dat$F) %in% names(flume_names)),
            length(dat$sigma_H) == 1L, is.finite(dat$sigma_H), dat$sigma_H >= 0)
  ft_to_cm <- 30.48
  cfs_to_lps <- (ft_to_cm / 100)^3 * 1000
  tails <- c((1 - prob) / 2, (1 + prob) / 2)
  # Paired parameter draws at observed stages: the same estimand as the
  # parameter ribbons. Do not add process/measurement noise or use Q_true.
  mean_draws <- vapply(seq_along(dat$H_obs), function(i) {
    post$a[, dat$F[i]] * dat$H_obs[i]^post$b[, dat$F[i]]
  }, numeric(nrow(post$a)))
  bounds <- apply(mean_draws, 2, stats::quantile, probs = tails)
  data.frame(
    observation_index = seq_along(dat$H_obs),
    flume = factor(as.character(dat$F), levels = names(flume_names),
                   labels = unname(flume_names)),
    stage_cm = dat$H_obs * ft_to_cm,
    stage_lower_cm = (dat$H_obs - dat$sigma_H) * ft_to_cm,
    stage_upper_cm = (dat$H_obs + dat$sigma_H) * ft_to_cm,
    measured_flow_lps = dat$Q_obs * cfs_to_lps,
    latent_true_flow_lps = apply(post$Q_true, 2, stats::median) * cfs_to_lps,
    estimated_flow_lps = apply(mean_draws, 2, stats::median) * cfs_to_lps,
    flow_lower_lps = bounds[1, ] * cfs_to_lps,
    flow_upper_lps = bounds[2, ] * cfs_to_lps,
    interval_probability = prob,
    posterior_draws = nrow(post$a),
    interval_kind = "parameter",
    point_estimand = "observation-specific latent true flow Q_true",
    estimand = "mean rating-curve flow a*H^b conditional on observed stage"
  )
}

plot_flume_points <- function(points, manufacturer, departures = FALSE,
                              extrapolation_fraction = 0.10) {
  stopifnot(nrow(manufacturer) == 1L,
            all(c("a_m", "b_m") %in% names(manufacturer)),
            is.finite(manufacturer$a_m), manufacturer$a_m > 0,
            is.finite(manufacturer$b_m), manufacturer$b_m > 0,
            length(extrapolation_fraction) == 1L,
            is.finite(extrapolation_fraction), extrapolation_fraction >= 0)
  ft_to_cm <- 30.48
  cfs_to_lps <- (ft_to_cm / 100)^3 * 1000
  reference <- function(h_cm) {
    manufacturer$a_m * (h_cm / ft_to_cm)^manufacturer$b_m * cfs_to_lps
  }
  x_max <- max(max(points$stage_cm) * (1 + extrapolation_fraction),
               max(points$stage_upper_cm))
  curve <- data.frame(stage_cm = seq(0, x_max, length.out = 300))
  curve$flow_lps <- reference(curve$stage_cm)
  d <- points
  d$center <- d$latent_true_flow_lps
  d$lower <- d$flow_lower_lps
  d$upper <- d$flow_upper_lps
  d$measured <- d$measured_flow_lps
  if (departures) {
    # The deterministic positive denominator makes this an increasing affine
    # transform, equivalent to transforming paired parameter draws first.
    q_m <- reference(d$stage_cm)
    # Relative differences are undefined at zero manufacturer discharge.
    keep <- is.finite(q_m) & q_m > 0
    d <- d[keep, , drop = FALSE]
    q_m <- q_m[keep]
    stopifnot(nrow(d) > 0L)
    d$center <- 100 * (d$center - q_m) / q_m
    d$lower <- 100 * (d$lower - q_m) / q_m
    d$upper <- 100 * (d$upper - q_m) / q_m
    d$measured <- 100 * (d$measured - q_m) / q_m
  }
  colors <- c("3D Printed 1" = "#0072B2", "3D Printed 2" = "#D55E00",
              "Fiberglass" = "#009E73")
  p <- ggplot2::ggplot(d, ggplot2::aes(stage_cm, center))
  reference_label <- if (departures) "Manufacturer reference (0%)" else "Manufacturer equation"
  if (departures) {
    p <- p + ggplot2::geom_hline(
      ggplot2::aes(yintercept = 0, linetype = "Manufacturer reference (0%)"),
      color = "black", linewidth = 0.7
    )
  } else {
    p <- p + ggplot2::geom_line(
      data = curve, ggplot2::aes(stage_cm, flow_lps, linetype = "Manufacturer equation"),
      inherit.aes = FALSE,
      color = "black", linewidth = 0.7
    )
  }
  y_limits <- if (departures) {
    bound <- max(abs(c(d$lower, d$upper, d$center, d$measured))) * 1.15
    c(-max(bound, .Machine$double.eps), max(bound, .Machine$double.eps))
  } else {
    c(0, max(d$upper, d$center, d$measured, curve$flow_lps) * 1.05)
  }
  p +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = lower, ymax = upper, color = flume),
      width = 0.22, linewidth = 0.65, show.legend = FALSE
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = stage_lower_cm, xmax = stage_upper_cm, color = flume),
      orientation = "y", width = diff(y_limits) * 0.015,
      linewidth = 0.65, show.legend = FALSE
    ) +
    # White latent points underneath smaller measured points expose both layers.
    ggplot2::geom_point(
      ggplot2::aes(color = flume, shape = "Estimated latent true flow"),
      fill = "white", size = 2.7, stroke = 0.7, show.legend = c(color = FALSE, shape = TRUE)
    ) +
    ggplot2::geom_point(
      ggplot2::aes(y = measured, fill = flume, shape = "Measured flow"),
      color = "black", size = 1.3, stroke = 0.4, show.legend = c(fill = FALSE, shape = TRUE)
    ) +
    ggplot2::facet_wrap(ggplot2::vars(flume), nrow = 1, drop = FALSE) +
    ggplot2::scale_color_manual(values = colors, guide = "none") +
    ggplot2::scale_fill_manual(values = colors, guide = "none") +
    ggplot2::scale_shape_manual(
      name = NULL, values = c("Estimated latent true flow" = 21, "Measured flow" = 21),
      breaks = c("Estimated latent true flow", "Measured flow")
    ) +
    ggplot2::scale_linetype_manual(name = NULL, values = setNames("solid", reference_label)) +
    ggplot2::guides(
      shape = ggplot2::guide_legend(order = 1,
        override.aes = list(fill = c("white", "grey60"), color = "black", size = c(3, 2))),
      linetype = ggplot2::guide_legend(order = 2, override.aes = list(color = "black"))
    ) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0)) +
    ggplot2::coord_cartesian(xlim = c(0, x_max), ylim = y_limits, expand = FALSE) +
    ggplot2::labs(x = "Stage (cm)",
                  y = if (departures) "Relative difference from manufacturer (%)" else "Flow rate (L/s)") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      legend.position = "bottom", legend.box = "horizontal",
      strip.background = ggplot2::element_rect(fill = "grey96"),
      panel.grid.major = ggplot2::element_line(color = "grey92", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank()
    )
}

save_flume_point_plots <- function(post, dat, manufacturer, output_dir,
                                  prefix = basename(output_dir), input_csv = NULL,
                                  fit_file = NULL, prob = 0.95) {
  points <- build_flume_point_data(post, dat, prob)
  q_m <- manufacturer$a_m * (points$stage_cm / 30.48)^manufacturer$b_m *
    (30.48 / 100)^3 * 1000
  points$manufacturer_flow_lps <- q_m
  points$flow_minus_manufacturer_lps <- points$estimated_flow_lps - q_m
  points$difference_lower_lps <- points$flow_lower_lps - q_m
  points$difference_upper_lps <- points$flow_upper_lps - q_m
  points$measured_minus_manufacturer_lps <- points$measured_flow_lps - q_m
  relative <- function(x) {
    result <- rep(NA_real_, length(q_m))
    positive <- is.finite(q_m) & q_m > 0
    result[positive] <- 100 * (x[positive] - q_m[positive]) / q_m[positive]
    result
  }
  points$relative_difference_pct <- relative(points$estimated_flow_lps)
  points$latent_relative_difference_pct <- relative(points$latent_true_flow_lps)
  points$relative_lower_pct <- relative(points$flow_lower_lps)
  points$relative_upper_pct <- relative(points$flow_upper_lps)
  points$measured_relative_difference_pct <- relative(points$measured_flow_lps)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  filename <- function(stem, ext) file.path(output_dir, paste0(prefix, "_", stem, ".", ext))
  p <- plot_flume_points(points, manufacturer)
  residual <- plot_flume_points(points, manufacturer, departures = TRUE)
  for (ext in c("png", "jpg")) {
    ggplot2::ggsave(filename("flume_points_with_xy_errors", ext), p,
                   width = 11, height = 6, dpi = 300, bg = "white")
    ggplot2::ggsave(filename("flume_points_relative_difference", ext), residual,
                   width = 11, height = 6, dpi = 300, bg = "white")
  }
  utils::write.csv(points, filename("flume_point_plot_data", "csv"), row.names = FALSE)
  provenance <- function(label, path) {
    if (is.null(path)) return(character())
    sprintf("- %s: `%s`; MD5 `%s`.", label, normalizePath(path, winslash = "/"),
            unname(tools::md5sum(path)))
  }
  writeLines(c(
    "# Point-plot definitions", "",
    sprintf("- %d calibration observations; %d posterior draws per observation.",
            nrow(points), nrow(post$a)),
    "- White-filled points: posterior medians of observation-specific latent true flow Q_true, drawn underneath solid colored measured-flow points.",
    sprintf("- Vertical bars: %g%% equal-tailed parameter credible intervals of a*H^b, conditional on observed stage; the same estimand as the parameter ribbons.", 100 * prob),
    "- The point and bar estimands differ: points summarize Q_true; vertical bars summarize a*H^b and are not credible intervals for Q_true. The parameter mean need not coincide with the latent point. No additional process or measurement scatter is added to the bars.",
    sprintf("- Horizontal bars: recorded gauge accuracy +/- %g ft (+/- %g cm); no probability coverage is assigned.", dat$sigma_H, dat$sigma_H * 30.48),
    "- The current likelihood fixes stage. Horizontal accuracy bars are not propagated into vertical credible intervals.",
    sprintf("- Manufacturer: Q (cfs) = %g * H (ft)^%g; converted to cm and L/s.", manufacturer$a_m, manufacturer$b_m),
    "- Relative difference (%) = 100 * (flow - manufacturer) / manufacturer at observed stage. Plotted white points use Q_true; vertical bars use a*H^b. Both parameter-mean and latent medians are retained separately in the CSV. Negative values indicate lower flow; zero-reference rows are NA in the CSV and omitted from the percent plot.",
    "- Gauge accuracy and manufacturer uncertainty are not propagated into relative differences.",
    "- All observations retain their measured stages; no jitter, dodge, aggregation, or refitting is applied. Facets share axis limits.",
    "- The manufacturer curve extends from zero to 10% above maximum observed stage.",
    provenance("Calibration input", input_csv), provenance("Saved fit", fit_file)
  ), filename("flume_point_plot_notes", "md"))
  write_flume_figure_captions(points, manufacturer, output_dir, prefix)
  invisible(list(points = points, plot = p, departures = residual))
}

# Reader-facing captions are exported separately and never embedded in images.
write_flume_figure_captions <- function(points, manufacturer, output_dir, prefix) {
  prob <- sprintf("%g%%", 100 * unique(points$interval_probability))
  gauge <- (points$stage_upper_cm[1] - points$stage_lower_cm[1]) / 2
  counts <- table(droplevels(points$flume))
  sample_text <- if (length(unique(counts)) == 1L) {
    sprintf("%d observations per flume.", counts[1])
  } else {
    paste0(paste(sprintf("%s: %d observations", names(counts), counts), collapse = "; "), ".")
  }
  point_text <- paste(
    "White-filled circles represent posterior medians of observation-specific latent true flow;",
    "smaller filled circles represent measured flow and are drawn above the latent estimates.")
  interval_text <- sprintf(paste(
    "Vertical bars show %s equal-tailed parameter credible intervals for mean rating-curve flow at each observed stage.",
    "These bars describe the mean rating curve, rather than uncertainty in the observation-specific latent points;",
    "the latent points need not coincide with the centers of these intervals."), prob)
  gauge_text <- sprintf(paste(
    "Horizontal bars show point-gage accuracy of +/- %g cm, without assigning probability coverage.",
    "Stage is fixed in the fitted model, so this accuracy is displayed separately and is not propagated into the vertical intervals."), gauge)
  relative_interval_text <- sub("for mean rating-curve flow", "for the percent relative difference in mean rating-curve flow", interval_text, fixed = TRUE)
  curves_text <- sprintf(paste(
    "Colored curves show posterior medians of the mean stage-discharge relationship; shaded bands show %s",
    "equal-tailed parameter credible intervals. %s %s The dotted black curve is the manufacturer reference.",
    "Intervals condition on observed stage and do not add predictive process or measurement scatter."),
    prob, point_text, interval_text)
  writeLines(c(
    "# Figure captions", "",
    "## Stage-discharge calibration points", "",
    paste("Stage-discharge calibration of two 3D-printed flumes and one fiberglass flume;", sample_text,
          point_text, "The solid black curve shows the manufacturer equation.", interval_text, gauge_text,
          "Panels share axis limits; the reference curve extends from zero to 10% above the maximum observed stage."), "",
    "## Relative differences from the manufacturer", "",
    paste("Percent relative differences between calibration flow and the manufacturer prediction at each observed stage,",
          "calculated as 100 x (flow - manufacturer prediction) / manufacturer prediction.",
          point_text, "The solid black zero line indicates agreement with the manufacturer prediction; negative values indicate lower flow.",
          relative_interval_text, gauge_text,
          "Manufacturer predictions are treated as fixed. Relative differences at zero manufacturer flow are undefined and omitted.",
          "Panels share axis limits."), "",
    "## Combined stage-discharge curves", "", curves_text, "",
    "## Stage-discharge curves by flume", "",
    paste(curves_text, "Separate panels share axis limits."), "",
    "## Equation and units", "",
    sprintf("The manufacturer equation is Q (cfs) = %g H (ft)^%g. Figures report stage in cm and flow in L/s.",
            manufacturer$a_m, manufacturer$b_m)
  ), file.path(output_dir, paste0(prefix, "_flume_figure_captions.md")))
}
