# Sourceable companion plots for a new observed discharge at fixed stage.
# Sourcing defines functions only; it does not read a fit, refit, or export files.

build_flume_uncertainty_components <- function(post, dat, manufacturer,
    flume_names = c("1" = "3D Printed 1", "2" = "3D Printed 2", "3" = "Fiberglass"),
    grid_points = 160L) {
  assert_flume_metric_units(dat)
  required <- c("a", "b", "sigma_proc", "sigma_Q_add", "cv_Q")
  stopifnot(all(required %in% names(post)),
            all(c("H_obs", "F", "J") %in% names(dat)),
            is.matrix(post$a), nrow(post$a) > 1L, ncol(post$a) == dat$J,
            identical(dim(post$a), dim(post$b)),
            identical(dim(post$a), dim(post$sigma_proc)),
            identical(dim(post$a), dim(post$sigma_Q_add)),
            length(post$cv_Q) == nrow(post$a),
            length(dat$H_obs) == length(dat$F),
            all(is.finite(dat$H_obs)), all(dat$H_obs > 0),
            all(dat$F %in% seq_len(dat$J)),
            all(as.character(seq_len(dat$J)) %in% names(flume_names)),
            length(grid_points) == 1L, grid_points >= 2L,
            nrow(manufacturer) == 1L,
            all(c("a_m", "b_m") %in% names(manufacturer)),
            is.finite(manufacturer$a_m), manufacturer$a_m > 0,
            is.finite(manufacturer$b_m), manufacturer$b_m > 0)
  stopifnot(all(vapply(post[required], function(x) all(is.finite(x)) && all(x > 0), logical(1))))
  parts <- c("Infiltration", "Flowmeter", "Process", "Rating curve")
  rows <- lapply(seq_len(dat$J), function(j) {
    observed_h <- dat$H_obs[dat$F == j]
    stopifnot(length(unique(observed_h)) >= 2L)
    h <- sort(unique(c(seq(min(observed_h), max(observed_h), length.out = grid_points), observed_h)))
    mu <- vapply(h, function(h_i) post$a[, j] * h_i^post$b[, j], numeric(nrow(post$a)))
    proc2 <- post$sigma_proc[, j]^2
    infil2 <- post$sigma_Q_add[, j]^2
    cv2 <- as.numeric(post$cv_Q)^2
    # Exact moments for the implemented likelihood, integrating a NEW latent Q:
    # E(Q_new|theta,H) = mu;
    # Var(Q_new|theta,H) = proc2 + infil2 + cv2*(mu^2 + proc2).
    # Law of total variance adds Var_posterior(mu).
    curve_var <- apply(mu, 2, stats::var)
    meter_var <- colMeans(sweep(sweep(mu^2, 1, proc2, "+"), 1, cv2, "*"))
    v <- cbind("Infiltration" = rep(mean(infil2), length(h)),
               "Flowmeter" = meter_var, "Process" = rep(mean(proc2), length(h)),
               "Rating curve" = curve_var)
    total <- rowSums(v)
    # Check the sum independently, retaining the cv2*proc2 interaction.
    expected_conditional <- colMeans(sweep(mu^2, 1, cv2, "*")) +
                            mean(proc2 + infil2 + cv2 * proc2)
    stopifnot(all(is.finite(v)), all(v >= 0), all(total > 0),
              max(abs(total - (curve_var + expected_conditional))) < 1e-10)
    ref <- manufacturer$a_m * h^manufacturer$b_m
    do.call(rbind, lapply(seq_along(parts), function(k) {
      data.frame(flume = unname(flume_names[as.character(j)]), stage_m = h,
                 reference_flow_lps = ref, mean_flow_lps = colMeans(mu),
                 component = parts[k], variance_lps2 = v[, k],
                 combined_variance_lps2 = total,
                 sd_lps = sqrt(v[, k]), relative_sd_pct = 100 * sqrt(v[, k]) / ref,
                 variance_share_pct = 100 * v[, k] / total,
                 posterior_draws = nrow(post$a))
    }))
  })
  budget <- do.call(rbind, rows)
  budget$flume <- factor(budget$flume, levels = unname(flume_names[as.character(seq_len(dat$J))]))
  budget$component <- factor(budget$component, levels = parts)
  shares <- stats::aggregate(variance_share_pct ~ flume + stage_m, budget, sum)
  stopifnot(max(abs(shares$variance_share_pct - 100)) < 1e-10)
  budget
}

plot_flume_uncertainty_components <- function(budget) {
  parts <- levels(budget$component)
  colors <- c("Infiltration" = "#CC79A7", "Flowmeter" = "#E69F00",
              "Process" = "#009E73", "Rating curve" = "#0072B2", "Combined" = "#111111")
  total_rows <- budget[budget$component == "Infiltration", ]
  total_rows$component <- "Combined"
  total_rows$sd_lps <- sqrt(total_rows$combined_variance_lps2)
  total_rows$relative_sd_pct <- 100 * total_rows$sd_lps / total_rows$reference_flow_lps
  lines <- rbind(transform(budget, component = as.character(component)), total_rows)
  lines$component <- factor(lines$component, levels = c(parts, "Combined"))
  display <- rbind(transform(lines, metric = "Absolute SD (L/s)", uncertainty = sd_lps),
                   transform(lines, metric = "Relative SD (%)", uncertainty = relative_sd_pct))
  display$metric <- factor(display$metric, levels = c("Absolute SD (L/s)", "Relative SD (%)"))
  plot_theme <- ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), legend.position = "bottom",
                   legend.title = ggplot2::element_blank(),
                   strip.background = ggplot2::element_rect(fill = "grey95"),
                   panel.spacing = grid::unit(0.8, "lines"))
  p_sd <- ggplot2::ggplot(display, ggplot2::aes(stage_m, uncertainty, color = component, linetype = component)) +
    ggplot2::geom_line(linewidth = 1.05) +
    ggplot2::facet_grid(metric ~ flume, scales = "free_y") +
    ggplot2::scale_color_manual(values = colors, drop = FALSE) +
    ggplot2::scale_linetype_manual(values = c("Infiltration" = "solid", "Flowmeter" = "longdash",
                                            "Process" = "dotted", "Rating curve" = "dotdash", "Combined" = "solid")) +
    ggplot2::scale_x_continuous(breaks = c(0.06, 0.09, 0.12, 0.15)) +
    ggplot2::scale_y_continuous(limits = c(0, NA), expand = ggplot2::expansion(mult = c(0, 0.06))) +
    ggplot2::labs(x = "Stage (m)", y = "Uncertainty scale (one standard deviation)") + plot_theme
  p_share <- ggplot2::ggplot(budget, ggplot2::aes(stage_m, variance_share_pct, fill = component)) +
    ggplot2::geom_area(position = ggplot2::position_stack(reverse = TRUE), linewidth = 0) +
    ggplot2::facet_wrap(ggplot2::vars(flume), nrow = 1) +
    ggplot2::scale_fill_manual(values = colors[parts], drop = FALSE) +
    ggplot2::scale_x_continuous(breaks = c(0.06, 0.09, 0.12, 0.15)) +
    # Clip display rather than dropping vertices just above 100 due to roundoff.
    ggplot2::scale_y_continuous(breaks = seq(0, 100, 25), expand = c(0, 0)) +
    ggplot2::coord_cartesian(ylim = c(0, 100), expand = FALSE) +
    ggplot2::labs(x = "Stage (m)", y = "Share of total predictive variance (%)") + plot_theme
  stopifnot(is.null(p_sd$labels$caption), is.null(p_share$labels$caption),
            is.null(p_sd$labels$title), is.null(p_share$labels$title))
  invisible(ggplot2::ggplot_build(p_sd))
  share_build <- ggplot2::ggplot_build(p_share)
  stopifnot(all(is.finite(share_build$data[[1]]$ymin)), all(is.finite(share_build$data[[1]]$ymax)))
  list(sd = p_sd, shares = p_share)
}

save_flume_uncertainty_components <- function(post, dat, manufacturer, output_dir,
    prefix = basename(output_dir), input_csv = NULL, fit_file = NULL,
    flume_names = c("1" = "3D Printed 1", "2" = "3D Printed 2", "3" = "Fiberglass")) {
  budget <- build_flume_uncertainty_components(post, dat, manufacturer, flume_names)
  plots <- plot_flume_uncertainty_components(budget)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  path <- function(stem, ext) file.path(output_dir, paste0(prefix, "_flume_uncertainty_", stem, ".", ext))
  for (ext in c("png", "jpg")) {
    ggplot2::ggsave(path("components_sd", ext), plots$sd, width = 11, height = 7, dpi = 300, bg = "white")
    ggplot2::ggsave(path("variance_shares", ext), plots$shares, width = 11, height = 4.4, dpi = 300, bg = "white")
  }
  utils::write.csv(budget, path("components_data", "csv"), row.names = FALSE)
  captions <- c(
    "# Uncertainty figure captions", "",
    "## Absolute and relative uncertainty components", "",
    "Components of predictive uncertainty for a new discharge measurement at fixed stage in two 3D-printed flumes and one fiberglass flume. Colored lines show the square roots of posterior variance contributions associated with infiltration, proportional flowmeter uncertainty, process variation around the mean rating curve, and uncertainty in the rating-curve parameters. The black line shows the combined predictive standard deviation, obtained from the sum of the four variance contributions. Upper panels report L/s; lower panels express these scales as percentages of the manufacturer prediction at the same stage. Lines represent one-standard-deviation scales, not 95% interval half-widths, and component standard deviations do not add linearly. Each panel spans the observed stage range for that specimen. Stage and the manufacturer equation are treated as fixed.",
    "",
    "## Shares of predictive variance", "",
    "Relative contributions of infiltration uncertainty, proportional flowmeter uncertainty, process variation, and rating-curve parameter uncertainty to predictive discharge variance at fixed stage. Stacked areas show each component as a percentage of total predictive variance and sum to 100% at each stage. A decreasing share indicates a smaller contribution relative to the other components and does not necessarily indicate a decrease in absolute variance. The flowmeter contribution includes its interaction with process variation, arising because measurement uncertainty depends on latent true discharge. Contributions are calculated from joint posterior draws under the specified calibration model. Each panel spans the observed stage range for that specimen; stage uncertainty is not propagated.")
  writeLines(captions, path("figure_captions", "md"))
  provenance <- function(label, file) {
    if (is.null(file)) return(character())
    sprintf("%s: `%s`; MD5 `%s`.", label, normalizePath(file, winslash = "/"), unname(tools::md5sum(file)))
  }
  diag_file <- file.path(output_dir, paste0(prefix, "_flume_model_3_chain_diagnostics.csv"))
  diagnostic_note <- "Sampling diagnostics were not supplied; consult the fitted model before interpreting these allocations."
  if (file.exists(diag_file)) {
    diagnostics <- utils::read.csv(diag_file)
    diagnostic_note <- sprintf("Current fit: %d divergences among %d draws; E-BFMI %.3f to %.3f. Component estimates are provisional until sampling diagnostics are resolved.",
       sum(diagnostics$divergent_transitions), sum(diagnostics$post_warmup_draws), min(diagnostics$ebfmi), max(diagnostics$ebfmi))
  }
  notes <- c(
    "# Predictive uncertainty component definitions", "",
    "These companion figures describe a NEW observed discharge at known, fixed stage.",
    "They do not decompose the existing 95% parameter credible bars for mean rating-curve flow.",
    "Infiltration retains the parameter's terminology; its SD does not estimate a directional loss.",
    "", "## Calculation", "",
    "Let theta denote one joint posterior draw and mu = a*H^b, in native m and L/s units.",
    "V_curve = Var_posterior(mu).",
    "V_process = E_posterior(sigma_proc^2).",
    "V_infiltration = E_posterior(sigma_Q_add^2).",
    "V_flowmeter = E_posterior[cv_Q^2 * (mu^2 + sigma_proc^2)].",
    "V_total = V_curve + V_process + V_infiltration + V_flowmeter.",
    "The law of total variance integrates over a new latent Q_true; no random simulation is needed.",
    "The cv_Q^2*sigma_proc^2 interaction is allocated to the flowmeter component, where it originates in the likelihood.",
    "Expectations use paired joint draws; products of marginal summaries are not used.",
    "Noise-scale parameter uncertainty is integrated through these expectations, not added as another independent variance term.",
    "Variances are already in (L/s)^2; each component line is sqrt(V_component) and the combined line is sqrt(V_total).",
    "Relative SD = 100*SD/manufacturer prediction at the same stage. These are not relative biases or 95% intervals.",
    "Stacked shares = 100*V_component/V_total. Shares sum to 100% at every stage.",
    if (!is.null(post$H_true))
      "The fitted coefficient draws integrate calibration-stage uncertainty. These predictions condition on a new known stage and do not add new stage-measurement uncertainty or manufacturer uncertainty; panels remain within each specimen's observed stage range."
    else
      "Neither stage nor manufacturer uncertainty is propagated; panels remain within each specimen's observed stage range.",
    "The allocations describe the fitted model; physical sources are not isolated experimentally.",
    if (!is.null(dat$cv_Q_prior_median))
      sprintf("The proportional measurement-SD coefficient is estimated; its prior median is %g%%. This empirical SD model does not explicitly reproduce the manufacturer's additive velocity-error specification.", dat$cv_Q_prior_median * 100)
    else
      "The proportional measurement-SD coefficient is estimated. Consult run metadata for its prior median; it is not a fixed measurement SD.",
    sprintf("Posterior mean cv_Q: %.4f%%; root mean square cv_Q: %.4f%%.", mean(post$cv_Q)*100, sqrt(mean(post$cv_Q^2))*100),
    "The flowmeter line also depends on the fitted mean, its uncertainty, and process scatter, so it need not equal either coefficient summary.",
    "", "## Provenance and diagnostics", "",
    provenance("Calibration input", input_csv), provenance("Saved fit", fit_file),
    sprintf("%d joint posterior draws; %d flume-stage combinations.", nrow(post$a), nrow(budget)/4),
    diagnostic_note,
    "Figure captions are exported separately in the companion uncertainty_figure_captions.md file; images contain no captions.")
  writeLines(notes, path("components_notes", "md"))
  invisible(list(sd = plots$sd, shares = plots$shares, data = budget,
                 captions = path("figure_captions", "md"), data_file = path("components_data", "csv")))
}
