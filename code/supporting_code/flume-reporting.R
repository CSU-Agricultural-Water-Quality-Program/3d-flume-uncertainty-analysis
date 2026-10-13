# Supporting reporting functions; sourced by code/bayes_flume.R.
# Defines summaries, diagnostics, rating-curve plots, and manuscript tables.
# Loading this file does not fit a model or export files.

# metric_precis
metric_precis <- function(fit, prob = 0.89) {
  original <- precis(fit, depth = 2, prob = prob)
  if (!is.null(attr(fit, "cstanfit"))) {
    draws <- attr(fit, "cstanfit")$draws(variables=rownames(original), format="draws_array")
  } else {
    draws <- posterior::as_draws_array(as.array(attr(fit, "stanfit")))
    draws <- draws[, , rownames(original), drop=FALSE]
  }
  names <- dimnames(draws)[[3]]
  names <- sub("^a\\[", "a_m[", names)
  names <- sub("^alpha\\[", "alpha_m[", names)
  names <- sub("^H_true\\[", "H_true_m[", names)
  names <- sub("^(Q_true|mu|sigma_proc|sigma_Q_add)(\\[|$)", "\\1_lps\\2", names)
  dimnames(draws)[[3]] <- names
  result <- as.data.frame(posterior::summarise_draws(
    draws, mean=base::mean, sd=stats::sd,
    interval=function(x) stats::quantile(x, probs=c((1-prob)/2,(1+prob)/2)),
    rhat=posterior::rhat, ess_bulk=posterior::ess_bulk))
  rownames(result) <- result$variable
  result$variable <- NULL
  methods::new("precis", result, digits=3)
}

# write_metric_precis
write_metric_precis <- function(summary, filename) {
  result <- as.data.frame(summary)
  v <- rownames(result)
  result$units <- ifelse(grepl("^H_true_m\\[",v), "m",
    ifelse(grepl("_lps",v), "L/s",
    ifelse(grepl("^a_m\\[",v), "(L/s)/m^b",
    ifelse(grepl("^alpha_m\\[",v), "log(a_m)", "dimensionless"))))
  write.csv(result,filename,row.names=TRUE)
}

# dashboard_flume
dashboard_flume <- function(fit) {
  backend <- attr(fit, "cstanfit")
  if (!is.null(backend)) {
    proxy <- as.list(backend)
    class(proxy) <- class(backend)
    proxy$sampler_diagnostics <- function(...) {
      x <- backend$sampler_diagnostics(...)
      required <- c("accept_stat__", "divergent__", "energy__")
      variables <- dimnames(x)[[3]]
      stopifnot(all(required %in% variables))
      x[, , c(required, setdiff(variables, required)), drop = FALSE]
    }
    attr(fit, "cstanfit") <- proxy
  }
  dashboard(fit)
}

# save_dashboard
save_dashboard <- function(fit, filename) {
  png(filename, width = 8, height = 8, units = "in", res = 300)
  on.exit(dev.off(), add = TRUE)
  dashboard_flume(fit)
  invisible(filename)
}

# export_stan_diagnostics
export_stan_diagnostics <- function(fit, output_dir, prefix, dat) {
  backend <- attr(fit, "cstanfit")
  stopifnot(!is.null(backend))
  x <- backend$sampler_diagnostics(inc_warmup = FALSE)
  lp <- backend$draws(variables = "lp__", format = "draws_array")
  metadata <- backend$metadata()
  required <- c("energy__", "divergent__", "treedepth__", "stepsize__",
                "accept_stat__", "n_leapfrog__")
  stopifnot(all(required %in% dimnames(x)[[3]]),
            identical(dim(x)[1:2], dim(lp)[1:2]))
  chain_ids <- metadata$id
  n_iter <- dim(x)[1]
  chain_rows <- lapply(seq_len(dim(x)[2]), function(j) {
    energy <- x[, j, "energy__"]
    # Match CmdStanR's finite-sample normalization (N energy draws).
    bfmi <- sum(diff(energy)^2) / (length(energy) * stats::var(energy))
    data.frame(
      chain = chain_ids[j], post_warmup_draws = n_iter,
      energy_mean = mean(energy), energy_sd = stats::sd(energy),
      energy_median = median(energy), energy_q025 = unname(quantile(energy, 0.025)),
      energy_q975 = unname(quantile(energy, 0.975)), ebfmi = bfmi,
      ebfmi_below_0_3 = bfmi < 0.3,
      divergent_transitions = sum(x[, j, "divergent__"]),
      divergent_fraction = mean(x[, j, "divergent__"]),
      max_treedepth = metadata$max_treedepth,
      max_treedepth_hits = sum(x[, j, "treedepth__"] >= metadata$max_treedepth),
      treedepth_mean = mean(x[, j, "treedepth__"]),
      acceptance_mean = mean(x[, j, "accept_stat__"]),
      step_size_mean = mean(x[, j, "stepsize__"]),
      leapfrog_steps_mean = mean(x[, j, "n_leapfrog__"])
    )
  })
  chain_summary <- do.call(rbind, chain_rows)
  # Validate named-column calculations against CmdStanR's own diagnostics.
  reference <- backend$diagnostic_summary(quiet = TRUE)
  stopifnot(identical(as.numeric(chain_summary$divergent_transitions),
                      as.numeric(reference$num_divergent)),
            identical(as.numeric(chain_summary$max_treedepth_hits),
                      as.numeric(reference$num_max_treedepth)),
            isTRUE(all.equal(chain_summary$ebfmi, reference$ebfmi, tolerance = 1e-10)))
  draw_rows <- lapply(seq_len(dim(x)[2]), function(j) {
    variables <- dimnames(x)[[3]]
    values <- as.data.frame(setNames(lapply(variables, function(v) as.numeric(x[, j, v])),
                                    variables), check.names = FALSE)
    data.frame(chain = chain_ids[j], iteration = seq_len(n_iter),
               phase = "sampling", lp__ = as.numeric(lp[, j, "lp__"]), values,
               check.names = FALSE)
  })
  write.csv(chain_summary, file.path(output_dir, paste0(prefix, "_chain_diagnostics.csv")),
            row.names = FALSE)
  write.csv(do.call(rbind, draw_rows),
            file.path(output_dir, paste0(prefix, "_sampler_diagnostics.csv")), row.names = FALSE)
  run_metadata <- data.frame(
    field = c("input_csv", "input_md5", "calibration_stage_units", "calibration_flow_units",
              "reporting_stage_units", "reporting_flow_units", "seed", "chains",
              "warmup_per_chain", "sampling_per_chain", "adapt_delta", "max_treedepth",
              "run_started_utc", "rethinking_version", "cmdstanr_version", "cmdstan_version",
              "stage_accuracy_half_width_m", "stage_error_model",
              "stage_accuracy_scope", "cv_Q_prior_median_fraction",
              "cv_Q_prior_log_sd", "sigma_Q_add_prior_median_lps",
              "measurement_error_interpretation", "a_prior_median", "b_prior_median",
              "rating_prior_log_sd", "sigma_proc_prior_rate_lps",
              "sigma_Q_add_prior_log_sd", "manufacturer_equation"),
    value = c(normalizePath(input_csv, winslash = "/"), unname(tools::md5sum(input_csv)),
              "m", "L/s", "m", "L/s", metadata$seed, dim(x)[2],
              metadata$iter_warmup, metadata$iter_sampling, metadata$adapt_delta,
              metadata$max_treedepth, metadata$start_datetime[1],
              as.character(packageVersion("rethinking")), as.character(packageVersion("cmdstanr")),
              as.character(cmdstanr::cmdstan_version()), dat$sigma_H,
              "H_true | H_obs ~ Uniform(H_obs - bound, H_obs + bound)",
              "reported accuracy applied to computed depth; independent across observations",
              dat$cv_Q_prior_median, dat$cv_Q_prior_log_sd, dat$sigma_Q_add_prior_median,
              "empirical proportional plus additive SD; manufacturer velocity term not explicit",
              dat$a_prior_median, dat$b_prior_median, dat$rating_prior_log_sd,
              dat$sigma_proc_prior_rate, dat$sigma_Q_add_prior_log_sd,
              "Q (L/s) = 941.2 * H (m)^2.58")
  )
  write.csv(run_metadata, file.path(output_dir, paste0(prefix, "_run_metadata.csv")),
            row.names = FALSE)
  invisible(chain_summary)
}

# plot_flume_posterior
plot_flume_posterior <- function(post, dat, prob = 0.95, ndraws = 4000,
                                 show_manufacturer = NULL,
                                 ribbon_kind = c("process", "parameter", "predictive"),
                                 show_qobs_errors = FALSE,
                                 extrapolation_fraction = 0.10,
                                 flume_names = c("1" = "3D Printed 1",
                                                 "2" = "3D Printed 2",
                                                 "3" = "Fiberglass")) {
  
  assert_flume_metric_units(dat)
  ribbon_kind <- match.arg(ribbon_kind)
  stopifnot(all(c("Q_obs","H_obs","F","J") %in% names(dat)))
  
  Q    <- as.numeric(dat$Q_obs)
  H    <- as.numeric(dat$H_obs)
  Fidx <- as.integer(dat$F)
  J    <- as.integer(dat$J)
  
  stopifnot(all(is.finite(H)), all(H >= 0), max(H) > 0,
            length(extrapolation_fraction) == 1L,
            is.finite(extrapolation_fraction), extrapolation_fraction >= 0,
            all(is.finite(post$a)), all(post$a > 0),
            all(is.finite(post$b)), all(post$b > 0))
  # Evaluate in calibration units, including exact zero. For positive a and b,
  # exp(log(a) + b * log(0)) = 0, the exact power-law limit at zero stage.
  Hgrid <- seq(0, max(H) * (1 + extrapolation_fraction), length.out = 200)
  logH  <- log(Hgrid)
  
  # subset posterior draws
  Sfull <- nrow(post$a)
  S     <- min(ndraws, Sfull)
  sel   <- if (S < Sfull) sample.int(Sfull, S) else seq_len(S)
  
  a_draw  <- post$a[sel, , drop=FALSE]
  b_draw  <- post$b[sel, , drop=FALSE]
  sp_draw <- post$sigma_proc[sel, , drop=FALSE]
  cv_draw <- as.numeric(post$cv_Q[sel])
  sq_draw <- post$sigma_Q_add[sel, , drop=FALSE]
  
  flumes_present <- sort(unique(Fidx))
  
  # central interval endpoints and z for prob
  p_lo  <- (1 - prob)/2
  p_hi  <- 1 - p_lo
  zcrit <- qnorm((1 + prob)/2)
  
  # ---------- Smooth ribbons for Q_true ----------
  # No per-point random noise; everything is built from smooth functions
  rib_list <- lapply(flumes_present, function(j){
    # posterior of mu(H) on the grid
    log_mu <- outer(b_draw[, j], logH, `*`) +
      matrix(log(a_draw[, j]), nrow=S, ncol=length(Hgrid))
    mu_mat <- exp(log_mu)
    
    mu_med <- apply(mu_mat, 2, median)
    mu_lo  <- apply(mu_mat, 2, quantile, p_lo)
    mu_hi  <- apply(mu_mat, 2, quantile, p_hi)
    
    if (ribbon_kind == "parameter") {
      Q_lo <- mu_lo
      Q_hi <- mu_hi
    } else if (ribbon_kind == "process") {
      # add smooth process width using the median sigma_proc
      sp_med <- median(sp_draw[, j])
      Q_lo <- mu_med - zcrit * sp_med
      Q_hi <- mu_med + zcrit * sp_med
    } else { # predictive: include measurement + process on top of mu_med
      cv_med  <- median(cv_draw)
      sq_med  <- median(sq_draw[, j])
      sp_med  <- median(sp_draw[, j])
      sd_pred <- sqrt((cv_med * mu_med)^2 + sq_med^2 + sp_med^2)
      Q_lo <- mu_med - zcrit * sd_pred
      Q_hi <- mu_med + zcrit * sd_pred
    }
    
    data.frame(
      flume = j,
      H     = Hgrid,
      Q_lo  = Q_lo,
      Q_hi  = Q_hi,
      Q_med = mu_med
    )
  })
  ribbon_df <- do.call(rbind, rib_list)
  
  # White points always summarize observation-specific latent Q_true. Parameter
  # bars summarize mean rating-curve flow, matching the parameter ribbons;
  # these two estimands are distinguished in the exported manuscript captions.
  Qinterval_draws <- if (ribbon_kind == "parameter") {
    vapply(seq_along(H), function(i) {
      a_draw[, Fidx[i]] * H[i]^b_draw[, Fidx[i]]
    }, numeric(S))
  } else {
    post$Q_true[sel, , drop = FALSE]
  }
  Qtrue_df <- data.frame(
    H    = H,
    Q    = apply(post$Q_true[sel, , drop = FALSE], 2, median),
    Q_parameter_med = if (ribbon_kind == "parameter") apply(Qinterval_draws, 2, median) else NA_real_,
    Q_lo = apply(Qinterval_draws, 2, quantile, p_lo),
    Q_hi = apply(Qinterval_draws, 2, quantile, p_hi),
    flume = Fidx,
    ptype = "Q_true"
  )
  
  # ---------- Q_obs points (optional error bars controlled by prob) ----------
  cv_med   <- median(post$cv_Q)
  sq_med_j <- apply(post$sigma_Q_add, 2, median)
  sp_med_j <- apply(post$sigma_proc,   2, median)
  se_obs   <- sqrt( (cv_med * Q)^2 + (sq_med_j[Fidx])^2 + (sp_med_j[Fidx])^2 )
  Qobs_df  <- data.frame(
    H = H, Q = Q,
    Q_lo = pmax(Q - zcrit*se_obs, 0),
    Q_hi = Q + zcrit*se_obs,
    flume = Fidx,
    ptype = "Q_obs"
  )
  
  # ---------- Manufacturer curve ----------
  manuf_df <- NULL
  if (!is.null(show_manufacturer)) {
    if ("flume" %in% names(show_manufacturer)) {
      tmp <- subset(show_manufacturer, flume %in% flumes_present)
      manuf_df <- do.call(rbind, lapply(tmp$flume, function(j){
        a_m <- tmp$a_m[tmp$flume == j][1]; b_m <- tmp$b_m[tmp$flume == j][1]
        data.frame(flume = j, H = Hgrid, Qm = a_m * Hgrid^b_m)
      }))
    } else {
      manuf_df <- data.frame(H = Hgrid,
                             Qm = show_manufacturer$a_m[1] * Hgrid^(show_manufacturer$b_m[1]))
    }
  }
  
  # Inputs, posterior draws, intervals, and manufacturer flow are already metric.
  # Use coordinate limits so negative Gaussian interval tails are clipped only
  # in the display, rather than deleting rows or changing posterior summaries.
  y_max <- max(ribbon_df$Q_hi, ribbon_df$Q_med, Qtrue_df$Q_hi, Qtrue_df$Q, Qobs_df$Q,
               if (show_qobs_errors) Qobs_df$Q_hi,
               if (!is.null(manuf_df)) manuf_df$Qm) * 1.10

  # ---------- aesthetics ----------
  color_map <- c("1"="#0072B2","2"="#D55E00","3"="#009E73")
  color_map <- color_map[intersect(names(color_map), as.character(flumes_present))]
  line_map <- c("1"="solid", "2"="longdash", "3"="dotdash")[names(color_map)]
  curve_labels <- flume_names[names(color_map)]
  if (!is.null(manuf_df)) {
    color_map <- c(color_map, "Manufacturer equation"="black")
    line_map <- c(line_map, "Manufacturer equation"="dotted")
    curve_labels <- c(curve_labels, "Manufacturer equation")
  }
  # White-filled estimated-flow circles contrast with colored measured-flow
  # circles. Thin intervals and light ribbons reduce visual crowding.
  shape_map <- c("Q_true" = 21, "Q_obs" = 21)
  

  ggplot() +
    # ribbons and median curve
    geom_ribbon(data = ribbon_df,
                aes(H, ymin = Q_lo, ymax = Q_hi, fill = factor(flume),
                    color = factor(flume)),
                alpha = 0.12, linewidth = 0.3, show.legend = FALSE) +
    geom_line(data = ribbon_df,
              aes(H, Q_med, color = factor(flume), linetype = factor(flume)),
              linewidth = 1.0) +
    
    # Estimated flow: white fill with flume-colored outline.
    geom_linerange(data = Qtrue_df,
                   aes(H, ymin = Q_lo, ymax = Q_hi, color = factor(flume)),
                   linewidth = 0.65, alpha = 0.7, show.legend = FALSE) +
    
    # Q_obs: colored fill with black outline.
    { if (show_qobs_errors)
      geom_linerange(data = Qobs_df,
                     aes(H, ymin = Q_lo, ymax = Q_hi, color = factor(flume)),
                     linewidth = 0.55, alpha = 0.55, show.legend = FALSE) } +
    
    # manufacturer equation
    { if (!is.null(manuf_df))
      geom_line(data = manuf_df,
                aes(H, Qm, color = "Manufacturer equation",
                    linetype = "Manufacturer equation"),
                linewidth = 1.8) } +
    
    # Larger white latent estimates beneath smaller measured-flow points.
    geom_point(data = Qtrue_df,
               aes(H, Q, color = factor(flume), shape = ptype),
               size = 2.5, stroke = 0.8, fill = "white",
               show.legend = c(color = FALSE, shape = TRUE)) +
    geom_point(data = Qobs_df,
               aes(H, Q, fill = factor(flume), shape = ptype),
               color = "black", size = 1.5, stroke = 0.45,
               show.legend = c(fill = FALSE, shape = TRUE)) +

    # scales and guides
    scale_color_manual(values = color_map,
                       name = "Flume",
                       breaks = names(color_map), labels = curve_labels) +
    scale_fill_manual(values = color_map, guide = "none") +
    scale_shape_manual(values = shape_map,
                       name = "Points",
                       breaks = c("Q_true","Q_obs"),
                       labels = c("Estimated latent true flow", "Measured flow")) +
    scale_linetype_manual(values = line_map, name = "Flume",
                          breaks = names(line_map), labels = curve_labels) +
    guides(
      color   = guide_legend(order = 1, nrow = 1),
      linetype= guide_legend(order = 1, nrow = 1),
      shape   = guide_legend(
        order = 2, nrow = 1,
        override.aes = list(
          # two keys, in order of breaks: Q_true then Q_obs
          shape = c(21, 21),
          size  = c(2.5, 2.2),
          stroke= c(0.45, 0.8),
          color = c("black", "black"),  # outline in legend
          fill  = c("white", "grey60")  # generic estimated versus measured flow
        )
      )
    ) +
    scale_x_continuous(expand = expansion(mult = 0)) +
    scale_y_continuous(expand = expansion(mult = 0)) +
    coord_cartesian(xlim = c(0, max(Hgrid)),
                    ylim = c(0, y_max), expand = FALSE) +
    labs(x = "Stage (m)", y = "Flow rate (L/s)") +
    theme_bw(base_size = 13) +
    theme(legend.position = "bottom", legend.box = "vertical",
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 11),
          legend.key.width = grid::unit(1.4, "cm"),
          panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
          panel.grid.minor = element_blank())
}

# generate_flume_table
generate_flume_table <- function(post, depths_m = seq.int(0L, 24L) / 100,
                                 prob = 0.95,
                                 kind = c("parameter", "process", "predictive"),
                                 ndraws = 4000,
                                 manuf = NULL,
                                 flume_names = c("1"="3D Printed 1","2"="3D Printed 2","3"="Fiberglass"),
                                 digits = 3) {
  
  kind <- match.arg(kind)
  stopifnot(is.matrix(post$a), is.matrix(post$b))
  
  # subset draws
  Sfull <- nrow(post$a)
  S     <- min(ndraws, Sfull)
  sel   <- if (S < Sfull) sample.int(Sfull, S) else seq_len(S)
  
  a  <- post$a[sel, , drop=FALSE]              # S x J
  b  <- post$b[sel, , drop=FALSE]              # S x J
  sp <- post$sigma_proc[sel, , drop=FALSE]     # S x J
  cv <- as.numeric(post$cv_Q[sel])             # S
  sq <- post$sigma_Q_add[sel, , drop=FALSE]    # S x J
  
  J <- ncol(a)
  stopifnot(length(depths_m) > 0, all(is.finite(depths_m)), all(depths_m >= 0),
            all(is.finite(a)), all(a > 0), all(is.finite(b)), all(b > 0))
  H <- depths_m  # evaluate the metric calibration directly
  G <- length(H)
  p_lo <- (1 - prob)/2
  p_hi <- 1 - p_lo
  
  # compute per-flume summaries
  out_list <- lapply(seq_len(J), function(j){
    # S x G mu draws
    mu_mat <- exp( outer(b[,j], log(H), `*`) + matrix(log(a[,j]), nrow=S, ncol=G) )
    
    if (kind == "parameter") {
      lo <- apply(mu_mat, 2, quantile, p_lo)
      md <- apply(mu_mat, 2, median)
      hi <- apply(mu_mat, 2, quantile, p_hi)
    } else if (kind == "process") {
      # add process noise with draw-specific sigma_proc (smooth summary via simulation)
      yrep <- mu_mat + matrix(rnorm(S*G), nrow=S, ncol=G) * matrix(sp[,j], nrow=S, ncol=G)
      lo <- apply(yrep, 2, quantile, p_lo)
      md <- apply(mu_mat, 2, median)  # center is mu
      hi <- apply(yrep, 2, quantile, p_hi)
    } else { # predictive: add process + measurement
      sd_pred <- sqrt( (cv^2) %o% (rep(1,G)) * (mu_mat^2) + (sq[,j]^2) %o% rep(1,G) + (sp[,j]^2) %o% rep(1,G) )
      yrep <- mu_mat + matrix(rnorm(S*G), nrow=S, ncol=G) * sd_pred
      lo <- apply(yrep, 2, quantile, p_lo)
      md <- apply(mu_mat, 2, median)
      hi <- apply(yrep, 2, quantile, p_hi)
    }
    
    data.frame(Depth_m = depths_m,
               flume = names(flume_names)[j],
               median_lps = md,
               lo_lps = lo, hi_lps = hi)
  })
  
  df_long <- do.call(rbind, out_list)
  
  # wide table with one set of columns per flume
  # value columns as "median [lo, hi]" strings for quick reporting
  library(dplyr); library(tidyr)
  df_fmt <- df_long %>%
    mutate(label = paste0(sprintf(paste0("%.",digits,"f"), median_lps),
                          " [", sprintf(paste0("%.",digits,"f"), lo_lps),
                          ", ", sprintf(paste0("%.",digits,"f"), hi_lps), "]"),
           Flume = paste0(flume_names[flume], " (L/s)")) %>%
    select(Depth_m, Flume, label) %>%
    pivot_wider(names_from = Flume, values_from = label) %>%
    arrange(Depth_m)
  
  # manufacturer column
  if (!is.null(manuf)) {
    if (all(c("a_m","b_m") %in% names(manuf)) && !("flume" %in% names(manuf))) {
      df_fmt$`Manufacturer Eq. (L/s)` <- sprintf(paste0("%.",digits,"f"),
                                           manuf$a_m[1] * df_fmt$Depth_m^manuf$b_m[1])
    } else if (all(c("flume","a_m","b_m") %in% names(manuf))) {
      # if provided per flume, add a column per flume with "Manu." suffix
      for (j in seq_len(J)) {
        if (as.character(j) %in% manuf$flume) {
          a_m <- manuf$a_m[manuf$flume == j][1]
          b_m <- manuf$b_m[manuf$flume == j][1]
          colname <- paste0(flume_names[as.character(j)], " (Manu., L/s)")
          df_fmt[[colname]] <- sprintf(paste0("%.",digits,"f"),
                                      a_m * df_fmt$Depth_m^b_m)
        }
      }
    }
  }
  
  # also return a numeric version with separate median/lo/hi columns per flume (for CSV/analysis)
  df_num <- df_long %>%
    mutate(Flume = flume_names[flume]) %>%
    select(Depth_m, Flume, median_lps, lo_lps, hi_lps) %>%
    pivot_wider(names_from = Flume,
                values_from = c(median_lps, lo_lps, hi_lps),
                names_sep = "_") %>%
    arrange(Depth_m)
  if (!is.null(manuf)) {
    if (!("flume" %in% names(manuf))) {
      df_num$Manufacturer_lps <- manuf$a_m[1] *
        df_num$Depth_m^manuf$b_m[1]
    } else {
      for (j in seq_len(J)) {
        m <- subset(manuf, flume == j)
        if (nrow(m) > 0) {
          df_num[[paste0("Manufacturer_lps_", flume_names[as.character(j)])]] <-
            m$a_m[1] * df_num$Depth_m^m$b_m[1]
        }
      }
    }
  }
  
  list(display = df_fmt, numeric = df_num)
}

# generate_manuscript_table4
generate_manuscript_table4 <- function(post, prob = 0.95,
                                       flume_names = c("1"="3D Printed 1",
                                                       "2"="3D Printed 2",
                                                       "3"="Fiberglass")) {
  J <- ncol(post$a)
  S <- nrow(post$a)
  stopifnot(is.matrix(post$a), is.matrix(post$b),
            identical(dim(post$a), dim(post$b)),
            identical(dim(post$a), dim(post$sigma_proc)),
            identical(dim(post$a), dim(post$sigma_Q_add)),
            length(post$cv_Q) == S,
            all(as.character(seq_len(J)) %in% names(flume_names)),
            length(prob) == 1L, is.finite(prob), prob > 0, prob < 1)
  flume_names <- flume_names[as.character(seq_len(J))]
  parameters <- c("a", "b", "sigma_proc", "sigma_Q_add", "cv_Q")
  parameter_labels <- c("a", "b", "\u03c3proc", "\u03c3Q,add", "cv")
  descriptions <- c(
    "Scale coefficient in Q (L/s) = a * H (m)^b.",
    "Dimensionless exponent controlling rating-curve steepness.",
    "Process standard deviation (L/s).",
    "Additive discharge measurement-error standard deviation (L/s).",
    "Shared proportional flow-meter error (% of reading)."
  )
  units <- c("(L/s)/m^b", "dimensionless", "L/s", "L/s", "% of reading")
  digits <- c(2L, 2L, 3L, 3L, 3L)
  # cv_Q is a fraction in the model; multiply by 100 for percent of reading.
  metric_draws <- list(
    a = post$a,
    b = post$b,
    sigma_proc = post$sigma_proc,
    sigma_Q_add = post$sigma_Q_add,
    cv_Q = matrix(as.numeric(post$cv_Q) * 100, nrow = S, ncol = J)
  )
  display <- data.frame(Parameter = parameter_labels, Description = descriptions,
                        stringsAsFactors = FALSE, check.names = FALSE)
  rows <- vector("list", length(parameters) * J)
  for (j in seq_len(J)) {
    formatted <- character(length(parameters))
    for (k in seq_along(parameters)) {
      values <- metric_draws[[parameters[k]]][, j]
      stopifnot(all(is.finite(values)))
      interval <- as.numeric(rethinking::HPDI(values, prob = prob))
      center <- mean(values)
      fmt <- paste0("%.", digits[k], "f (%.", digits[k], "f \u2013 %.", digits[k], "f)")
      formatted[k] <- sprintf(fmt, center, interval[1], interval[2])
      rows[[(k - 1L) * J + j]] <- data.frame(
        parameter = parameters[k], flume = unname(flume_names[j]), units = units[k],
        posterior_mean = center, hpdi_lower = interval[1], hpdi_upper = interval[2],
        probability = prob, interval_type = "HPDI", posterior_draws = S,
        scope = if (parameters[k] == "cv_Q") "shared" else "flume-specific"
      )
    }
    display[[flume_names[j]]] <- formatted
  }
  list(display = display, numeric = do.call(rbind, rows))
}

# write_manuscript_csv
write_manuscript_csv <- function(table, filename) {
  connection <- file(filename, open = "wt", encoding = "UTF-8")
  on.exit(close(connection), add = TRUE)
  cat("\ufeff", file = connection)
  write.csv(table, connection, row.names = FALSE)
}
