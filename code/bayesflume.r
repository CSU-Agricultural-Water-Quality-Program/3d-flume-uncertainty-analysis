# bayesflume.r
# Bayesian Calibration of 3D-Printed Trapezoidal Flume
# AJ Brown
# Agricultural Water Quality Program
# Colorado State University


# Load necessary libraries
library(rethinking)
library(ggplot2)
library(dplyr)
library(tidyr)

# Calibration inputs and priors remain in ft and cfs. Convert only reporting
# quantities, using the international foot (exactly 0.3048 m).
# Source: NIST SP 811, Appendix B.8:
# https://www.nist.gov/pml/special-publication-811/nist-guide-si-appendix-b-conversion-factors/nist-guide-si-appendix-b8
FT_TO_CM <- 30.48
CFS_TO_LPS <- (FT_TO_CM / 100)^3 * 1000  # exactly 28.316846592 L/s per cfs

# Transform complete chains before summarizing: a_cm = a * CFS_TO_LPS /
# FT_TO_CM^b, so Q in L/s = a_cm * (H in cm)^b. The sampled b must be used
# in each draw; transforming marginal means or interval endpoints is invalid.
metric_precis <- function(fit, prob = 0.89) {
  original <- precis(fit, depth = 2, prob = prob)
  if (!is.null(attr(fit, "cstanfit"))) {
    draws <- attr(fit, "cstanfit")$draws(
      variables = rownames(original), format = "draws_array"
    )
  } else {
    draws <- posterior::as_draws_array(as.array(attr(fit, "stanfit")))
    draws <- draws[, , rownames(original), drop = FALSE]
  }
  variables <- dimnames(draws)[[3]]
  metric_names <- variables
  for (k in seq_along(variables)) {
    v <- variables[k]
    if (grepl("^(a|alpha)\\[", v)) {
      b_name <- sub("^(a|alpha)", "b", v)
      if (grepl("^a\\[", v)) {
        draws[, , k] <- draws[, , k] * CFS_TO_LPS / FT_TO_CM^draws[, , b_name]
        metric_names[k] <- sub("^a", "a_cm", v)
      } else {
        draws[, , k] <- draws[, , k] + log(CFS_TO_LPS) -
          draws[, , b_name] * log(FT_TO_CM)
        metric_names[k] <- sub("^alpha", "alpha_cm", v)
      }
    } else if (grepl("^(Q_true|mu|sigma_proc|sigma_Q_add)(\\[|$)", v)) {
      draws[, , k] <- draws[, , k] * CFS_TO_LPS
      metric_names[k] <- sub("^(Q_true|mu|sigma_proc|sigma_Q_add)", "\\1_lps", v)
    }
  }
  dimnames(draws)[[3]] <- metric_names
  result <- as.data.frame(posterior::summarise_draws(
    draws, mean = base::mean, sd = stats::sd,
    interval = function(x) stats::quantile(x, probs = c((1 - prob)/2, (1 + prob)/2)),
    rhat = posterior::rhat, ess_bulk = posterior::ess_bulk
  ))
  rownames(result) <- result$variable
  result$variable <- NULL
  methods::new("precis", result, digits = 3)
}

# Keep units explicit in CSVs; b and cv_Q are dimensionless. alpha_cm is the
# log of the numerical coefficient a_cm when stage is cm and flow is L/s.
write_metric_precis <- function(summary, filename) {
  result <- as.data.frame(summary)
  v <- rownames(result)
  result$units <- ifelse(grepl("_lps", v), "L/s",
                        ifelse(grepl("^a_cm\\[", v), "(L/s)/cm^b",
                               ifelse(grepl("^alpha_cm\\[", v), "log(a_cm)", "dimensionless")))
  write.csv(result, filename, row.names = TRUE)
}

# rethinking 2.42's dashboard() assumes divergent__ and energy__ occupy columns
# 2 and 3. CmdStanR's restored CSV fits return them in columns 5 and 6. Supply a reporting-only
# backend proxy with named diagnostics reordered; keep the fitted object intact.
# Source: https://github.com/rmcelreath/rethinking/blob/master/R/dashboard.r
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

# dashboard() draws with base graphics. Always close its file device, even
# if plotting fails; the on-screen dashboard is still drawn separately.
save_dashboard <- function(fit, filename) {
  png(filename, width = 8, height = 8, units = "in", res = 300)
  on.exit(dev.off(), add = TRUE)
  dashboard_flume(fit)
  invisible(filename)
}

# Export post-warmup diagnostics by chain and by iteration. Energy is an HMC
# diagnostic with an arbitrary reference level; E-BFMI measures energy movement
# relative to the chain's marginal energy variance, not whether its mean is low.
# https://mc-stan.org/cmdstanr/reference/fit-method-diagnostic_summary.html
export_stan_diagnostics <- function(fit, output_dir, prefix) {
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
              "run_started_utc", "rethinking_version", "cmdstanr_version", "cmdstan_version"),
    value = c(normalizePath(input_csv, winslash = "/"), unname(tools::md5sum(input_csv)),
              "ft", "cfs", "cm", "L/s", metadata$seed, dim(x)[2],
              metadata$iter_warmup, metadata$iter_sampling, metadata$adapt_delta,
              metadata$max_treedepth, metadata$start_datetime[1],
              as.character(packageVersion("rethinking")), as.character(packageVersion("cmdstanr")),
              as.character(cmdstanr::cmdstan_version()))
  )
  write.csv(run_metadata, file.path(output_dir, paste0(prefix, "_run_metadata.csv")),
            row.names = FALSE)
  invisible(chain_summary)
}

# set working directory
setwd("C:/Users/ansle/OneDrive/Documents/GitHub/3d-flume-uncertainty-analysis/data")

# Current calibration products are separate from older Bayes and Python outputs.
# output_dir <- "../output/bayes/sep2026"
# output_prefix <- "sep2026"
output_dir <- "../output/bayes/oct2026"
output_prefix <- "oct2026"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
output_file <- function(filename) file.path(output_dir, paste0(output_prefix, "_", filename))

# Load the flume calibration data
# df_raw <- read.csv("../data/flume_data.csv")
# df_raw <- read.csv("../data/sim_data.csv") # for testing with simulated data
# input_csv <- "../data/New Data Sep 2026/flume_data_sep2026.csv"
input_csv <- "../data/New Data Oct 2026/flume_data_oct2026.csv"
df_raw <- read.csv(input_csv, fileEncoding = "UTF-8-BOM")

# Drop flume calibration rows (i.e., where flume_type == EQ)
df <- df_raw %>% filter(flume_type != "EQ")

# create data list for Stan model
dat <- list(
  N       = nrow(df),
  Q_obs   = as.numeric(df$flow_cfs),                 # observed discharge (cfs)
  H_obs   = as.numeric(df$depth_ft),                 # observed head (ft)
  F       = as.numeric(as.factor(df$flume_type)),    # flume index 1..J
  J       = length(unique(df$flume_type)),           # number of flume groups
  # sigma_H = 0.005                                  # previous gauge accuracy, ft
  sigma_H = 0.01 / FT_TO_CM                          # October accuracy: +/- 0.01 cm, reported separately
)


# Define the Bayesian models

# Model 1.0 
# - basic model, not looking at flume types separately
# - True flow Q is latent variable
# flume_model_1 <- ulam(
#     alist(
#       # Observation model for discharge
#       Q_obs ~ dnorm(Q_true, sigma_true),
#       
#       # Process model for true discharge (power-law)
#       vector[N]:Q_true ~ dnorm(a * H_obs^b, sigma_proc),
#       
#       # Priors
#       a ~ dlnorm(log(1.55) , 0.5),
#       b ~ dlnorm(log(2.58) , 0.5),
#       sigma_true ~ dexp(1),
#       sigma_proc ~ dexp(1)
#     ),
#     data = dat, chains = 4, cores = 4, iter = 2000, log_lik = TRUE
#   )

# Summarize the model
# precis(flume_model_1, depth = 2)
# plot(precis(flume_model_1, depth = 2))


# Model 2.0
# - hierarchical model, looking at flume types separately
# - True flow Q is latent variable
# flume_model_2 <- ulam(
#     alist(
#       # Observation model for discharge (sigma varies by flume via F)
#       Q_obs ~ dnorm( Q_true , sigma_true[F] ),
#       
#       # Process model for latent true discharge (sigma varies by flume)
#       vector[N]:Q_true ~ dnorm( mu , sigma_proc[F] ),
#       
#       # Vectorized mean, per row using flume-indexed a and b
#       mu <- a[F] * H_obs^b[F],
#       
#       # Per-flume priors
#       vector[J]:a ~ dlnorm( log(1.55) , 0.5 ),
#       vector[J]:b ~ dlnorm( log(2.58) , 0.5 ),
#       
#       # Per-flume scatter
#       vector[J]:sigma_true ~ dexp(1),
#       vector[J]:sigma_proc ~ dexp(1)
#     ),
#     data = dat,
#     chains = 4, cores = 4, iter = 2000, log_lik = TRUE
#   )


# Summarize the model
# precis(flume_model_2, depth = 2)
# plot(precis(flume_model_2, depth = 2))


# Model 3.0
# - hierarchical model, looking at flume types separately
# - True flow Q is observed with error
# - includes measurement error explicitly

# H_obs is treated as fixed. sigma_H is recorded in dat but is not currently
# used in this likelihood; adding latent stage is a separate model change.

fit_seed <- 1815795276
flume_model_3 <- ulam(
  alist(
    # Latent true discharge with per-flume process scatter (vectorized)
    vector[N]:Q_true ~ dnorm( mu , sigma_proc[F] ),
    
    # Vectorized mean by flume
    mu <- a[F] * H_obs^b[F],
    
    # Observation model: heteroscedastic SD
    # global percent-of-reading + per-flume additive infiltration
    Q_obs ~ dnorm(
      Q_true,
      sqrt(
        square( cv_Q * Q_true ) +          # global % of reading
          square( sigma_Q_add[F] )           # per-flume infiltration
      )
    ),
    
    # Per-flume rating-curve parameters
    vector[J]:a ~ dlnorm( log(1.55) , 0.5 ),
    vector[J]:b ~ dlnorm( log(2.58) , 0.5 ),
    
    # Per-flume process scatter
    vector[J]:sigma_proc ~ dexp(1),
    
    # Error hyperparameters
    # Previous September 2026 proportional-error prior, retained for reference:
    # cv_Q ~ dlnorm( log(2e-4) , 0.5 ),        # median 0.02% of reading
    # October 2026: center the prior at the reported 2% of reading.
    # cv_Q is a dimensionless fraction: 2% = 0.02, not 2 or 0.0002.
    # The lognormal median is 0.02; retain the previous log-scale SD of 0.5.
    cv_Q ~ dlnorm( log(0.02) , 0.5 ),          # shared across all flumes
    # This additive prior is in cfs, not percent; retain its existing value.
    vector[J]:sigma_Q_add ~ dlnorm( log(0.02) , 0.5 )  # per flume
  ),
  data = dat, chains = 4, cores = 4, iter = 4000, seed = fit_seed
)

# Preserve complete draws and diagnostics for later exports without refitting.
attr(flume_model_3, "cstanfit")$save_object(
  output_file("flume_model_3_cmdstanfit.rds")
)
chain_diagnostics <- export_stan_diagnostics(flume_model_3, output_dir,
                                            paste0(output_prefix, "_flume_model_3"))

# model convergence check
traceplot_ulam(flume_model_3, ask=F, window=c(2100,4000))
trankplot(flume_model_3, ask=F, window=c(2100,4000))
dashboard_flume(flume_model_3)
save_dashboard(flume_model_3, output_file("flume_model_3_dashboard.png"))

# Summarize the model and save output
p3 <- metric_precis(flume_model_3)
write_metric_precis(p3, output_file("flume_model_3_precis.csv"))
p3

# plot model summary and save as jpg
jpeg(output_file("flume_model_3_precis.jpg"), width = 6.5, height = 6.5, units = "in", res = 600)
plot(p3)
dev.off()
par(mfrow=c(1,1))
plot(p3)


post3 <- extract.samples(flume_model_3)
# Reproduce the posterior-draw selections used in the reporting products.
set.seed(fit_seed)


# Model 4 is disabled and retained for reference. It has consistently shown
# poorer convergence in the project's comparisons (more divergences and poorer
# chain diagnostics). This is an empirical observation, not a general claim
# about log parameterization: its flow likelihood remains on the original scale.
# # Model 4: same as Model 3, but reparameterize a = exp(alpha)
# # so the mean is computed via log space:
# #   mu = exp(alpha[F] + b[F] * log(H))
# #
# # Likelihood and uncertainty structure remain on original Q scale:
# #   Q_true ~ Normal(mu, sigma_proc[F])
# #   Q_obs  ~ Normal(Q_true, sqrt((cv_Q*Q_true)^2 + sigma_Q_add[F]^2))
#
# # Model 4.0
# # - identical likelihood / uncertainty structure to Model 3
# # - reparameterizes a[F] as exp(alpha[F]) so the power law mean is computed in log space
# # - goal: improve sampling geometry for a/b without changing your manuscript’s uncertainty story
#
# flume_model_4 <- ulam(
#   alist(
#     # Latent true discharge with per-flume process scatter (vectorized)
#     vector[N]:Q_true ~ dnorm( mu , sigma_proc[F] ),
#
#     # Mean computed in log space, then exponentiated back to original scale
#     # (this is mathematically identical to mu <- a[F] * H_obs^b[F] when a[F] = exp(alpha[F]))
#     mu <- exp( alpha[F] + b[F] * log(H_obs + 1e-9) ),
#
#     # Observation model: SAME as Model 3 (heteroscedastic SD)
#     Q_obs ~ dnorm(
#       Q_true,
#       sqrt(
#         square( cv_Q * Q_true ) +
#           square( sigma_Q_add[F] )
#       )
#     ),
#
#     # Priors: alpha is the log of a, consistent with your original a prior
#     # If a ~ LogNormal(log(1.55), 0.5), then alpha = log(a) ~ Normal(log(1.55), 0.5)
#     vector[J]:alpha ~ dnorm( log(1.55) , 0.5 ),
#
#     # Same b prior as Model 3
#     vector[J]:b ~ dlnorm( log(2.58) , 0.5 ),
#
#     # Same per-flume process scatter
#     vector[J]:sigma_proc ~ dexp(1),
#
#     # Same error hyperparameters
#     cv_Q ~ dlnorm( log(2e-4) , 0.5 ),
#     vector[J]:sigma_Q_add ~ dlnorm( log(0.02) , 0.5 )
#   ),
#   data = dat,
#   chains = 4, cores = 4, iter = 4000
#   #log_lik = TRUE
# )
#
#
# # model convergence check
# traceplot_ulam(flume_model_4, ask=F, window=c(2100,4000))
# trankplot(flume_model_4, ask=F, window=c(2100,4000))
# dashboard_flume(flume_model_4)
# save_dashboard(flume_model_4, output_file("flume_model_4_dashboard.png"))
#
# # Summarize the model and save output
# p4 <- metric_precis(flume_model_4)
# write_metric_precis(p4, output_file("flume_model_4_precis.csv"))
# p4
#
# # plot model summary and save as jpg
# jpeg(output_file("flume_model_4_precis.jpg"), width = 6.5, height = 6.5, units = "in", res = 600)
# plot(p4)
# dev.off()
# par(mfrow=c(1,1))
# plot(p4)
#
#
#
######################### Graphing Outputs #####################################

# Stage discharge curves by flume type with uncertainty intervals
#' Plot Bayesian stage–discharge curves with posterior ribbons
#'
#' @param post  list returned by rethinking::extract.samples(ulam_fit)
#'              expects components: a[J], b[J], sigma_proc[J], cv_Q, sigma_Q_add[J]
#' @param dat   list you passed to ulam, needs: Q_obs[N], H_obs[N], F[N], J
#' @param prob  probability mass for ribbons (default 0.95)
#' @param ndraws number of posterior draws to use for ribbons (default 4000)
#' @param extrapolation_fraction extend the grid above max observed stage (default 0.10)
#' @param show_manufacturer optional data.frame with columns flume (factor or int 1..J),
#'              a_m (manufacturer a), b_m (manufacturer b) to add dotted curves
#' @return      ggplot object with stage in cm and flow in L/s
# --- Dependencies ---
# install.packages(c("ggplot2","dplyr","tidyr"))

plot_flume_posterior <- function(post, dat, prob = 0.95, ndraws = 4000,
                                 show_manufacturer = NULL,
                                 ribbon_kind = c("process", "parameter", "predictive"),
                                 show_qobs_errors = FALSE,
                                 extrapolation_fraction = 0.10,
                                 flume_names = c("1" = "3D Printed 1",
                                                 "2" = "3D Printed 2",
                                                 "3" = "Fiberglass")) {
  
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
  
  # Convert every plotted quantity, including intervals and manufacturer flow.
  to_metric <- function(x) {
    x$H <- x$H * FT_TO_CM
    q_columns <- grep("^Q", names(x), value = TRUE)
    x[q_columns] <- lapply(x[q_columns], function(q) q * CFS_TO_LPS)
    x
  }
  ribbon_df <- to_metric(ribbon_df)
  Qtrue_df <- to_metric(Qtrue_df)
  Qobs_df <- to_metric(Qobs_df)
  if (!is.null(manuf_df)) manuf_df <- to_metric(manuf_df)
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
    coord_cartesian(xlim = c(0, max(Hgrid) * FT_TO_CM),
                    ylim = c(0, y_max), expand = FALSE) +
    labs(x = "Stage (cm)", y = "Flow rate (L/s)") +
    theme_bw(base_size = 13) +
    theme(legend.position = "bottom", legend.box = "vertical",
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 11),
          legend.key.width = grid::unit(1.4, "cm"),
          panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
          panel.grid.minor = element_blank())
}






# Assuming you already have:
#   post3 <- extract.samples(flume_model_3)
#   dat   <- list(N=nrow(df), Q_obs=df$flow_cfs, H_obs=df$depth_ft,
#                 F=as.integer(as.factor(df$flume_type)), J=length(unique(df$flume_type)))

# Optional manufacturer coefficients by flume index (edit as needed)
# Equation: Q (cfs) = a_m * H (ft)^b_m; converted only for reporting.
manuf <- data.frame(
  a_m   = c(1.55),
  b_m   = c(2.58)
)

# Generate the plot
g <- plot_flume_posterior(
  post  = post3,              # posterior samples from ulam
  dat   = dat,                # data list used for the model
  prob  = 0.95,               # credible interval probability
  ndraws = nrow(post3$a),     # all draws, shared with point-only summaries
  show_manufacturer = manuf,  # add manufacturer equation line
  ribbon_kind = "parameter",  # options: "parameter", "process", or "predictive"
  show_qobs_errors = FALSE,   # measurement error bars for Q_obs off (recommended)
  flume_names = c("1" = "3D Printed 1",
                  "2" = "3D Printed 2",
                  "3" = "Fiberglass")
)

# Print to screen
print(g)

# Save a high-resolution version
ggsave(output_file("flume_posterior_curves_final.jpg"), g,
       width = 9, height = 6, dpi = 600)

# Separate panels reveal each uncertainty band without changing the intervals,
# measured coordinates, or the shared zero-to-extrapolation axis limits.
g_panels <- g + facet_wrap(
  vars(flume), nrow = 1,
  labeller = labeller(flume = c("1" = "3D Printed 1", "2" = "3D Printed 2",
                               "3" = "Fiberglass"))
)
ggsave(output_file("flume_posterior_curves_by_flume.png"), g_panels,
       width = 11, height = 5, dpi = 300)

# Additive point-only comparison and departures from the manufacturer curve.
# Gauge accuracy is displayed separately; this does not change the likelihood.
source("../code/flume-point-plots.r")
point_products <- save_flume_point_plots(
  post3, dat, manuf, output_dir, prefix = output_prefix, input_csv = input_csv,
  fit_file = output_file("flume_model_3_cmdstanfit.rds")
)
print(point_products$plot)
print(point_products$departures)

# Permanent predictive-uncertainty components and variance shares.
# The helper only defines functions; use the same joint draws as all other plots.
source("../code/flume-uncertainty-components.r")
uncertainty_products <- save_flume_uncertainty_components(
  post3, dat, manuf, output_dir, prefix = output_prefix, input_csv = input_csv,
  fit_file = output_file("flume_model_3_cmdstanfit.rds")
)
print(uncertainty_products$sd)
print(uncertainty_products$shares)


# Generative predictions output
# Generate posterior summaries of Q(H) on a grid for each flume
# post  : extract.samples(ulam_fit)
# depths_cm: numeric vector of stage (cm) at which to report results
# prob  : central interval mass (e.g., 0.95)
# kind  : "parameter" (default), "process", or "predictive"
# manuf : data.frame with either global a_m,b_m columns, or per-flume (flume,a_m,b_m)
# flume_names: named character vector mapping "1","2","3",... -> labels to use
generate_flume_table <- function(post, depths_cm = seq.int(0L, 24L),
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
  stopifnot(length(depths_cm) > 0, all(is.finite(depths_cm)), all(depths_cm >= 0),
            all(is.finite(a)), all(a > 0), all(is.finite(b)), all(b > 0))
  H <- depths_cm / FT_TO_CM  # evaluate the unchanged ft/cfs calibration
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
    
    data.frame(Depth_cm = depths_cm,
               flume = names(flume_names)[j],
               median_lps = md * CFS_TO_LPS,
               lo_lps = lo * CFS_TO_LPS, hi_lps = hi * CFS_TO_LPS)
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
    select(Depth_cm, Flume, label) %>%
    pivot_wider(names_from = Flume, values_from = label) %>%
    arrange(Depth_cm)
  
  # manufacturer column
  if (!is.null(manuf)) {
    if (all(c("a_m","b_m") %in% names(manuf)) && !("flume" %in% names(manuf))) {
      df_fmt$`Manufacturer Eq. (L/s)` <- sprintf(paste0("%.",digits,"f"),
                                           manuf$a_m[1] * (df_fmt$Depth_cm / FT_TO_CM)^manuf$b_m[1] * CFS_TO_LPS)
    } else if (all(c("flume","a_m","b_m") %in% names(manuf))) {
      # if provided per flume, add a column per flume with "Manu." suffix
      for (j in seq_len(J)) {
        if (as.character(j) %in% manuf$flume) {
          a_m <- manuf$a_m[manuf$flume == j][1]
          b_m <- manuf$b_m[manuf$flume == j][1]
          colname <- paste0(flume_names[as.character(j)], " (Manu., L/s)")
          df_fmt[[colname]] <- sprintf(paste0("%.",digits,"f"),
                                      a_m * (df_fmt$Depth_cm / FT_TO_CM)^b_m * CFS_TO_LPS)
        }
      }
    }
  }
  
  # also return a numeric version with separate median/lo/hi columns per flume (for CSV/analysis)
  df_num <- df_long %>%
    mutate(Flume = flume_names[flume]) %>%
    select(Depth_cm, Flume, median_lps, lo_lps, hi_lps) %>%
    pivot_wider(names_from = Flume,
                values_from = c(median_lps, lo_lps, hi_lps),
                names_sep = "_") %>%
    arrange(Depth_cm)
  if (!is.null(manuf)) {
    if (!("flume" %in% names(manuf))) {
      df_num$Manufacturer_lps <- manuf$a_m[1] *
        (df_num$Depth_cm / FT_TO_CM)^manuf$b_m[1] * CFS_TO_LPS
    } else {
      for (j in seq_len(J)) {
        m <- subset(manuf, flume == j)
        if (nrow(m) > 0) {
          df_num[[paste0("Manufacturer_lps_", flume_names[as.character(j)])]] <-
            m$a_m[1] * (df_num$Depth_cm / FT_TO_CM)^m$b_m[1] * CFS_TO_LPS
        }
      }
    }
  }
  
  list(display = df_fmt, numeric = df_num)
}

# grid and manufacturer equation
# Evaluate predictions at integer centimetres, rather than rounding labels for
# predictions evaluated at converted feet. Retain approximately the former
# 0-0.8 ft reporting extent (0-24.4 cm), now at 0, 1, ..., 24 cm.
depth_grid_cm <- seq.int(0L, 24L)
manuf <- data.frame(a_m = 1.55, b_m = 2.58)  # global manufacturer curve

res <- generate_flume_table(
  post = post3,
  depths_cm = depth_grid_cm,
  prob = 0.95,
  kind = "parameter",     # "parameter" for comparing flume calibration curves
  ndraws = 4000,
  manuf = manuf,
  flume_names = c("1"="3D Printed 1","2"="3D Printed 2","3"="Fiberglass"),
  digits = 3
)

# Quick view for reporting (one column per flume with "median [lo, hi]")
res$display

# Numeric table for export or further analysis
res$numeric

write.csv(res$display, output_file("flume_generated_display.csv"), row.names = FALSE)
write.csv(res$numeric, output_file("flume_generated_numeric.csv"), row.names = FALSE)


# Manuscript Table 4: posterior means with 95% HPDIs, using every saved draw.
# This differs from the prediction table's medians and equal-tailed intervals,
# and from precis()'s default 89% interval. Transform paired draws before
# computing means or intervals; do not convert marginal summary endpoints.
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
    "Scale coefficient in Q (L/s) = a * H (cm)^b.",
    "Dimensionless exponent controlling rating-curve steepness.",
    "Process standard deviation (L/s).",
    "Additive discharge measurement-error standard deviation (L/s).",
    "Shared proportional flow-meter error (% of reading)."
  )
  units <- c("(L/s)/cm^b", "dimensionless", "L/s", "L/s", "% of reading")
  digits <- c(5L, 2L, 3L, 3L, 3L)
  # cv_Q is a fraction in the model; multiply by 100 for percent of reading.
  metric_draws <- list(
    a = post$a * CFS_TO_LPS / FT_TO_CM^post$b,
    b = post$b,
    sigma_proc = post$sigma_proc * CFS_TO_LPS,
    sigma_Q_add = post$sigma_Q_add * CFS_TO_LPS,
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

# UTF-8 BOM keeps Greek symbols and interval dashes readable when the CSV is
# opened in Excel for copying into the manuscript.
write_manuscript_csv <- function(table, filename) {
  connection <- file(filename, open = "wt", encoding = "UTF-8")
  on.exit(close(connection), add = TRUE)
  cat("\ufeff", file = connection)
  write.csv(table, connection, row.names = FALSE)
}

table4 <- generate_manuscript_table4(post3, prob = 0.95)
write_manuscript_csv(table4$display, output_file("manuscript_table4.csv"))
write.csv(table4$numeric, output_file("manuscript_table4_numeric.csv"),
          row.names = FALSE)

