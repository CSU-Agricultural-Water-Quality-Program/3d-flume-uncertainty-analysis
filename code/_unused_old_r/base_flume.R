# base_flume.R — main analysis entry point
# Bayesian Calibration of 3D-Printed Trapezoidal Flume
# AJ Brown
# Agricultural Water Quality Program
# Colorado State University


# Load necessary libraries
library(rethinking)
library(ggplot2)
library(dplyr)
library(tidyr)

# Native calibration and reporting units: stage in meters, discharge in L/s.
# Prior settings are fixed inputs, not estimated hyperparameters.
# Resolve from this file when sourced, or from getwd() when run line by line.
.source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(x) x$ofile))
.cli_file <- grep("^--file=", commandArgs(), value=TRUE)
.entry_file <- if (length(.source_files)) tail(unlist(.source_files),1L) else
  if (length(.cli_file)) sub("^--file=","",.cli_file[1]) else NULL
.project_starts <- c(if (!is.null(.entry_file))
  dirname(normalizePath(.entry_file,winslash="/",mustWork=TRUE)), getwd())
flume_repo_dir <- NULL
for (.start in unique(.project_starts)) {
  .candidate <- normalizePath(.start,winslash="/",mustWork=TRUE)
  repeat {
    if (file.exists(file.path(.candidate,"code/base_flume.R")) &&
        file.exists(file.path(.candidate,"code/supporting_code/project-paths.R"))) {
      flume_repo_dir <- .candidate
      break
    }
    .parent <- dirname(.candidate)
    if (identical(.parent,.candidate)) break
    .candidate <- .parent
  }
  if (!is.null(flume_repo_dir)) break
}
if (is.null(flume_repo_dir))
  stop("Cannot locate the flume project. Source code/base_flume.R as a file, ",
       "or set your working directory to the repository or one of its subfolders.")
source(file.path(flume_repo_dir,"code/supporting_code/project-paths.R"))
flume_repo_dir <- flume_project_root(flume_repo_dir)

# 1. Load fixed settings and supporting functions (no fitting here).
source(file.path(flume_repo_dir,"code/supporting_code/flume-metric-units.R"))
source(file.path(flume_repo_dir,"code/supporting_code/flume-reporting.R"))
source(file.path(flume_repo_dir,"code/supporting_code/flume-point-plots.R"))
source(file.path(flume_repo_dir,"code/supporting_code/flume-uncertainty-components.R"))
# Sampling requires an explicit opt-in; syntax/data checks do not sample.
RUN_CALIBRATION <- identical(Sys.getenv("FLUME_RUN_CALIBRATION"), "true")

# 2. Set input/output paths without changing the working directory.

# Current calibration products are separate from older Bayes and Python outputs.
# output_dir <- "../output/bayes/sep2026"
# output_prefix <- "sep2026"
# Keep the previous October fit intact while other uncertainty terms are reviewed.
output_dir <- file.path(flume_repo_dir,"output/bayes/oct2026_metric")
output_prefix <- "oct2026_metric"

output_file <- function(filename) file.path(output_dir, paste0(output_prefix, "_", filename))

# Load the flume calibration data
# Historical imperial inputs/fits are preserved separately. Prepare this CSV
# with code/tools/prepare-flume-metric-data.R; do not silently convert observations here.
input_csv <- file.path(flume_repo_dir,"data/New Data Oct 2026/flume_data_oct2026_metric.csv")
df_raw <- read.csv(input_csv, fileEncoding = "UTF-8-BOM")

# Drop flume calibration rows (i.e., where flume_type == EQ)
df <- df_raw %>% filter(flume_type != "EQ")

# create data list for Stan model
dat <- list(
  N       = nrow(df),
  Q_obs   = as.numeric(df$flow_lps),                 # observed discharge (L/s)
  H_obs   = as.numeric(df$depth_m),                  # observed head (m)
  F       = as.numeric(as.factor(df$flume_type)),    # flume index 1..J
  J       = length(unique(df$flume_type)),           # number of flume groups
  sigma_H = STAGE_ACCURACY_M,                       # stage-accuracy half-width, m; not an SD
  stage_units = "m", flow_units = "L/s"
)
assert_flume_metric_units(dat)
stopifnot(all(c("depth_m", "flow_lps", "flume_type") %in% names(df)),
          nrow(df) == 30L, all(table(df$flume_type) == 10L))
dat$a_prior_median <- MANUFACTURER_A
dat$b_prior_median <- MANUFACTURER_B
dat$rating_prior_log_sd <- RATING_PRIOR_LOG_SD
dat$sigma_proc_prior_rate <- SIGMA_PROC_PRIOR_RATE_LPS
dat$sigma_Q_add_prior_log_sd <- SIGMA_Q_ADD_PRIOR_LOG_SD
dat$cv_Q_prior_median <- CV_Q_PRIOR_MEDIAN
dat$cv_Q_prior_log_sd <- CV_Q_PRIOR_LOG_SD
dat$sigma_Q_add_prior_median <- SIGMA_Q_ADD_PRIOR_MEDIAN_LPS
dat$H_lower <- dat$H_obs - dat$sigma_H
dat$H_upper <- dat$H_obs + dat$sigma_H
stopifnot(all(is.finite(dat$H_obs)), all(dat$H_lower > 0),
          all(dat$H_upper > dat$H_lower))


# 3. Define and fit the Bayesian model; archived alternatives stay commented.

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
#       a ~ dlnorm(log(a_prior_median) , rating_prior_log_sd),
#       b ~ dlnorm(log(b_prior_median) , rating_prior_log_sd),
#       sigma_true ~ dexp(sigma_proc_prior_rate),
#       sigma_proc ~ dexp(sigma_proc_prior_rate)
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
#       vector[J]:a ~ dlnorm( log(a_prior_median) , rating_prior_log_sd ),
#       vector[J]:b ~ dlnorm( log(b_prior_median) , rating_prior_log_sd ),
#       
#       # Per-flume scatter
#       vector[J]:sigma_true ~ dexp(sigma_proc_prior_rate),
#       vector[J]:sigma_proc ~ dexp(sigma_proc_prior_rate)
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

# Conditional stage-error model: true stage is uniform within the reported
# +/- 1.5-mm bound around observed depth. This assumes the reported bound
# applies to computed depth (not separately to two gauge readings), no shared
# datum error, and a locally flat stage prior. Stage uncertainty now enters
# the discharge likelihood through H_true. Other priors remain under review.

if (!RUN_CALIBRATION) {
  stop("Metric inputs and model are prepared. Sampling is disabled; explicitly set FLUME_RUN_CALIBRATION=true to fit.")
}
dir.create(output_dir, showWarnings=FALSE, recursive=TRUE)
fit_seed <- 1815795276
flume_model_3 <- ulam(
  alist(
    vector[N]:H_true ~ dunif(H_lower, H_upper),

    # Latent true discharge with per-flume process scatter (vectorized)
    vector[N]:Q_true ~ dnorm( mu , sigma_proc[F] ),
    
    # Vectorized mean by flume
    mu <- a[F] * H_true[i]^b[F],
    
    # Observation model: heteroscedastic SD
    # empirical proportional SD + per-flume additive uncertainty attributable
    # to infiltration and potentially other unresolved measurement effects
    Q_obs ~ dnorm(
      Q_true,
      sqrt(
        square( cv_Q * Q_true ) +          # global % of reading
          square( sigma_Q_add[F] )           # per-flume additive uncertainty
      )
    ),
    
    # Per-flume rating-curve parameters
    vector[J]:a ~ dlnorm( log(a_prior_median) , rating_prior_log_sd ),
    vector[J]:b ~ dlnorm( log(b_prior_median) , rating_prior_log_sd ),
    
    # Per-flume process scatter
    vector[J]:sigma_proc ~ dexp(sigma_proc_prior_rate),
    
    # Error hyperparameters
    # Previous September 2026 proportional-error prior, retained for reference:
    # cv_Q ~ dlnorm( log(2e-4) , 0.5 ),        # median 0.02% of reading
    # October revision: median 0.2% = 0.002; prior width remains 0.5.
    # This is an empirical SD prior informed by the proportional specification,
    # not an assertion that the manufacturer's accuracy limit equals 1 SD.
    cv_Q ~ dlnorm( log(cv_Q_prior_median) , cv_Q_prior_log_sd ),
    # Retain the lab allowance as the existing additive-SD prior median.
    # It does not encode a strict bound or estimate a directional loss.
    vector[J]:sigma_Q_add ~ dlnorm( log(sigma_Q_add_prior_median) , sigma_Q_add_prior_log_sd )
  ),
  data = dat[setdiff(names(dat), c("stage_units", "flow_units"))],
  constraints = list(H_true = "lower=H_lower,upper=H_upper"),
  chains = 4, cores = 4, iter = 4000, seed = fit_seed
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
#     vector[N]:H_true ~ dunif(H_lower, H_upper),
#
#     # Latent true discharge with per-flume process scatter (vectorized)
#     vector[N]:Q_true ~ dnorm( mu , sigma_proc[F] ),
#
#     # Mean computed in log space, then exponentiated back to original scale
#     # (this is mathematically identical to mu <- a[F] * H_obs^b[F] when a[F] = exp(alpha[F]))
#     mu <- exp( alpha[F] + b[F] * log(H_true[i]) ),
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
#     # With native metric a, alpha = log(a) has the same fixed prior center and width.
#     vector[J]:alpha ~ dnorm( log(a_prior_median) , rating_prior_log_sd ),
#
#     # Same b prior as Model 3
#     vector[J]:b ~ dlnorm( log(b_prior_median) , rating_prior_log_sd ),
#
#     # Same per-flume process scatter
#     vector[J]:sigma_proc ~ dexp(sigma_proc_prior_rate),
#
#     # Same error hyperparameters
#     cv_Q ~ dlnorm( log(cv_Q_prior_median) , cv_Q_prior_log_sd ),
#     vector[J]:sigma_Q_add ~ dlnorm( log(sigma_Q_add_prior_median) , sigma_Q_add_prior_log_sd )
#   ),
#   data = dat[setdiff(names(dat), c("stage_units", "flow_units"))],
#   constraints = list(H_true = "lower=H_lower,upper=H_upper"),
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
# 4. Export rating-curve, point, and uncertainty-component figures.
# Published metric equation: Q (L/s) = a_m * H (m)^b_m.
manuf <- data.frame(
  a_m   = c(MANUFACTURER_A),
  b_m   = c(MANUFACTURER_B)
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
# Horizontal bars display the reported stage bound; the fit uses latent stage.
point_products <- save_flume_point_plots(
  post3, dat, manuf, output_dir, prefix = output_prefix, input_csv = input_csv,
  fit_file = output_file("flume_model_3_cmdstanfit.rds")
)
print(point_products$plot)
print(point_products$departures)

# Permanent predictive-uncertainty components and variance shares.
# The helper only defines functions; use the same joint draws as all other plots.
uncertainty_products <- save_flume_uncertainty_components(
  post3, dat, manuf, output_dir, prefix = output_prefix, input_csv = input_csv,
  fit_file = output_file("flume_model_3_cmdstanfit.rds")
)
print(uncertainty_products$sd)
print(uncertainty_products$shares)


# 5. Export stage-discharge predictions and manuscript Table 4.
# Report stage in metres at 0.01-m increments, retaining the physical extent.
depth_grid_m <- seq.int(0L, 24L) / 100
manuf <- data.frame(a_m=MANUFACTURER_A, b_m=MANUFACTURER_B)

res <- generate_flume_table(
  post = post3,
  depths_m = depth_grid_m,
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


# Parameter summaries use all saved draws.
table4 <- generate_manuscript_table4(post3, prob = 0.95)
write_manuscript_csv(table4$display, output_file("manuscript_table4.csv"))
write.csv(table4$numeric, output_file("manuscript_table4_numeric.csv"),
          row.names = FALSE)

