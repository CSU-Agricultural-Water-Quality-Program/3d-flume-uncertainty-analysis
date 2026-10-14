# bayes_flume.R — main calibration and reporting script
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
# Locate the project from this file, or from the working directory for Run.
.source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(x) x$ofile))
.cli_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE))
.entry_file <- if (length(.source_files)) tail(unlist(.source_files), 1L) else
  if (length(.cli_file)) .cli_file[1] else NULL
flume_repo_dir <- if (is.null(.entry_file)) getwd() else
  dirname(normalizePath(.entry_file, winslash="/"))
while (!file.exists(file.path(flume_repo_dir, "code/supporting_code/project-paths.R"))) {
  .parent <- dirname(flume_repo_dir)
  if (identical(.parent, flume_repo_dir))
    stop("Set the working directory to the flume repository or source code/bayes_flume.R.")
  flume_repo_dir <- .parent
}
source(file.path(flume_repo_dir, "code/supporting_code/project-paths.R"))
flume_repo_dir <- flume_project_root(flume_repo_dir)

# 1. Load supporting functions; all calibration settings are defined below.
source(file.path(flume_repo_dir,"code/supporting_code/flume-metric-units.R"))
source(file.path(flume_repo_dir,"code/supporting_code/flume-reporting.R"))
source(file.path(flume_repo_dir,"code/supporting_code/flume-point-plots.R"))
source(file.path(flume_repo_dir,"code/supporting_code/flume-uncertainty-components.R"))

# 2. Set input/output paths without changing the working directory.

# Save the metric calibration separately from historical fits.
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

# Data and editable uncertainty/prior settings for Stan.
# To change a setting, edit and rerun this data/settings section (including
# H_lower/H_upper below for a stage-bound change), then rerun the model.
# Prior medians and widths below are fixed inputs, not fitted hyperparameters.
dat <- list(
  N       = nrow(df),
  Q_obs   = as.numeric(df$flow_lps),                 # observed discharge (L/s)
  H_obs   = as.numeric(df$depth_m),                  # observed head (m)
  F       = as.numeric(as.factor(df$flume_type)),    # flume index 1..J
  J       = length(unique(df$flume_type)),           # number of flume groups
  stage_units = "m", flow_units = "L/s",

  # POINT GAUGE: stage accuracy supplied by the study team, +/-1.5 mm.
  # Source: reported calibration allowance; a gauge manual was not verified.
  # Model assumption: uniform true depth within this bound, independently per
  # observation. Applied to computed depth, not separately to two readings.
  sigma_H = 0.0015,                                # half-width in m, not an SD

  # RATING CURVE: Openchannelflow, Large 60-Degree-V Trapezoidal Flume
  # Discharge Table, pp. 1-2: Q (L/s) = 941.2 * H (m)^2.58.
  # Published coefficients set lognormal prior medians; width is our choice.
  a_prior_median = 941.2,                          # scale: (L/s)/m^b
  b_prior_median = 2.58,                           # dimensionless exponent
  rating_prior_log_sd = 0.5,                       # dimensionless; chosen prior width

  # FLOW METER: Badger Meter ModMAG M2000 User Manual, MAG-UM-01272-EN-43
  # (August 2026), pp. 66 and 75: +/-0.2% of reading plus +/-1 mm/s,
  # under the stated reference/installation conditions.
  # Model choice: the current empirical proportional-SD prior median is 1.5%.
  # The manual does not identify this accuracy limit as 1 SD. Its velocity
  # term is not modeled explicitly in this empirical error specification.
  # cv_Q_prior_median = 0.002,                      # SD/flow fraction; 0.2% = 0.002
  cv_Q_prior_median = 0.015,                        # SD/flow fraction; 1.5% = 0.015
  cv_Q_prior_log_sd = 0.5,                          # dimensionless; chosen prior width

  # FLUME INFILTRATION: CSU Hydraulics Laboratory allowance reported by the
  # study team, +/-0.02 cfs; this is a lab allowance, not a meter specification.
  # 0.02 cfs * 28.316846592 L/s per cfs = 0.56633693184 L/s.
  # Model choice: this magnitude centers the additive-SD prior. It is neither
  # a hard error bound nor a fitted mean loss; unresolved additive error may
  # also contribute. The prior width is our choice, not supplied by the lab.
  sigma_Q_add_prior_median = 0.56633693184,        # additive SD in L/s
  sigma_Q_add_prior_log_sd = 0.5,                  # dimensionless; chosen prior width

  # PROCESS SCATTER: variation in true discharge around the rating curve.
  # Model choice, not an instrument accuracy. Preserve the previous Exp(1)
  # prior in cfs: mean 1 cfs = 28.316846592 L/s, rate = 1/28.316846592.
  sigma_proc_prior_rate = 0.0353146667214886       # rate in (L/s)^-1
)
assert_flume_metric_units(dat)
stopifnot(all(c("depth_m", "flow_lps", "flume_type") %in% names(df)),
          nrow(df) == 30L, all(table(df$flume_type) == 10L))
dat$H_lower <- dat$H_obs - dat$sigma_H
dat$H_upper <- dat$H_obs + dat$sigma_H
stopifnot(all(is.finite(dat$H_obs)), all(dat$H_lower > 0),
          all(dat$H_upper > dat$H_lower))

# 3. Fit the Bayesian calibration model.

# Point-gauge uncertainty enters the rating curve through latent H_true.
# Its uniform bounds use dat$sigma_H; no shared datum error is modeled.

dir.create(output_dir, showWarnings=FALSE, recursive=TRUE)
fit_seed <- 1815795276
flume_model <- ulam(
  alist(
    vector[N]:H_true ~ dunif(H_lower, H_upper),  # point-gauge stage uncertainty

    # Latent true discharge with per-flume process scatter (vectorized)
    vector[N]:Q_true ~ dnorm(mu , sigma_proc[F]),
    
    # Vectorized mean by flume
    mu <- a[F] * H_true[i]^b[F],
    
    # Observation model: heteroscedastic SD
    # empirical proportional SD + per-flume additive uncertainty attributable
    # to infiltration and potentially other unresolved measurement effects
    Q_obs ~ dnorm(
      Q_true,
      sqrt(
        square( cv_Q * Q_true ) +          # ModMAG-informed proportional SD
          square( sigma_Q_add[F] )         # lab infiltration-related additive SD
      )
    ),
    
    # Per-flume rating-curve parameters
    vector[J]:a ~ dlnorm( log(a_prior_median) , rating_prior_log_sd ),
    vector[J]:b ~ dlnorm( log(b_prior_median) , rating_prior_log_sd ),
    
    # Per-flume process scatter
    vector[J]:sigma_proc ~ dexp(sigma_proc_prior_rate),
    
    # Measurement-error parameter priors
    # Flow meter: empirical proportional SD, shared by all flumes.
    cv_Q ~ dlnorm( log(cv_Q_prior_median) , cv_Q_prior_log_sd ),
    # Flume infiltration: symmetric additive SD, estimated separately by flume.
    vector[J]:sigma_Q_add ~ dlnorm( log(sigma_Q_add_prior_median) , sigma_Q_add_prior_log_sd )
  ),
  data = dat[setdiff(names(dat), c("stage_units", "flow_units"))],
  constraints = list(H_true = "lower=H_lower,upper=H_upper"),
  chains = 4, cores = 4, iter = 4000, seed = fit_seed
)

# Preserve complete draws and diagnostics for later exports without refitting.
attr(flume_model, "cstanfit")$save_object(
  output_file("flume_model_3_cmdstanfit.rds")
)
chain_diagnostics <- export_stan_diagnostics(flume_model, output_dir,
                                            paste0(output_prefix, "_flume_model_3"), dat)

# model convergence check
traceplot_ulam(flume_model, ask=F, window=c(2100,4000))
trankplot(flume_model, ask=F, window=c(2100,4000))
dashboard_flume(flume_model)
save_dashboard(flume_model, output_file("flume_model_3_dashboard.png"))

# Summarize the model and save output
model_summary <- metric_precis(flume_model)
write_metric_precis(model_summary, output_file("flume_model_3_precis.csv"))
model_summary

# plot model summary and save as jpg
jpeg(output_file("flume_model_3_precis.jpg"), width = 6.5, height = 6.5, units = "in", res = 600)
plot(model_summary)
dev.off()
par(mfrow=c(1,1))
plot(model_summary)

post <- extract.samples(flume_model)
# Reproduce the posterior-draw selections used in the reporting products.
set.seed(fit_seed)

# 4. Export rating-curve, point, and uncertainty-component figures.
# Published metric equation: Q (L/s) = a_m * H (m)^b_m.
# Openchannelflow discharge table, pp. 1-2; fixed reference even if priors change.
manuf <- data.frame(
  a_m   = 941.2,
  b_m   = 2.58
)

# Generate the plot
g <- plot_flume_posterior(
  post  = post,               # posterior samples from ulam
  dat   = dat,                # data list used for the model
  prob  = 0.95,               # credible interval probability
  ndraws = nrow(post$a),       # all draws, shared with point-only summaries
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
  post, dat, manuf, output_dir, prefix = output_prefix, input_csv = input_csv,
  fit_file = output_file("flume_model_3_cmdstanfit.rds")
)
print(point_products$plot)
print(point_products$departures)

# Permanent predictive-uncertainty components and variance shares.
# The helper only defines functions; use the same joint draws as all other plots.
uncertainty_products <- save_flume_uncertainty_components(
  post, dat, manuf, output_dir, prefix = output_prefix, input_csv = input_csv,
  fit_file = output_file("flume_model_3_cmdstanfit.rds")
)
print(uncertainty_products$sd)
print(uncertainty_products$shares)

# 5. Export stage-discharge predictions and manuscript Table 4.
# Report stage in metres at 0.01-m increments, retaining the physical extent.
depth_grid_m <- seq.int(0L, 24L) / 100

res <- generate_flume_table(
  post = post,
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
table4 <- generate_manuscript_table4(post, prob = 0.95)
write_manuscript_csv(table4$display, output_file("manuscript_table4.csv"))
write.csv(table4$numeric, output_file("manuscript_table4_numeric.csv"),
          row.names = FALSE)

