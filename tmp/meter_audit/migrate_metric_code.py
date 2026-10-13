from pathlib import Path

ROOT = Path.cwd()
def read(name): return (ROOT/name).read_text(encoding='utf-8-sig')
def write(name,text): (ROOT/name).write_text(text,encoding='utf-8')
def replace(text,old,new,count=None):
    found=text.count(old)
    if not found or (count is not None and found!=count):
        raise RuntimeError(f'Replacement mismatch: {old[:100]!r}, count={found}')
    return text.replace(old,new)

s=read('code/bayesflume.r')
start=s.index('# Calibration inputs and priors remain')
end=s.index('# rethinking 2.42')
s=s[:start]+'''# Native calibration and reporting units: stage in metres, discharge in L/s.
# Prior settings are fixed inputs, not estimated hyperparameters.
flume_repo_dir <- if (file.exists("code/flume-metric-units.r")) {
  normalizePath(".", winslash="/")
} else normalizePath("..", winslash="/")
source(file.path(flume_repo_dir, "code/flume-metric-units.r"))
# Sampling requires an explicit opt-in; syntax/data checks do not sample.
RUN_CALIBRATION <- identical(Sys.getenv("FLUME_RUN_CALIBRATION"), "true")

# Summaries use native draws unchanged; suffixes record units, not conversions.
metric_precis <- function(fit, prob = 0.89) {
  original <- precis(fit, depth = 2, prob = prob)
  if (!is.null(attr(fit, "cstanfit"))) {
    draws <- attr(fit, "cstanfit")$draws(variables=rownames(original), format="draws_array")
  } else {
    draws <- posterior::as_draws_array(as.array(attr(fit, "stanfit")))
    draws <- draws[, , rownames(original), drop=FALSE]
  }
  names <- dimnames(draws)[[3]]
  names <- sub("^a\\\\[", "a_m[", names)
  names <- sub("^alpha\\\\[", "alpha_m[", names)
  names <- sub("^H_true\\\\[", "H_true_m[", names)
  names <- sub("^(Q_true|mu|sigma_proc|sigma_Q_add)(\\\\[|$)", "\\\\1_lps\\\\2", names)
  dimnames(draws)[[3]] <- names
  result <- as.data.frame(posterior::summarise_draws(
    draws, mean=base::mean, sd=stats::sd,
    interval=function(x) stats::quantile(x, probs=c((1-prob)/2,(1+prob)/2)),
    rhat=posterior::rhat, ess_bulk=posterior::ess_bulk))
  rownames(result) <- result$variable
  result$variable <- NULL
  methods::new("precis", result, digits=3)
}

write_metric_precis <- function(summary, filename) {
  result <- as.data.frame(summary)
  v <- rownames(result)
  result$units <- ifelse(grepl("^H_true_m\\\\[",v), "m",
    ifelse(grepl("_lps",v), "L/s",
    ifelse(grepl("^a_m\\\\[",v), "(L/s)/m^b",
    ifelse(grepl("^alpha_m\\\\[",v), "log(a_m)", "dimensionless"))))
  write.csv(result,filename,row.names=TRUE)
}

'''+s[end:]
s=replace(s,'"stage_accuracy_half_width_cm"','"stage_accuracy_half_width_m"')
s=replace(s,'"sigma_Q_add_prior_median_cfs"','"sigma_Q_add_prior_median_lps"')
s=replace(s,'"ft", "cfs", "cm", "L/s"','"m", "L/s", "m", "L/s"')
s=replace(s,'STAGE_ACCURACY_CM','STAGE_ACCURACY_M')
s=replace(s,'SIGMA_Q_ADD_PRIOR_MEDIAN_CFS','SIGMA_Q_ADD_PRIOR_MEDIAN_LPS')
s=replace(s,'"measurement_error_interpretation"),','"measurement_error_interpretation", "a_prior_median", "b_prior_median",\n              "rating_prior_log_sd", "sigma_proc_prior_rate_lps",\n              "sigma_Q_add_prior_log_sd", "manufacturer_equation"),')
s=replace(s,'"empirical proportional plus additive SD; manufacturer velocity term not explicit")',
'''"empirical proportional plus additive SD; manufacturer velocity term not explicit",
              MANUFACTURER_A, MANUFACTURER_B, RATING_PRIOR_LOG_SD,
              SIGMA_PROC_PRIOR_RATE_LPS, SIGMA_Q_ADD_PRIOR_LOG_SD,
              "Q (L/s) = 941.2 * H (m)^2.58")''')
s=replace(s,'setwd("C:/Users/ansle/OneDrive/Documents/GitHub/3d-flume-uncertainty-analysis/data")','setwd(file.path(flume_repo_dir, "data"))')
s=replace(s,'oct2026_empirical_error','oct2026_metric')
s=replace(s,'dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)','')
s=replace(s,'input_csv <- "../data/New Data Oct 2026/flume_data_oct2026.csv"','input_csv <- "New Data Oct 2026/flume_data_oct2026_metric.csv"')
s=replace(s,'# df_raw <- read.csv("../data/flume_data.csv")\n# df_raw <- read.csv("../data/sim_data.csv") # for testing with simulated data\n# input_csv <- "../data/New Data Sep 2026/flume_data_sep2026.csv"',
'''# Historical imperial inputs/fits are preserved separately. Prepare this CSV
# with code/prepare-flume-metric-data.r; do not silently convert observations here.''')
s=replace(s,'Q_obs   = as.numeric(df$flow_cfs),                 # observed discharge (cfs)','Q_obs   = as.numeric(df$flow_lps),                 # observed discharge (L/s)')
s=replace(s,'H_obs   = as.numeric(df$depth_ft),                 # observed head (ft)','H_obs   = as.numeric(df$depth_m),                  # observed head (m)')
s=replace(s,'sigma_H = STAGE_ACCURACY_M / FT_TO_CM             # stage-accuracy half-width, ft; not an SD',
'''sigma_H = STAGE_ACCURACY_M,                       # stage-accuracy half-width, m; not an SD
  stage_units = "m", flow_units = "L/s"''')
s=replace(s,'dat$cv_Q_prior_median <- CV_Q_PRIOR_MEDIAN',
'''assert_flume_metric_units(dat)
stopifnot(all(c("depth_m", "flow_lps", "flume_type") %in% names(df)),
          nrow(df) == 30L, all(table(df$flume_type) == 10L))
dat$a_prior_median <- MANUFACTURER_A
dat$b_prior_median <- MANUFACTURER_B
dat$rating_prior_log_sd <- RATING_PRIOR_LOG_SD
dat$sigma_proc_prior_rate <- SIGMA_PROC_PRIOR_RATE_LPS
dat$sigma_Q_add_prior_log_sd <- SIGMA_Q_ADD_PRIOR_LOG_SD
dat$cv_Q_prior_median <- CV_Q_PRIOR_MEDIAN''')
s=replace(s,'log(1.55)','log(a_prior_median)')
s=replace(s,'log(2.58)','log(b_prior_median)')
s=replace(s,'log(a_prior_median) , 0.5','log(a_prior_median) , rating_prior_log_sd')
s=replace(s,'log(b_prior_median) , 0.5','log(b_prior_median) , rating_prior_log_sd')
s=replace(s,'log(a_prior_median) , 0.5','log(a_prior_median) , rating_prior_log_sd') if 'log(a_prior_median) , 0.5' in s else s
s=replace(s,'~ dexp(1)','~ dexp(sigma_proc_prior_rate)')
s=replace(s,'log(sigma_Q_add_prior_median) , 0.5','log(sigma_Q_add_prior_median) , sigma_Q_add_prior_log_sd')
s=replace(s,'fit_seed <- 1815795276',
'''if (!RUN_CALIBRATION) {
  stop("Metric inputs and model are prepared. Sampling is disabled; explicitly set FLUME_RUN_CALIBRATION=true to fit.")
}
dir.create(output_dir, showWarnings=FALSE, recursive=TRUE)
fit_seed <- 1815795276''')
s=replace(s,'data = dat, constraints = list(H_true = "lower=H_lower,upper=H_upper"),',
'''data = dat[setdiff(names(dat), c("stage_units", "flow_units"))],
  constraints = list(H_true = "lower=H_lower,upper=H_upper"),''')
s=replace(s,'#     cv_Q ~ dlnorm( log(2e-4) , 0.5 ),','#     cv_Q ~ dlnorm( log(cv_Q_prior_median) , cv_Q_prior_log_sd ),')
s=replace(s,'#     vector[J]:sigma_Q_add ~ dlnorm( log(0.02) , 0.5 )',
 '#     vector[J]:sigma_Q_add ~ dlnorm( log(sigma_Q_add_prior_median) , sigma_Q_add_prior_log_sd )')
s=replace(s,'#     # Latent true discharge with per-flume process scatter (vectorized)',
 '#     vector[N]:H_true ~ dunif(H_lower, H_upper),\n#\n#     # Latent true discharge with per-flume process scatter (vectorized)')
s=replace(s,'log(H_obs + 1e-9)','log(H_true[i])')
s=replace(s,'#   data = dat,\n#   chains = 4, cores = 4, iter = 4000',
 '#   data = dat[setdiff(names(dat), c("stage_units", "flow_units"))],\n#   constraints = list(H_true = "lower=H_lower,upper=H_upper"),\n#   chains = 4, cores = 4, iter = 4000')
s=replace(s,'#   dat   <- list(N=nrow(df), Q_obs=df$flow_cfs, H_obs=df$depth_ft,',
 '#   dat   <- list(N=nrow(df), Q_obs=df$flow_lps, H_obs=df$depth_m,')
start=s.index('  # Convert every plotted quantity')
end=s.index('  # Use coordinate limits',start)
s=s[:start]+'  # Inputs, posterior draws, intervals, and manufacturer flow are already metric.\n'+s[end:]
s=replace(s,'max(Hgrid) * FT_TO_CM','max(Hgrid)')
s=replace(s,'stage in cm','stage in m')
s=replace(s,'Stage (cm)','Stage (m)')
s=replace(s,'# Equation: Q (cfs) = a_m * H (ft)^b_m; converted only for reporting.',
 '# Published metric equation: Q (L/s) = a_m * H (m)^b_m.')
s=replace(s,'a_m   = c(1.55)','a_m   = c(MANUFACTURER_A)')
s=replace(s,'b_m   = c(2.58)','b_m   = c(MANUFACTURER_B)')
s=s.replace('depths_cm','depths_m').replace('Depth_cm','Depth_m').replace('depth_grid_cm','depth_grid_m')
s=replace(s,'depths_m = seq.int(0L, 24L)','depths_m = seq.int(0L, 24L) / 100')
s=replace(s,'H <- depths_m / FT_TO_CM  # evaluate the unchanged ft/cfs calibration','H <- depths_m  # evaluate the metric calibration directly')
s=replace(s,'median_lps = md * CFS_TO_LPS','median_lps = md')
s=replace(s,'lo_lps = lo * CFS_TO_LPS, hi_lps = hi * CFS_TO_LPS','lo_lps = lo, hi_lps = hi')
s=s.replace('(df_fmt$Depth_m / FT_TO_CM)','df_fmt$Depth_m').replace('(df_num$Depth_m / FT_TO_CM)','df_num$Depth_m')
s=s.replace(' * CFS_TO_LPS','')
start=s.index('# Evaluate predictions at integer centimetres')
end=s.index('\nres <- generate_flume_table',start)
s=s[:start]+'''# Report stage in metres at 0.01-m increments, retaining the physical extent.
depth_grid_m <- seq.int(0L, 24L) / 100
manuf <- data.frame(a_m=MANUFACTURER_A, b_m=MANUFACTURER_B)
'''+s[end:]
s=replace(s,'H (cm)^b','H (m)^b')
s=replace(s,'(L/s)/cm^b','(L/s)/m^b')
s=replace(s,'a = post$a / FT_TO_CM^post$b','a = post$a')
s=replace(s,'digits <- c(5L, 2L, 3L, 3L, 3L)','digits <- c(2L, 2L, 3L, 3L, 3L)')
# The commented alternate model remains available in the same metric units.
s=s.replace('# If a ~ LogNormal(log(a_prior_median), 0.5), then alpha = log(a) ~ Normal(log(a_prior_median), 0.5)',
 '# With native metric a, alpha = log(a) has the same fixed prior center and width.')
# Always make data units explicit in plotting functions.
s=replace(s,'  ribbon_kind <- match.arg(ribbon_kind)','  assert_flume_metric_units(dat)\n  ribbon_kind <- match.arg(ribbon_kind)')
write('code/bayesflume.r',s)

s=read('code/flume-point-plots.r')
s=replace(s,'  stopifnot(all(c("Q_obs", "H_obs", "F", "sigma_H") %in% names(dat)),',
 '  assert_flume_metric_units(dat)\n  stopifnot(all(c("Q_obs", "H_obs", "F", "sigma_H") %in% names(dat)),')
s=replace(s,'  ft_to_cm <- 30.48\n  cfs_to_lps <- (ft_to_cm / 100)^3 * 1000\n','',count=2)
s=s.replace('stage_cm','stage_m').replace('stage_lower_cm','stage_lower_m').replace('stage_upper_cm','stage_upper_m')
s=s.replace(' * ft_to_cm','').replace(' * cfs_to_lps','')
s=s.replace('h_cm','h_m').replace('(h_m / ft_to_cm)','h_m')
s=replace(s,'width = 0.22','width = 0.0022')
s=replace(s,'Stage (cm)','Stage (m)')
s=replace(s,'  q_m <- manufacturer$a_m * (points$stage_m / 30.48)^manufacturer$b_m *\n    (30.48 / 100)^3 * 1000',
 '  q_m <- manufacturer$a_m * points$stage_m^manufacturer$b_m')
s=replace(s,'sprintf("- Horizontal bars: recorded gauge accuracy +/- %g ft (+/- %g cm); no probability coverage is assigned.", dat$sigma_H, dat$sigma_H * 30.48),',
 'sprintf("- Horizontal bars: recorded gauge accuracy bound +/- %g m; no probability coverage is assigned.", dat$sigma_H),')
s=replace(s,'sprintf("- Manufacturer: Q (cfs) = %g * H (ft)^%g; converted to cm and L/s.", manufacturer$a_m, manufacturer$b_m),',
 'sprintf("- Manufacturer: Q (L/s) = %g * H (m)^%g.", manufacturer$a_m, manufacturer$b_m),')
s=replace(s,'accuracy bounds of +/- %g cm','accuracy bounds of +/- %g m')
s=replace(s,'sprintf("The manufacturer equation is Q (cfs) = %g H (ft)^%g. Figures report stage in cm and flow in L/s.",',
 'sprintf("The manufacturer equation is Q (L/s) = %g H (m)^%g. Calibration and figures use m and L/s.",')
write('code/flume-point-plots.r',s)

s=read('code/flume-uncertainty-components.r')
s=replace(s,'  required <- c("a", "b", "sigma_proc", "sigma_Q_add", "cv_Q")',
 '  assert_flume_metric_units(dat)\n  required <- c("a", "b", "sigma_proc", "sigma_Q_add", "cv_Q")')
s=replace(s,'  conversion <- 0.3048^3 * 1000\n','')
s=s.replace(' * conversion^2','').replace(' * conversion','')
s=s.replace('stage_cm','stage_m').replace('h * 30.48','h')
s=replace(s,'breaks = c(6, 9, 12, 15)','breaks = c(0.06, 0.09, 0.12, 0.15)',count=2)
s=replace(s,'Stage (cm)','Stage (m)')
s=replace(s,'in internal ft/cfs units','in native m and L/s units')
s=replace(s,'Variances are converted to (L/s)^2','Variances are already in (L/s)^2')
write('code/flume-uncertainty-components.r',s)
print('Updated native metric calibration, summaries, plots, predictions, and uncertainty components; no sampling.')
