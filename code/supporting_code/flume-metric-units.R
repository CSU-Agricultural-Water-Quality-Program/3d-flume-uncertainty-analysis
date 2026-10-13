# Unit helpers only. Calibration settings are in code/bayes_flume.R's dat list.

# Exact conversions used only to prepare inputs or inspect historical fits.
# The active metric model and reporting functions do not convert fitted draws.
# International foot: https://www.nist.gov/pml/special-publication-811
FT_TO_M <- 0.3048
CFS_TO_LPS <- FT_TO_M^3 * 1000

assert_flume_metric_units <- function(dat) {
  stopifnot(identical(dat$stage_units, "m"), identical(dat$flow_units, "L/s"))
  invisible(TRUE)
}

convert_legacy_flume_draws <- function(post) {
  stopifnot(is.matrix(post$a), identical(dim(post$a), dim(post$b)))
  # Preserve paired a,b draws when transforming a historical imperial fit.
  post$a <- post$a * CFS_TO_LPS / FT_TO_M^post$b
  for (name in intersect(c("Q_true", "sigma_proc", "sigma_Q_add"), names(post))) {
    post[[name]] <- post[[name]] * CFS_TO_LPS
  }
  if (!is.null(post$H_true)) post$H_true <- post$H_true * FT_TO_M
  post
}
