# Regenerate point, curve, and uncertainty-component figures without fitting.
# Run from any directory: Rscript code/regenerate-flume-point-plots.r
# Optional positional arguments: input_csv fit_rds output_dir.
# Uses a saved model 3 fit and checks its recorded input hash before plotting.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) %in% c(0L, 3L))
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/")
repo_dir <- dirname(dirname(script_path))
if (!length(args)) {
  args <- c(file.path(repo_dir, "data/New Data Oct 2026/flume_data_oct2026.csv"),
            file.path(repo_dir, "output/bayes/oct2026/oct2026_flume_model_3_cmdstanfit.rds"),
            file.path(repo_dir, "output/bayes/oct2026"))
}
input_csv <- args[1]
fit_file <- args[2]
output_dir <- args[3]
metadata_file <- sub("_cmdstanfit\\.rds$", "_run_metadata.csv", fit_file)
stopifnot(file.exists(input_csv), file.exists(fit_file), file.exists(metadata_file))
metadata <- read.csv(metadata_file)
stopifnot(identical(metadata$value[metadata$field == "input_md5"],
                    unname(tools::md5sum(input_csv))))
# Hash every existing result to mechanically check that generation is additive.
existing <- list.files(output_dir, full.names = TRUE)
before <- tools::md5sum(existing[file.info(existing)$isdir %in% FALSE])
library(cmdstanr)
library(ggplot2)
source(file.path(repo_dir, "code/flume-point-plots.r"))
source(file.path(repo_dir, "code/flume-uncertainty-components.r"))
df <- read.csv(input_csv, fileEncoding = "UTF-8-BOM")
df <- df[df$flume_type != "EQ", , drop = FALSE]
stopifnot(identical(sort(unique(df$flume_type)), c("3D1", "3D2", "FIBER")))
dat <- list(N = nrow(df), J = length(unique(df$flume_type)),
            Q_obs = df$flow_cfs, H_obs = df$depth_ft,
            F = as.integer(factor(df$flume_type)), sigma_H = 0.01 / 30.48)
fit <- readRDS(fit_file)
draws <- unclass(as.matrix(fit$draws(
  variables = c("a", "b", "sigma_proc", "sigma_Q_add", "cv_Q", "Q_true"),
  format = "draws_matrix"
)))
extract_vector <- function(name, n) {
  expected <- sprintf("%s[%d]", name, seq_len(n))
  stopifnot(setequal(grep(paste0("^", name, "\\["), colnames(draws), value = TRUE), expected))
  draws[, expected, drop = FALSE]
}
post <- list(a = extract_vector("a", dat$J), b = extract_vector("b", dat$J),
             sigma_proc = extract_vector("sigma_proc", dat$J),
             sigma_Q_add = extract_vector("sigma_Q_add", dat$J),
             cv_Q = draws[, "cv_Q"], Q_true = extract_vector("Q_true", dat$N))
products <- save_flume_point_plots(post, dat, data.frame(a_m = 1.55, b_m = 2.58),
                                  output_dir, input_csv = input_csv, fit_file = fit_file)
uncertainty_products <- save_flume_uncertainty_components(
  post, dat, data.frame(a_m = 1.55, b_m = 2.58), output_dir,
  input_csv = input_csv, fit_file = fit_file
)
# Independently check all exported medians and interval endpoints against draws.
csv_file <- file.path(output_dir, paste0(basename(output_dir), "_flume_point_plot_data.csv"))
exported <- read.csv(csv_file)
# Compute the reference independently from paired a,b columns in the saved fit.
metric_draws <- vapply(seq_len(nrow(df)), function(i) {
  draws[, sprintf("a[%d]", dat$F[i])] *
    df$depth_ft[i]^draws[, sprintf("b[%d]", dat$F[i])] * (0.3048^3 * 1000)
}, numeric(nrow(draws)))
stopifnot(nrow(exported) == nrow(df),
          identical(exported$observation_index, seq_len(nrow(df))),
          all(exported$posterior_draws == nrow(draws)),
          all(exported$interval_kind == "parameter"),
          isTRUE(all.equal(exported$latent_true_flow_lps,
                           unname(apply(post$Q_true * (0.3048^3 * 1000), 2, median)))),
          isTRUE(all.equal(exported$estimated_flow_lps, unname(apply(metric_draws, 2, median)))),
          isTRUE(all.equal(exported$flow_lower_lps, unname(apply(metric_draws, 2, quantile, 0.025)))),
          isTRUE(all.equal(exported$flow_upper_lps, unname(apply(metric_draws, 2, quantile, 0.975)))),
          isTRUE(all.equal(exported$stage_cm, df$depth_ft * 30.48)),
          isTRUE(all.equal(exported$stage_upper_cm - exported$stage_cm, rep(0.01, nrow(df)))))
for (i in seq_len(nrow(df))) {
  reference <- exported$manufacturer_flow_lps[i]
  if (reference > 0) {
    relative <- 100 * (metric_draws[, i] - reference) / reference
    latent_relative <- 100 * (post$Q_true[, i] * (0.3048^3 * 1000) - reference) / reference
    stopifnot(isTRUE(all.equal(exported$relative_difference_pct[i], median(relative))),
              isTRUE(all.equal(exported$latent_relative_difference_pct[i], median(latent_relative))),
              isTRUE(all.equal(exported$relative_lower_pct[i], unname(quantile(relative, 0.025)))),
              isTRUE(all.equal(exported$relative_upper_pct[i], unname(quantile(relative, 0.975)))))
  } else {
    stopifnot(is.na(exported$relative_difference_pct[i]),
              is.na(exported$relative_lower_pct[i]), is.na(exported$relative_upper_pct[i]))
  }
}
# Reuse only the plotting function definition from the main analysis. No model
# fitting, diagnostics, or top-level source side effects are executed here.
expressions <- parse(file.path(repo_dir, "code/bayesflume.r"))
definition <- Filter(function(x) {
  is.call(x) && identical(x[[1]], as.name("<-")) &&
    identical(x[[2]], as.name("plot_flume_posterior"))
}, as.list(expressions))
stopifnot(length(definition) == 1L)
eval(definition[[1]])
FT_TO_CM <- 30.48
CFS_TO_LPS <- 0.3048^3 * 1000
g <- plot_flume_posterior(post, dat, prob = 0.95, ndraws = nrow(draws),
                          show_manufacturer = data.frame(a_m = 1.55, b_m = 2.58),
                          ribbon_kind = "parameter", show_qobs_errors = FALSE)
# Main plot point intervals must equal the new point-only parameter summaries.
curve_points <- g$layers[[3]]$data
stopifnot(isTRUE(all.equal(unname(curve_points$Q), exported$latent_true_flow_lps)),
          isTRUE(all.equal(unname(curve_points$Q_parameter_med), exported$estimated_flow_lps)),
          isTRUE(all.equal(unname(curve_points$Q_lo), exported$flow_lower_lps)),
          isTRUE(all.equal(unname(curve_points$Q_hi), exported$flow_upper_lps)))
prefix <- basename(output_dir)
ggsave(file.path(output_dir, paste0(prefix, "_flume_posterior_curves_final.jpg")),
       g, width = 9, height = 6, dpi = 600)
g_panels <- g + facet_wrap(vars(flume), nrow = 1,
                           labeller = labeller(flume = c("1" = "3D Printed 1",
                             "2" = "3D Printed 2", "3" = "Fiberglass")))
ggsave(file.path(output_dir, paste0(prefix, "_flume_posterior_curves_by_flume.png")),
       g_panels, width = 11, height = 5, dpi = 300)
# Regeneration may update our own products, but must preserve all prior products.
ours <- paste0(prefix, "_flume_(points_|point_plot_|posterior_curves_|figure_captions|uncertainty_)")
protected <- before[!grepl(ours, basename(names(before)))]
stopifnot(identical(unname(protected), unname(tools::md5sum(names(protected)))))
cat(sprintf("Validated %d observations and %d paired parameter draws; point bars match curve intervals; fit, diagnostics, and tables unchanged.\nOutputs: %s\n",
            nrow(df), nrow(draws), normalizePath(output_dir, winslash = "/")))
