# Compatibility driver: regenerate permanent component figures from saved draws.
# Run: Rscript code/preview-flume-uncertainty-components.r
# Normal use: source bayesflume.r, which now generates these automatically.
library(cmdstanr)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
stopifnot(length(script_arg) == 1L)
repo <- dirname(dirname(normalizePath(sub("^--file=", "", script_arg), winslash = "/")))
source(file.path(repo, "code/flume-uncertainty-components.r"))
output_dir <- file.path(repo, "output/bayes/oct2026")
fit_path <- file.path(output_dir, "oct2026_flume_model_3_cmdstanfit.rds")
input_path <- file.path(repo, "data/New Data Oct 2026/flume_data_oct2026.csv")
metadata <- read.csv(file.path(output_dir, "oct2026_flume_model_3_run_metadata.csv"))
stopifnot(identical(metadata$value[metadata$field == "input_md5"],
                    unname(tools::md5sum(input_path))))
existing <- list.files(output_dir, full.names = TRUE, recursive = FALSE)
existing <- existing[!file.info(existing)$isdir]
protected <- existing[!grepl("^oct2026_flume_uncertainty_", basename(existing))]
before <- tools::md5sum(protected)
df <- read.csv(input_path, fileEncoding = "UTF-8-BOM")
df <- df[df$flume_type != "EQ", , drop = FALSE]
dat <- list(H_obs = df$depth_ft, F = as.integer(factor(df$flume_type)),
            J = length(unique(df$flume_type)))
d <- unclass(as.matrix(readRDS(fit_path)$draws(
  variables = c("a", "b", "sigma_proc", "sigma_Q_add", "cv_Q"), format = "draws_matrix")))
extract_vector <- function(name) d[, sprintf("%s[%d]", name, seq_len(dat$J)), drop = FALSE]
post <- list(a = extract_vector("a"), b = extract_vector("b"),
             sigma_proc = extract_vector("sigma_proc"), sigma_Q_add = extract_vector("sigma_Q_add"),
             cv_Q = d[, "cv_Q"])
products <- save_flume_uncertainty_components(post, dat, data.frame(a_m = 1.55, b_m = 2.58),
              output_dir, prefix = "oct2026", input_csv = input_path, fit_file = fit_path)
stopifnot(identical(unname(before), unname(tools::md5sum(protected))))
cat(sprintf("Validated %d flume-stage combinations from %d saved draws. Other outputs unchanged.\nCaptions: %s\n",
            nrow(products$data)/4, nrow(d), products$captions))
