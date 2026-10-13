library(rethinking)
library(dplyr)
repo <- normalizePath(getwd(), winslash = "/")
files <- c("code/bayesflume.r", "code/flume-point-plots.r",
           "code/flume-uncertainty-components.r", "code/regenerate-flume-point-plots.r")
invisible(lapply(files, parse))
expressions <- parse(files[1])
assignment <- function(x, name) is.call(x) && identical(x[[1]], as.name("<-")) &&
  identical(x[[2]], as.name(name))
for (name in c("FT_TO_CM", "CFS_TO_LPS", "STAGE_ACCURACY_CM",
              "CV_Q_PRIOR_MEDIAN", "CV_Q_PRIOR_LOG_SD", "SIGMA_Q_ADD_PRIOR_MEDIAN_CFS")) {
  eval(expressions[[which(vapply(expressions, assignment, logical(1), name=name))]])
}
df_raw <- read.csv("data/New Data Oct 2026/flume_data_oct2026.csv", fileEncoding="UTF-8-BOM")
df <- df_raw[df_raw$flume_type != "EQ", ]
eval(expressions[[which(vapply(expressions, assignment, logical(1), name="dat"))]])
dat$H_lower <- dat$H_obs - dat$sigma_H
dat$H_upper <- dat$H_obs + dat$sigma_H
dat$cv_Q_prior_median <- CV_Q_PRIOR_MEDIAN
dat$cv_Q_prior_log_sd <- CV_Q_PRIOR_LOG_SD
dat$sigma_Q_add_prior_median <- SIGMA_Q_ADD_PRIOR_MEDIAN_CFS
stopifnot(dat$N == 30, all(dat$H_lower > 0),
          isTRUE(all.equal(dat$sigma_H * FT_TO_CM, 0.15)),
          isTRUE(all.equal((dat$H_upper-dat$H_lower)*FT_TO_CM, rep(0.3,30))))
fit_seed <- 1815795276
model_call <- expressions[[which(vapply(expressions, assignment, logical(1), name="flume_model_3"))]][[3]]
model_call$sample <- FALSE
model <- eval(model_call)
stan_code <- model$model
stan_file <- file.path(repo, "tmp/meter_audit/stage_error_model.stan")
writeLines(stan_code, stan_file)
stopifnot(is.character(stan_code),
          grepl("vector<lower=H_lower,upper=H_upper>\\[N\\] H_true", stan_code),
          grepl("H_true ~ uniform", stan_code, fixed=TRUE),
          grepl("H_true[i]^", stan_code, fixed=TRUE))
stanc <- file.path(cmdstanr::cmdstan_path(), "bin/stanc.exe")
status <- system2(stanc, shQuote(stan_file), stdout=TRUE, stderr=TRUE)
stopifnot(is.null(attr(status,"status")) || attr(status,"status") == 0)
source("code/flume-point-plots.r")
fit <- readRDS("output/bayes/oct2026/oct2026_flume_model_3_cmdstanfit.rds")
stopifnot(!("H_true" %in% fit$metadata()$stan_variables))
draws <- as.matrix(fit$draws(variables=c("a","b","Q_true"), format="draws_matrix"))
post <- list(a=draws[,sprintf("a[%d]",1:3)], b=draws[,sprintf("b[%d]",1:3)],
             Q_true=draws[,sprintf("Q_true[%d]",1:30)])
old_points <- build_flume_point_data(post, dat)
stopifnot(!any(old_points$stage_uncertainty_fitted))
post$H_true <- matrix(rep(dat$H_obs, each=nrow(draws)),ncol=30)
new_points <- build_flume_point_data(post,dat)
stopifnot(all(new_points$stage_uncertainty_fitted),
          isTRUE(all.equal(new_points$stage_cm-new_points$stage_lower_cm,rep(0.15,30))))
write_flume_figure_captions(new_points,data.frame(a_m=1.55,b_m=2.58),
                          "tmp/meter_audit","stage_validation")
captions <- readLines("tmp/meter_audit/stage_validation_flume_figure_captions.md")
stopifnot(any(grepl("fits latent stage",captions,fixed=TRUE)),
          !any(grepl("Stage is fixed",captions,fixed=TRUE)))
cat("PASS: R syntax; 30 positive stage intervals at +/-0.15 cm; uniform latent-stage model; Stan semantic validation. No sampling.\n")
cat("PASS: saved fixed-stage fit remains correctly identified; new stage bounds and caption semantics validated using a reporting fixture.\n")
