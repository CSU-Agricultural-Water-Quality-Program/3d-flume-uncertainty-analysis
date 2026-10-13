# Validate native units, model generation/compilation, and reporting without sampling.
# Run from repo root: Rscript code/tools/validate-flume-metric.R
options(cmdstanr_no_ver_check=TRUE)
library(rethinking)
library(ggplot2)
library(dplyr)
library(tidyr)
script_arg <- grep("^--file=",commandArgs(),value=TRUE)
repo <- dirname(dirname(dirname(normalizePath(sub("^--file=","",script_arg[1]),winslash="/"))))
setwd(repo)
source("code/supporting_code/flume-metric-units.R")
validation_dir <- "tmp/metric-validation"
dir.create(validation_dir,recursive=TRUE,showWarnings=FALSE)
protected <- list.files("output/bayes",full.names=TRUE,recursive=TRUE)
protected <- protected[!file.info(protected)$isdir]
before <- tools::md5sum(protected)
original_csv <- "data/New Data Oct 2026/flume_data_oct2026.csv"
original_hash <- tools::md5sum(original_csv)
original <- read.csv(original_csv,fileEncoding="UTF-8-BOM")
input_csv <- "data/New Data Oct 2026/flume_data_oct2026_metric.csv"
df <- read.csv(input_csv)
stopifnot(nrow(df)==30L,identical(df$flume_type,original$flume_type),
          isTRUE(all.equal(df$depth_m/FT_TO_M,original$depth_ft,tolerance=1e-12)),
          isTRUE(all.equal(df$flow_lps/CFS_TO_LPS,original$flow_cfs,tolerance=1e-12)))
files <- c("code/bayes_flume.R","code/supporting_code/flume-metric-units.R","code/tools/prepare-flume-metric-data.R",
           "code/supporting_code/flume-point-plots.R","code/supporting_code/flume-uncertainty-components.R",
           "code/tools/regenerate-flume-figures.R","code/_unused_old_r/preview-flume-uncertainty-components.R",
           "code/supporting_code/flume-reporting.R","code/supporting_code/project-paths.R")
invisible(lapply(files,parse))
expressions <- parse(files[1])
assignment <- function(x,name) is.call(x) && identical(x[[1]],as.name("<-")) &&
  identical(x[[2]],as.name(name))
get_assignment <- function(name) {
  idx <- which(vapply(expressions,assignment,logical(1),name=name))
  stopifnot(length(idx)==1L)
  expressions[[idx]]
}
eval(get_assignment("dat"))
for (x in expressions) if (is.call(x) && identical(x[[1]],as.name("<-")) &&
  is.call(x[[2]]) && identical(x[[2]][[1]],as.name("$")) &&
  identical(x[[2]][[2]],as.name("dat"))) eval(x)
assert_flume_metric_units(dat)
# Use the editable settings in the main script; do not duplicate their values.
prior_fields <- c("a_prior_median", "b_prior_median", "rating_prior_log_sd",
                  "cv_Q_prior_median", "cv_Q_prior_log_sd",
                  "sigma_Q_add_prior_median", "sigma_Q_add_prior_log_sd",
                  "sigma_proc_prior_rate")
prior_settings <- unlist(dat[prior_fields])
stopifnot(length(prior_settings)==length(prior_fields),
          all(is.finite(prior_settings)),all(prior_settings>0),
          length(dat$sigma_H)==1L,is.finite(dat$sigma_H),dat$sigma_H>0)
stopifnot(all(dat$H_lower>0),all(table(dat$F)==10),
          isTRUE(all.equal(dat$H_upper-dat$H_lower,rep(2*dat$sigma_H,dat$N))))
fit_seed <- 1815795276
model_call <- get_assignment("flume_model")[[3]]
model_call$sample <- FALSE
model <- eval(model_call)
stan_file <- file.path(validation_dir,"native_metric_model.stan")
if (!file.exists(stan_file) || !identical(paste(readLines(stan_file),collapse="\n"),
                                         sub("\n$","",model$model)))
  writeLines(model$model,stan_file)
parameter_block <- strsplit(strsplit(model$model,"parameters{",fixed=TRUE)[[1]][2],
                           "model{",fixed=TRUE)[[1]][1]
stopifnot(!grepl("prior",parameter_block),
          grepl("H_true[i]^b[F[i]]",model$model,fixed=TRUE),
          model$data$a_prior_median==dat$a_prior_median,
          model$data$b_prior_median==dat$b_prior_median,
          model$data$sigma_Q_add_prior_median==dat$sigma_Q_add_prior_median)
# Compilation only: do not call $sample() or execute the calibration script.
compiled <- cmdstanr::cmdstan_model(stan_file,compile=TRUE,quiet=TRUE)
stopifnot(file.exists(compiled$exe_file()))

# Import reporting functions without executing the main workflow.
source("code/supporting_code/flume-reporting.R")
source("code/supporting_code/flume-point-plots.R")
source("code/supporting_code/flume-uncertainty-components.R")
saved <- readRDS("output/bayes/oct2026/oct2026_flume_model_3_cmdstanfit.rds")
draws <- unclass(as.matrix(saved$draws(variables=c("a","b","sigma_proc","sigma_Q_add","cv_Q","Q_true"),
                                      format="draws_matrix")))
extract_vector <- function(name,n) draws[,sprintf("%s[%d]",name,seq_len(n)),drop=FALSE]
legacy <- list(a=extract_vector("a",3),b=extract_vector("b",3),
  sigma_proc=extract_vector("sigma_proc",3),sigma_Q_add=extract_vector("sigma_Q_add",3),
  cv_Q=draws[,"cv_Q"],Q_true=extract_vector("Q_true",30))
post <- convert_legacy_flume_draws(legacy)
legacy_flow <- vapply(seq_len(dat$N),function(i)
  legacy$a[,dat$F[i]]*original$depth_ft[i]^legacy$b[,dat$F[i]]*CFS_TO_LPS,numeric(nrow(draws)))
metric_flow <- vapply(seq_len(dat$N),function(i)
  post$a[,dat$F[i]]*dat$H_obs[i]^post$b[,dat$F[i]],numeric(nrow(draws)))
stopifnot(isTRUE(all.equal(metric_flow,legacy_flow,tolerance=1e-12)))
manufacturer <- data.frame(a_m=941.2,b_m=2.58) # published reference, independent of edited priors
points <- build_flume_point_data(post,dat)
stopifnot(isTRUE(all.equal(points$stage_m,df$depth_m)),
          isTRUE(all.equal(points$stage_upper_m-points$stage_m,rep(dat$sigma_H,dat$N))),
          isTRUE(all.equal(unname(points$estimated_flow_lps),unname(apply(legacy_flow,2,median)))))
table <- generate_flume_table(post,ndraws=nrow(draws),manuf=manufacturer)
stopifnot(identical(table$numeric$Depth_m,seq.int(0L,24L)/100),
          isTRUE(all.equal(table$numeric$Manufacturer_lps,
                           941.2*table$numeric$Depth_m^2.58)),
          abs(table$numeric$Manufacturer_lps[11]-2.475608233971933)<1e-12)
table4 <- generate_manuscript_table4(post)
stopifnot(nrow(table4$numeric)==15L,
          all(table4$numeric$units[table4$numeric$parameter=="a"]=="(L/s)/m^b"))
for (j in 1:3) for (parameter in c("a","b","sigma_proc","sigma_Q_add","cv_Q")) {
  row <- table4$numeric[table4$numeric$parameter==parameter,][j,]
  values <- if (parameter=="cv_Q") post$cv_Q*100 else post[[parameter]][,j]
  stopifnot(abs(row$posterior_mean-mean(values))<1e-10,
    isTRUE(all.equal(unname(c(row$hpdi_lower,row$hpdi_upper)),
                     unname(as.numeric(rethinking::HPDI(values,prob=0.95))))))
}
write_manuscript_csv(table4$numeric,file.path(validation_dir,"table4_reporting_fixture.csv"))
budget <- build_flume_uncertainty_components(post,dat,manufacturer)
shares <- aggregate(variance_share_pct~flume+stage_m,budget,sum)
stopifnot(max(abs(shares$variance_share_pct-100))<1e-10,
          min(budget$stage_m)>=min(df$depth_m),max(budget$stage_m)<=max(df$depth_m))
g <- plot_flume_posterior(post,dat,prob=0.95,ndraws=nrow(draws),
                         show_manufacturer=manufacturer,ribbon_kind="parameter")
stopifnot(identical(g$labels$x,"Stage (m)"),
          isTRUE(all.equal(unname(g$layers[[3]]$data$Q),unname(points$latent_true_flow_lps))))
# Validation previews use historical draws converted to metric; they are not a new fit.
post$H_true <- matrix(rep(dat$H_obs,each=nrow(draws)),ncol=30)
points <- build_flume_point_data(post,dat)
ggsave(file.path(validation_dir,"point_units.png"),plot_flume_points(points,manufacturer),
       width=11,height=6,dpi=120,bg="white")
ggsave(file.path(validation_dir,"curve_units.png"),g,width=9,height=6,dpi=120,bg="white")
component_plots <- plot_flume_uncertainty_components(budget)
ggsave(file.path(validation_dir,"component_units.png"),component_plots$sd,
       width=11,height=7,dpi=120,bg="white")
write_flume_figure_captions(points,manufacturer,validation_dir,"metric_validation")
captions <- readLines(file.path(validation_dir,"metric_validation_flume_figure_captions.md"))
stopifnot(any(grepl("941.2 H (m)^2.58",captions,fixed=TRUE)),
          !any(grepl("Stage is fixed in the fitted model",captions,fixed=TRUE)),
          identical(before,tools::md5sum(names(before))),
          identical(original_hash,tools::md5sum(original_csv)))
report <- c("PASS: 30 metric inputs round-trip; source CSV unchanged.",
            "PASS: fixed metric prior centers/scales; no prior hyperparameters estimated.",
            "PASS: active model C++ compilation.",
            "PASS: draw-level physical equivalence; native table/plot values; 15 Table 4 summaries.",
            "PASS: metric axes/captions; component variances; historical output hashes unchanged.",
            "No posterior sampling. Validation previews use converted historical draws.")
writeLines(report,file.path(validation_dir,"validation_report.txt"))
cat(paste(report,collapse="\n"),"\n")
