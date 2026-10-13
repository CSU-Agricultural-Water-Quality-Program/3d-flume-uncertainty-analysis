from pathlib import Path

Path('code/regenerate-flume-point-plots.r').write_text('''# Regenerate metric figures from a saved fit; never compile or sample.
# Optional arguments: input_csv fit_rds output_dir.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) %in% c(0L,3L))
script_arg <- grep("^--file=",commandArgs(),value=TRUE)
repo_dir <- dirname(dirname(normalizePath(sub("^--file=","",script_arg[1]),winslash="/")))
source(file.path(repo_dir,"code/flume-metric-units.r"))
if (!length(args)) args <- c(
  file.path(repo_dir,"data/New Data Oct 2026/flume_data_oct2026_metric.csv"),
  file.path(repo_dir,"output/bayes/oct2026_metric/oct2026_metric_flume_model_3_cmdstanfit.rds"),
  file.path(repo_dir,"output/bayes/oct2026_metric"))
input_csv <- args[1]; fit_file <- args[2]; output_dir <- args[3]
metadata_file <- sub("_cmdstanfit\\\\.rds$","_run_metadata.csv",fit_file)
if (!all(file.exists(c(input_csv,fit_file,metadata_file))))
  stop("A matching saved fit and metadata are required. This driver does not fit a model.")
metadata <- read.csv(metadata_file,stringsAsFactors=FALSE)
value <- function(field,default=NULL) {
  x <- metadata$value[metadata$field==field]
  if (!length(x)) default else x
}
stopifnot(identical(value("input_md5"),unname(tools::md5sum(input_csv))))
legacy <- identical(value("calibration_stage_units"),"ft") &&
          identical(value("calibration_flow_units"),"cfs")
native <- identical(value("calibration_stage_units"),"m") &&
          identical(value("calibration_flow_units"),"L/s")
stopifnot(legacy || native)
if (legacy && normalizePath(output_dir,mustWork=FALSE,winslash="/")==
              normalizePath(dirname(fit_file),winslash="/"))
  stop("Use a separate output directory for metric displays of a historical imperial fit.")
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
existing <- list.files(output_dir,full.names=TRUE)
before <- tools::md5sum(existing[!file.info(existing)$isdir])
library(cmdstanr)
library(ggplot2)
source(file.path(repo_dir,"code/flume-point-plots.r"))
source(file.path(repo_dir,"code/flume-uncertainty-components.r"))
df <- read.csv(input_csv,fileEncoding="UTF-8-BOM")
df <- df[df$flume_type!="EQ",,drop=FALSE]
dat <- list(N=nrow(df),J=length(unique(df$flume_type)),
            F=as.integer(factor(df$flume_type)),stage_units="m",flow_units="L/s")
if (legacy) {
  stopifnot(all(c("depth_ft","flow_cfs") %in% names(df)))
  dat$H_obs <- df$depth_ft*FT_TO_M
  dat$Q_obs <- df$flow_cfs*CFS_TO_LPS
  dat$sigma_H <- as.numeric(value("stage_accuracy_half_width_cm","0.01"))/100
} else {
  stopifnot(all(c("depth_m","flow_lps") %in% names(df)))
  dat$H_obs <- df$depth_m; dat$Q_obs <- df$flow_lps
  dat$sigma_H <- as.numeric(value("stage_accuracy_half_width_m"))
}
stopifnot(length(dat$sigma_H)==1L,is.finite(dat$sigma_H),dat$sigma_H>0)
if (!is.null(value("cv_Q_prior_median_fraction")))
  dat$cv_Q_prior_median <- as.numeric(value("cv_Q_prior_median_fraction"))
fit <- readRDS(fit_file)
has_latent_stage <- "H_true" %in% fit$metadata()$stan_variables
draws <- unclass(as.matrix(fit$draws(
  variables=c("a","b","sigma_proc","sigma_Q_add","cv_Q","Q_true",
              if (has_latent_stage) "H_true"),format="draws_matrix")))
extract_vector <- function(name,n) draws[,sprintf("%s[%d]",name,seq_len(n)),drop=FALSE]
post <- list(a=extract_vector("a",dat$J),b=extract_vector("b",dat$J),
  sigma_proc=extract_vector("sigma_proc",dat$J),sigma_Q_add=extract_vector("sigma_Q_add",dat$J),
  cv_Q=draws[,"cv_Q"],Q_true=extract_vector("Q_true",dat$N))
if (has_latent_stage) post$H_true <- extract_vector("H_true",dat$N)
if (legacy) post <- convert_legacy_flume_draws(post)
manufacturer <- if (legacy) data.frame(a_m=1.55*CFS_TO_LPS/FT_TO_M^2.58,b_m=2.58) else
  data.frame(a_m=as.numeric(value("a_prior_median")),b_m=as.numeric(value("b_prior_median")))
products <- save_flume_point_plots(post,dat,manufacturer,output_dir,
                                  input_csv=input_csv,fit_file=fit_file)
components <- save_flume_uncertainty_components(post,dat,manufacturer,output_dir,
                                              input_csv=input_csv,fit_file=fit_file)
exported <- read.csv(file.path(output_dir,paste0(basename(output_dir),"_flume_point_plot_data.csv")))
reference <- vapply(seq_len(dat$N),function(i)
  post$a[,dat$F[i]]*dat$H_obs[i]^post$b[,dat$F[i]],numeric(nrow(post$a)))
stopifnot(nrow(exported)==dat$N,
  isTRUE(all.equal(exported$stage_m,dat$H_obs)),
  isTRUE(all.equal(exported$stage_upper_m-exported$stage_m,rep(dat$sigma_H,dat$N))),
  isTRUE(all.equal(exported$latent_true_flow_lps,unname(apply(post$Q_true,2,median)))),
  isTRUE(all.equal(exported$estimated_flow_lps,unname(apply(reference,2,median)))),
  isTRUE(all.equal(exported$flow_lower_lps,unname(apply(reference,2,quantile,0.025)))),
  isTRUE(all.equal(exported$flow_upper_lps,unname(apply(reference,2,quantile,0.975)))))
expressions <- parse(file.path(repo_dir,"code/bayesflume.r"))
definition <- Filter(function(x) is.call(x) && identical(x[[1]],as.name("<-")) &&
  identical(x[[2]],as.name("plot_flume_posterior")),as.list(expressions))
stopifnot(length(definition)==1L)
eval(definition[[1]])
g <- plot_flume_posterior(post,dat,prob=0.95,ndraws=nrow(draws),
                         show_manufacturer=manufacturer,ribbon_kind="parameter")
curve_points <- g$layers[[3]]$data
stopifnot(isTRUE(all.equal(unname(curve_points$Q),exported$latent_true_flow_lps)),
          isTRUE(all.equal(unname(curve_points$Q_lo),exported$flow_lower_lps)),
          isTRUE(all.equal(unname(curve_points$Q_hi),exported$flow_upper_lps)))
prefix <- basename(output_dir)
ggsave(file.path(output_dir,paste0(prefix,"_flume_posterior_curves_final.jpg")),g,
       width=9,height=6,dpi=600)
g_panels <- g+facet_wrap(vars(flume),nrow=1,labeller=labeller(flume=c(
  "1"="3D Printed 1","2"="3D Printed 2","3"="Fiberglass")))
ggsave(file.path(output_dir,paste0(prefix,"_flume_posterior_curves_by_flume.png")),g_panels,
       width=11,height=5,dpi=300)
ours <- paste0(prefix,"_flume_(points_|point_plot_|posterior_curves_|figure_captions|uncertainty_)")
protected <- before[!grepl(ours,basename(names(before)))]
stopifnot(identical(unname(protected),unname(tools::md5sum(names(protected)))))
cat(sprintf("Validated %d observations and %d saved draws in m and L/s; no sampling.\\n",dat$N,nrow(draws)))
''',encoding='utf-8')
Path('code/preview-flume-uncertainty-components.r').write_text('''# Compatibility alias for saved-fit metric figure regeneration; no fitting.
# Accepts the same optional input_csv fit_rds output_dir arguments.
script_arg <- grep("^--file=",commandArgs(),value=TRUE)
repo <- dirname(dirname(normalizePath(sub("^--file=","",script_arg[1]),winslash="/")))
source(file.path(repo,"code/regenerate-flume-point-plots.r"))
''',encoding='utf-8')
