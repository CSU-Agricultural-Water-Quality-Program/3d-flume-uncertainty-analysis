# Read-only reorganization checks: no model sampling or result exports.
options(cmdstanr_no_ver_check=TRUE)
repo <- normalizePath(getwd(),winslash="/")
entry <- file.path(repo,"code/base_flume.R")
Sys.unsetenv("FLUME_RUN_CALIBRATION")
files <- list.files(file.path(repo,"code"),pattern="[.][rR]$",recursive=TRUE,full.names=TRUE)
invisible(lapply(files,parse))
for (folder in c(repo,file.path(repo,"code"),tempdir())) {
  setwd(folder)
  before <- getwd()
  env <- new.env()
  result <- tryCatch(source(entry,local=env),error=function(e) conditionMessage(e))
  stopifnot(is.character(result),grepl("Sampling is disabled",result,fixed=TRUE),
            identical(getwd(),before),env$dat$N==30L,
            env$dat$a_prior_median==941.2,env$dat$b_prior_median==2.58,
            identical(env$dat$stage_units,"m"),identical(env$dat$flow_units,"L/s"),
            identical(env$flume_repo_dir,repo),file.exists(env$input_csv),
            identical(env$output_dir,file.path(repo,"output/bayes/oct2026_metric")))
}
setwd(repo)
# Reproduce RStudio's Run/selection execution: no source frames or --file flag.
expressions <- parse(entry)
for (folder in c(repo,file.path(repo,"code"),file.path(repo,"code/supporting_code"),
                 file.path(repo,"data/New Data Oct 2026"))) {
  setwd(folder)
  before <- getwd()
  e <- new.env()
  e$commandArgs <- function(...) character()
  e$sys.frames <- function() list()
  result <- tryCatch({
    for (expression in expressions) eval(expression,envir=e)
  },error=function(err) conditionMessage(err))
  stopifnot(grepl("Sampling is disabled",result,fixed=TRUE),
            is.null(e$.entry_file),identical(e$flume_repo_dir,repo),
            e$dat$N==30L,identical(getwd(),before),
            !exists("flume_model_3",envir=e,inherits=FALSE))
}
# Outside the project, interactive execution must give an actionable path error.
setwd(tempdir())
e <- new.env()
e$commandArgs <- function(...) character()
e$sys.frames <- function() list()
result <- tryCatch({
  for (expression in expressions) eval(expression,envir=e)
},error=function(err) conditionMessage(err))
stopifnot(grepl("set your working directory",result,fixed=TRUE),
          is.null(e$flume_repo_dir),!exists("dat",envir=e,inherits=FALSE))
setwd(repo)
for (name in c("metric_precis","write_metric_precis","dashboard_flume","save_dashboard",
               "export_stan_diagnostics","plot_flume_posterior","generate_flume_table",
               "generate_manuscript_table4","write_manuscript_csv","save_flume_point_plots",
               "save_flume_uncertainty_components")) stopifnot(is.function(get(name)))
expressions <- parse(entry)
model <- Filter(function(x) is.call(x) && identical(x[[1]],as.name("<-")) &&
                        identical(x[[2]],as.name("flume_model_3")),as.list(expressions))[[1]][[3]]
model$sample <- FALSE
env$fit_seed <- 1815795276
generated <- eval(model,envir=env)
previous <- paste(readLines(file.path(repo,"tmp/metric-validation/native_metric_model.stan")),collapse="\n")
stopifnot(identical(trimws(generated$model),trimws(previous)))

# Check moved standalone tools' bootstrap paths without executing their actions.
real_command_args <- commandArgs
for (relative in c("code/tools/prepare-flume-metric-data.R","code/tools/validate-flume-metric.R",
                   "code/tools/regenerate-flume-figures.R",
                   "code/_unused_old_r/preview-flume-uncertainty-components.R")) {
  file <- file.path(repo,relative)
  e <- new.env()
  e$commandArgs <- function(trailingOnly=FALSE) if (trailingOnly) character() else
    c("Rscript",paste0("--file=",file))
  x <- parse(file)
  for (expression in x) {
    if (is.call(expression) && identical(expression[[1]],as.name("source"))) break
    if (is.call(expression) && identical(expression[[1]],as.name("setwd"))) break
    eval(expression,envir=e)
  }
  candidate <- if (exists("repo_dir",envir=e,inherits=FALSE)) e$repo_dir else e$repo
  stopifnot(identical(candidate,repo))
}
cat("PASS: all R files parse; main resolves paths from repo, code, and external directories.\n")
cat("PASS: interactive execution resolves four working directories; external paths give clear guidance.\n")
cat("PASS: supporting functions load; standalone-tool roots resolve; working directory is preserved.\n")
cat("PASS: generated Stan model is identical to the previously validated metric model.\n")
cat("No sampling, compilation, or scientific output regeneration.\n")
