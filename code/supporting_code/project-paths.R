# Supporting path resolver; defines functions only and never changes getwd().
flume_project_root <- function(start=getwd()) {
  candidate <- normalizePath(start,winslash="/",mustWork=TRUE)
  if (!dir.exists(candidate)) candidate <- dirname(candidate)
  repeat {
    if (file.exists(file.path(candidate,"code/bayes_flume.R")) &&
        file.exists(file.path(candidate,"code/supporting_code/flume-metric-units.R")) &&
        dir.exists(file.path(candidate,"data"))) return(candidate)
    parent <- dirname(candidate)
    if (identical(parent,candidate)) stop("Cannot locate the flume project from: ",start)
    candidate <- parent
  }
}
