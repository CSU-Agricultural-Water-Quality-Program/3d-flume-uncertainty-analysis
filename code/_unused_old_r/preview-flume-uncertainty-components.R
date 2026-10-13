# Compatibility alias for saved-fit metric figure regeneration; no fitting.
# Accepts the same optional input_csv fit_rds output_dir arguments.
script_arg <- grep("^--file=",commandArgs(),value=TRUE)
repo <- dirname(dirname(dirname(normalizePath(sub("^--file=","",script_arg[1]),winslash="/"))))
source(file.path(repo,"code/tools/regenerate-flume-figures.R"))
