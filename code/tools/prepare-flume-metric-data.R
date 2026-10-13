# Convert the original October CSV once; preserve its values and row order.
# Run from the repository root: Rscript code/tools/prepare-flume-metric-data.R
script_arg <- grep("^--file=", commandArgs(), value=TRUE)
repo <- if (length(script_arg) == 1L) {
  dirname(dirname(dirname(normalizePath(sub("^--file=", "", script_arg), winslash="/"))))
} else normalizePath(getwd(), winslash="/")
source(file.path(repo,"code/supporting_code/flume-metric-units.R"))
source_csv <- file.path(repo,"data/New Data Oct 2026/flume_data_oct2026.csv")
metric_csv <- file.path(repo,"data/New Data Oct 2026/flume_data_oct2026_metric.csv")
source_hash <- tools::md5sum(source_csv)
original <- read.csv(source_csv,fileEncoding="UTF-8-BOM",stringsAsFactors=FALSE)
stopifnot(identical(names(original),c("depth_ft","flow_cfs","flume_type")),
          all(is.finite(original$depth_ft)), all(original$depth_ft > 0),
          all(is.finite(original$flow_cfs)), all(original$flow_cfs >= 0),
          !anyNA(original$flume_type),
          setequal(original$flume_type,c("3D1","3D2","FIBER")))
metric <- data.frame(depth_m=original$depth_ft * FT_TO_M,
                     flow_lps=original$flow_cfs * CFS_TO_LPS,
                     flume_type=original$flume_type)
write.csv(metric,metric_csv,row.names=FALSE,fileEncoding="UTF-8")
reopened <- read.csv(metric_csv,stringsAsFactors=FALSE)
stopifnot(identical(original$flume_type,reopened$flume_type),
          isTRUE(all.equal(original$depth_ft,reopened$depth_m/FT_TO_M,tolerance=1e-12)),
          isTRUE(all.equal(original$flow_cfs,reopened$flow_lps/CFS_TO_LPS,tolerance=1e-12)),
          identical(source_hash,tools::md5sum(source_csv)))
manifest <- data.frame(field=c("source_csv","source_md5","metric_csv","metric_md5",
                                "stage_units","flow_units","rows","stage_factor","flow_factor",
                                "conversion_scope"),
                       value=c(normalizePath(source_csv,winslash="/"),unname(source_hash),
                               normalizePath(metric_csv,winslash="/"),unname(tools::md5sum(metric_csv)),
                               "m","L/s",nrow(metric),FT_TO_M,CFS_TO_LPS,
                               "exact unit conversion of existing rounded CSV; no measurements re-estimated"))
write.csv(manifest,sub("\\.csv$","_provenance.csv",metric_csv),row.names=FALSE)
cat(sprintf("Validated %d metric rows; original CSV unchanged.\n%s\n",nrow(metric),metric_csv))
