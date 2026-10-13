# Headless driver: preserve old outputs and capture the complete run log.
options(cmdstanr_no_ver_check=TRUE)
repo <- normalizePath(getwd(),winslash="/")
protected_dirs <- file.path(repo,"output/bayes",c("oct2026","sep2026"))
protected <- unlist(lapply(protected_dirs,list.files,full.names=TRUE,recursive=TRUE))
protected <- protected[!file.info(protected)$isdir]
before <- tools::md5sum(protected)
saveRDS(before,file.path(repo,"tmp/meter_audit/protected_outputs_before_empirical.rds"))
grDevices::pdf(file=NULL)
source(file.path(repo,"code/bayesflume.r"), chdir=FALSE)
grDevices::dev.off()
stopifnot(identical(before,tools::md5sum(names(before))))
saveRDS(list(data=dat,posterior=post3),
        output_file("flume_model_3_data_and_posterior.rds"))
saveRDS(flume_model_3,output_file("flume_model_3_ulam.rds"))
cat("PASS: full calibration and reporting completed; all protected outputs unchanged.\n")
