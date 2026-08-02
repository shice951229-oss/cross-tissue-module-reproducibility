########################################
# Script: 12_finalize_registry_and_metadata.R
# Purpose: Finalize preprocessing registry and environment/code manifests.
# Input: Locked registry, installed package metadata, and analysis scripts.
# Output: Final registry, package versions, and code manifest.
# Software: R
# Version: 4.5.2
# Random seed: Not applicable (deterministic)
# Author: Study authors
########################################
#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 200)
root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
out <- file.path(root, "13_methodological_reanalysis_v9")
reg_path <- file.path(out, "dataset_preprocessing_registry.csv")
reg <- read.csv(reg_path, check.names = FALSE)

strict_files <- c(
  GSE98793 = "MDD_GSE98793.rds",
  GSE48350 = "AD_GSE48350.rds",
  GSE5281 = "AD_GSE5281.rds",
  GSE33000 = "AD_GSE33000.rds"
)
for (ds in names(strict_files)) {
  x <- readRDS(file.path(out, "preprocessed_strict_mapping", strict_files[[ds]]))
  i <- match(ds, reg$GEO_accession)
  reg$final_sample_n[i] <- ncol(x)
  reg$final_gene_n[i] <- nrow(x)
  reg$mapping_strategy[i] <- "final main: one-to-one official gene symbols only; ambiguous/unmapped probes removed; highest cross-sample IQR probe retained per gene"
}
strict_note <- "Strict mapping was promoted to the final main analysis because the prespecified stability rule failed in GSE33000; legacy first-symbol/highest-mean mapping is sensitivity only."
reg$notes <- vapply(strsplit(reg$notes, strict_note, fixed = TRUE), function(z) trimws(z[1]), character(1))
i33000 <- match("GSE33000", reg$GEO_accession)
reg$final_sample_n[i33000] <- 467L
reg$removed_samples_and_reason[i33000] <- "157 Huntington disease samples excluded from AD-control modeling (624 downloaded/processed; 467 analyzed)"
reg$notes[i33000] <- paste(reg$notes[i33000], strict_note)
reg$additional_quantile_normalization_in_main <- FALSE
write.csv(reg, reg_path, row.names = FALSE, na = "")

writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
ip <- as.data.frame(installed.packages()[, c("Package", "Version", "LibPath", "Built")], stringsAsFactors = FALSE)
ip <- ip[order(tolower(ip$Package)), ]
write.csv(ip, file.path(out, "package_versions.csv"), row.names = FALSE)

writeLines(c(
  "20260802 | matched random-module negative control, AD | B=2000 per module",
  "20260802 | matched random-module negative control, MDD blood and brain | B=2000 per module and tissue",
  "20260802 | WGCNA modulePreservation | 500 permutations per external liver test",
  "GSVA/ssGSEA, mixed models, linear models, Welch tests, BH adjustment, and figure generation were deterministic for fixed inputs.",
  "No stochastic sample deletion, preprocessing selection, probe-mapping selection, or result-dependent model selection was used."
), file.path(out, "random_seeds.txt"))

scripts <- list.files(file.path(out, "scripts"), full.names = TRUE)
scripts <- scripts[file.info(scripts)$isdir %in% FALSE]
workbook_builder <- file.path(out, "workbook_build", "build_workbooks.mjs")
if (file.exists(workbook_builder)) scripts <- c(scripts, workbook_builder)
info <- file.info(scripts)
sha <- vapply(scripts, function(f) digest::digest(file = f, algo = "sha256", serialize = FALSE), character(1))
relative_paths <- sub(paste0("^", out, "/?"), "13_methodological_reanalysis_v9/", normalizePath(scripts, winslash = "/"), fixed = FALSE)
code <- data.frame(
  script = basename(scripts), relative_path = relative_paths,
  language = ifelse(grepl("\\.R$", scripts, ignore.case = TRUE), "R",
                    ifelse(grepl("\\.py$", scripts, ignore.case = TRUE), "Python",
                           ifelse(grepl("\\.ps1$", scripts, ignore.case = TRUE), "PowerShell", "other"))),
  size_bytes = info$size, last_write_time = format(info$mtime, "%Y-%m-%d %H:%M:%S %z"),
  sha256 = sha, stringsAsFactors = FALSE
)
write.csv(code, file.path(out, "code_manifest.csv"), row.names = FALSE)
