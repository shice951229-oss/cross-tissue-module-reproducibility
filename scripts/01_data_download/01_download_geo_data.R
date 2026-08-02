########################################
# Script: 01_download_geo_data.R
# Purpose: Download public GEO series matrices and supplementary-file listings into an untracked local staging directory.
# Input: GEO accessions in config.yaml; internet connection.
# Output: Local GEO files under the configured staging directory; no files are committed.
# Software: R
# Version: 4.5.2
# Random seed: Not applicable (deterministic download)
# Author: Study authors
########################################
suppressPackageStartupMessages(library(GEOquery))
args <- commandArgs(trailingOnly = TRUE)
stage <- if (length(args) >= 1L) args[[1L]] else file.path("data", "raw")
dir.create(stage, recursive = TRUE, showWarnings = FALSE)
accessions <- c("GSE135251", "GSE167523", "GSE126848", "GSE98793", "GSE102556",
                "GSE144136", "GSE48350", "GSE5281", "GSE33000")
for (gse in accessions) {
  dest <- file.path(stage, gse)
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  GEOquery::getGEO(gse, GSEMatrix = TRUE, getGPL = FALSE, destdir = dest)
  try(GEOquery::getGEOSuppFiles(gse, makeDirectory = FALSE, baseDir = dest), silent = TRUE)
}
message("Public GEO download completed. These third-party files remain outside version control.")
