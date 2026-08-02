########################################
# Script: 04_score_fixed_modules_ssgsea.R
# Purpose: Apply the locked GSVA/ssGSEA scoring method to a dataset-specific expression matrix.
# Input: Expression RDS (genes x samples), module-definition CSV, output RDS path.
# Output: Dataset-specific ssGSEA score matrix RDS.
# Software: R
# Version: 4.5.2
# Random seed: Not applicable (deterministic)
# Author: Study authors
########################################
suppressPackageStartupMessages({library(GSVA); library(BiocParallel)})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: Rscript 04_score_fixed_modules_ssgsea.R expression.rds module_definitions.csv scores.rds")
expr <- readRDS(args[[1L]])
defs <- read.csv(args[[2L]], stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("module", "gene") %in% names(defs))) stop("module_definitions.csv must contain module and gene columns")
sets <- split(toupper(trimws(defs$gene)), defs$module)
sets <- lapply(sets, function(x) intersect(unique(x[nzchar(x)]), rownames(expr)))
sets <- sets[lengths(sets) >= 3L]
param <- GSVA::ssgseaParam(as.matrix(expr), sets, minSize = 3, normalize = TRUE)
scores <- as.matrix(GSVA::gsva(param, verbose = FALSE, BPPARAM = BiocParallel::SerialParam()))
saveRDS(scores, args[[3L]])
