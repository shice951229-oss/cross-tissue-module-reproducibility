########################################
# Script: 10_hpa_liver_specificity.R
# Purpose: Describe HPA liver-expression attributes of the 24 fixed modules.
# Input: HPA tissue consensus file and fixed module definitions.
# Output: Module liver-specificity table and figure.
# Software: R
# Version: 4.5.2
# Random seed: Not applicable (deterministic)
# Author: Study authors
########################################
#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 180)
suppressPackageStartupMessages({library(ggplot2); library(patchwork)})
root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
out <- file.path(root, "13_methodological_reanalysis_v9")
hpa_zip <- file.path(root, "08_donor_aware_reanalysis", "external_sources", "HPA", "rna_tissue_consensus.tsv.zip")
if (!file.exists(hpa_zip)) stop("HPA consensus archive unavailable")

modules <- readRDS(file.path(root, "04_intermediate", "gene_sets", "liver_modules_combined.rds"))
modules <- lapply(modules, function(x) unique(toupper(trimws(as.character(x)))))
x <- read.delim(unz(hpa_zip, "rna_tissue_consensus.tsv"), check.names = FALSE)
x$symbol <- toupper(trimws(x[["Gene name"]]))
x$nTPM <- suppressWarnings(as.numeric(x$nTPM))
x <- x[!is.na(x$symbol) & nzchar(x$symbol) & is.finite(x$nTPM), ]

# One row per HPA symbol and tissue; use the maximum if an annotation duplication exists.
wide <- reshape(aggregate(nTPM ~ symbol + Tissue, x, max), idvar = "symbol", timevar = "Tissue", direction = "wide")
rownames(wide) <- wide$symbol
mat <- as.matrix(wide[, setdiff(names(wide), "symbol"), drop = FALSE])
storage.mode(mat) <- "numeric"
colnames(mat) <- sub("^nTPM\\.", "", colnames(mat))
if (!"liver" %in% colnames(mat)) stop("Liver column not found in HPA consensus data")
other <- mat[, colnames(mat) != "liver", drop = FALSE]
liver <- mat[, "liver"]
max_other <- apply(other, 1, max, na.rm = TRUE)
mean_other <- rowMeans(other, na.rm = TRUE)
median_other <- apply(other, 1, median, na.rm = TRUE)
max_all <- apply(mat, 1, max, na.rm = TRUE)
second_all <- apply(mat, 1, function(v) sort(v, decreasing = TRUE, na.last = NA)[2])
top_tissue <- colnames(mat)[max.col(mat, ties.method = "first")]

ann <- data.frame(
  symbol = rownames(mat), liver_nTPM = liver, non_liver_median_nTPM = median_other,
  liver_vs_nonliver_log2_ratio = log2((liver + 1) / (median_other + 1)),
  liver_enriched = liver >= 1 & liver >= 4 * pmax(max_other, .Machine$double.eps),
  liver_enhanced = liver >= 1 & liver >= 4 * pmax(mean_other, .Machine$double.eps),
  tissue_enriched_any = max_all >= 1 & max_all >= 4 * pmax(second_all, .Machine$double.eps),
  top_tissue = top_tissue, stringsAsFactors = FALSE
)
rownames(ann) <- ann$symbol
universe <- ann$symbol

hyper <- function(k, n, K, N) phyper(k - 1, K, N - K, n, lower.tail = FALSE)
rows <- lapply(names(modules), function(m) {
  g <- intersect(modules[[m]], universe); a <- ann[g, , drop = FALSE]
  data.frame(
    module = m, original_gene_n = length(modules[[m]]), HPA_gene_n = length(g),
    HPA_coverage_fraction = length(g) / length(modules[[m]]),
    liver_enriched_n = sum(a$liver_enriched), liver_enriched_fraction = mean(a$liver_enriched),
    liver_enhanced_n = sum(a$liver_enhanced), liver_enhanced_fraction = mean(a$liver_enhanced),
    tissue_enriched_any_n = sum(a$tissue_enriched_any), tissue_enriched_any_fraction = mean(a$tissue_enriched_any),
    median_liver_nTPM = median(a$liver_nTPM),
    median_liver_vs_nonliver_log2_ratio = median(a$liver_vs_nonliver_log2_ratio),
    P_liver_enriched = hyper(sum(a$liver_enriched), nrow(a), sum(ann$liver_enriched), nrow(ann)),
    P_liver_enhanced = hyper(sum(a$liver_enhanced), nrow(a), sum(ann$liver_enhanced), nrow(ann)),
    P_tissue_enriched_any = hyper(sum(a$tissue_enriched_any), nrow(a), sum(ann$tissue_enriched_any), nrow(ann)),
    specificity_definition = "log2((liver nTPM + 1)/(median non-liver tissue nTPM + 1))",
    enrichment_definition = ">=4-fold versus every other tissue; enhanced >=4-fold versus non-liver mean",
    background_universe = paste0(nrow(ann), " HPA consensus genes"), stringsAsFactors = FALSE
  )
})
res <- do.call(rbind, rows)
res$BH_FDR_liver_enriched <- p.adjust(res$P_liver_enriched, "BH")
res$BH_FDR_liver_enhanced <- p.adjust(res$P_liver_enhanced, "BH")
res$BH_FDR_tissue_enriched_any <- p.adjust(res$P_tissue_enriched_any, "BH")
write.csv(res, file.path(out, "module_liver_specificity_HPA_GTEx.csv"), row.names = FALSE, na = "")

res$module_label <- sub("^(hepatic_|WGCNA_)", "", res$module)
res$module_label <- gsub("_module$|_", " ", res$module_label)
res$module_label <- factor(res$module_label, levels = res$module_label[order(res$median_liver_vs_nonliver_log2_ratio)])
p1 <- ggplot(res, aes(median_liver_vs_nonliver_log2_ratio, module_label)) +
  geom_vline(xintercept = 0, color = "#BDBDBD") + geom_point(aes(color = BH_FDR_liver_enhanced < .05), size = 2) +
  scale_color_manual(values = c(`FALSE` = "#767676", `TRUE` = "#D24B40"), name = "Liver-enhanced\nBH-FDR <0.05") +
  theme_classic(base_size = 7, base_family = "Arial") +
  labs(x = "Median liver vs non-liver expression specificity", y = NULL,
       title = "HPA consensus tissue specificity of the 24 modules")
p2 <- ggplot(res, aes(liver_enriched_fraction, module_label)) +
  geom_point(aes(size = HPA_gene_n), color = "#3182BD") + scale_x_continuous(labels = scales::percent) +
  theme_classic(base_size = 7, base_family = "Arial") +
  labs(x = "Liver-enriched genes", y = NULL, size = "HPA genes")
p <- p1 + p2 + patchwork::plot_layout(widths = c(1.25, 1))
ggsave(file.path(out, "Figure_module_liver_specificity.pdf"), p, width = 180 / 25.4, height = 125 / 25.4,
       device = grDevices::cairo_pdf, family = "Arial")
ggsave(file.path(out, "Figure_module_liver_specificity.png"), p, width = 180 / 25.4, height = 125 / 25.4,
       dpi = 600, device = ragg::agg_png)

writeLines(c(
  "# HPA consensus tissue-specificity audit", "",
  "HPA consensus nTPM values were used as an optional descriptive tissue-specificity analysis.",
  "A liver-enriched gene had liver nTPM at least fourfold higher than every other tissue; a liver-enhanced gene had liver nTPM at least fourfold above the non-liver mean.",
  "The module-level background was all HPA consensus genes. BH correction was applied to 24 modules separately for each enrichment definition.",
  "This analysis describes annotation and expression specificity only. It cannot identify the tissue origin of brain or blood module scores."
), file.path(out, "HPA_liver_specificity_report.md"))
