########################################
# Script: 08b_postprocess_liver_preservation.R
# Purpose: Postprocess locked module-preservation objects without redefining modules.
# Input: Saved preservation objects and module definitions.
# Output: Combined preservation table, size-effect table, report inputs, and figure.
# Software: R
# Version: 4.5.2
# Random seed: Not applicable (postprocessing only)
# Author: Study authors
########################################
#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 180)
suppressPackageStartupMessages({library(ggplot2); library(patchwork)})
root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
out <- file.path(root, "13_methodological_reanalysis_v9")
seed <- 20260802L; n_perm <- 500L
modules_all <- readRDS(file.path(root, "04_intermediate", "gene_sets", "liver_modules_combined.rds"))
modules <- modules_all[grepl("^WGCNA_", names(modules_all))]
modules <- lapply(modules, function(x) unique(toupper(trimws(as.character(x)))))
ref0 <- readRDS(file.path(out, "preprocessed_main_old_mapping", "liver_GSE135251.rds"))
tests <- list(
  GSE126848 = readRDS(file.path(out, "preprocessed_main_old_mapping", "liver_GSE126848.rds")),
  GSE167523 = readRDS(file.path(out, "preprocessed_main_old_mapping", "liver_GSE167523.rds"))
)
ref_var <- apply(ref0, 1, var, na.rm = TRUE)
var_threshold <- quantile(ref_var[ref_var > 0], 0.10, na.rm = TRUE)
ref_eligible <- names(ref_var)[is.finite(ref_var) & ref_var > var_threshold]
module_union <- unique(unlist(modules, use.names = FALSE))

extract_pair <- function(container, test_name) {
  ref_node <- container[[1L]]
  hit <- grep(test_name, names(ref_node), fixed = TRUE)
  if (length(hit)) return(ref_node[[hit[1L]]])
  if (length(ref_node) >= 2L) return(ref_node[[2L]])
  stop("Could not locate preservation results for ", test_name)
}

rows <- lapply(names(tests), function(test_name) {
  rds <- file.path(out, paste0("liver_module_preservation_", test_name, "_full.rds"))
  if (!file.exists(rds)) stop("Missing formal preservation object: ", rds)
  mp <- readRDS(rds)
  ztab <- as.data.frame(extract_pair(mp$preservation$Z, test_name))
  mtab <- as.data.frame(extract_pair(mp$preservation$observed, test_name))
  ztab$color <- rownames(ztab); mtab$color <- rownames(mtab)
  z_col <- grep("^Zsummary", names(ztab), value = TRUE)[1L]
  rank_col <- grep("^medianRank", names(mtab), value = TRUE)[1L]
  common <- Reduce(intersect, list(ref_eligible, rownames(tests[[test_name]]), module_union))
  color <- rep("grey", length(common)); names(color) <- common
  for (m in names(modules)) color[intersect(common, modules[[m]])] <- sub("^WGCNA_", "", m)
  color <- color[color != "grey"]; counts <- table(color)
  d <- data.frame(
    module = paste0("WGCNA_", names(counts)), test_dataset = test_name, reference_dataset = "GSE135251",
    common_gene_n = as.integer(counts),
    reference_original_gene_n = vapply(paste0("WGCNA_", names(counts)), function(m) length(modules[[m]]), integer(1)),
    stringsAsFactors = FALSE
  )
  d$coverage_fraction <- d$common_gene_n / d$reference_original_gene_n
  d$Zsummary <- ztab[[z_col]][match(sub("^WGCNA_", "", d$module), ztab$color)]
  d$medianRank <- mtab[[rank_col]][match(sub("^WGCNA_", "", d$module), mtab$color)]
  d$common_gene_adequacy <- ifelse(d$common_gene_n >= 20L, "adequate for interpretation", "too few; descriptive only")
  d$preservation_class <- as.character(cut(d$Zsummary, breaks = c(-Inf, 2, 10, Inf),
                                           labels = c("no clear evidence", "moderate", "strong"), right = FALSE))
  d$preservation_class[d$common_gene_n < 20L] <- "low common-gene count"
  d$network_type <- "unsigned"
  d$reference_gene_filter <- paste0("variance > 10th percentile of positive variance (", signif(var_threshold, 7), ")")
  d$permutations <- n_perm; d$seed <- seed
  d
})
res <- do.call(rbind, rows); res <- res[order(res$test_dataset, res$medianRank), ]
write.csv(res, file.path(out, "liver_module_preservation_results.csv"), row.names = FALSE, na = "")
size_cor <- do.call(rbind, lapply(split(res, res$test_dataset), function(d) {
  ct <- cor.test(log10(d$common_gene_n), d$Zsummary, method = "spearman", exact = FALSE)
  data.frame(test_dataset = d$test_dataset[1], spearman_rho_log_size_vs_Zsummary = unname(ct$estimate), P = ct$p.value)
}))
write.csv(size_cor, file.path(out, "liver_module_preservation_size_effect.csv"), row.names = FALSE)

res$module_label <- sub("^WGCNA_", "", res$module)
pz <- ggplot(res, aes(Zsummary, reorder(module_label, Zsummary))) +
  geom_vline(xintercept = 2, linetype = "dashed", color = "#767676") +
  geom_vline(xintercept = 10, linetype = "dotted", color = "#4D4D4D") +
  geom_point(aes(size = common_gene_n, color = preservation_class), alpha = .9) + facet_wrap(~test_dataset, scales = "free_y") +
  scale_color_manual(values = c("no clear evidence" = "#BDBDBD", moderate = "#3182BD", strong = "#D24B40", "low common-gene count" = "#7B3294")) +
  scale_size_continuous(range = c(1.5, 5)) + theme_classic(base_size = 7, base_family = "Arial") +
  theme(legend.position = "bottom") + labs(x = "WGCNA preservation Zsummary", y = NULL, color = "Evidence", size = "Common genes",
                                           title = "A  Preservation strength")
pr <- ggplot(res, aes(medianRank, reorder(module_label, -medianRank))) +
  geom_point(aes(size = common_gene_n), color = "#4D4D4D", alpha = .85) + facet_wrap(~test_dataset, scales = "free_y") +
  scale_size_continuous(range = c(1.5, 5)) + theme_classic(base_size = 7, base_family = "Arial") +
  theme(legend.position = "bottom") + labs(x = "Preservation medianRank (lower is better)", y = NULL, size = "Common genes",
                                           title = "B  Size-aware relative rank")
p <- pz / pr + plot_annotation(title = "External preservation of GSE135251 WGCNA modules")
ggsave(file.path(out, "Figure_liver_module_preservation.pdf"), p, width = 180/25.4, height = 200/25.4,
       device = grDevices::cairo_pdf, family = "Arial")
ggsave(file.path(out, "Figure_liver_module_preservation.png"), p, width = 180/25.4, height = 200/25.4,
       dpi = 600, device = ragg::agg_png)

summary_line <- function(ds) {
  ok <- res$test_dataset == ds & res$common_gene_n >= 20L
  s <- sum(ok & res$Zsummary >= 10); m <- sum(ok & res$Zsummary >= 2 & res$Zsummary < 10)
  n <- sum(ok & res$Zsummary < 2); low <- sum(res$test_dataset == ds & res$common_gene_n < 20L)
  paste0("- ", ds, ": ", s, " strong, ", m, " moderate, ", n,
         " with no clear evidence, and ", low, " low-common-gene modules retained only descriptively.")
}
writeLines(c(
  "# External liver module-preservation report", "", paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")), "",
  "## Design", "", "- Reference: GSE135251; tests: GSE126848 and GSE167523.",
  "- Sixteen pre-defined WGCNA modules; unsigned WGCNA modulePreservation; 500 permutations per test; seed 20260802.",
  "- Zsummary is interpreted jointly with medianRank and common module size. Modules with fewer than 20 common genes are descriptive only.", "",
  "## Results", "", vapply(names(tests), summary_line, character(1)), "",
  "## Inferential boundary", "",
  "Preservation supports recurrence of within-module co-expression structure in external liver data. It does not establish the tissue origin of projected brain or blood scores, circulation-mediated transfer, mediation, mechanism, or causality."
), file.path(out, "liver_module_preservation_report.md"))
message("Formal module-preservation objects combined and plotted.")
