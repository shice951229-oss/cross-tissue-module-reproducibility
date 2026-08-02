########################################
# Script: 08_liver_module_preservation.R
# Purpose: Run formal WGCNA module-preservation tests in two external liver datasets.
# Input: Reference/test liver matrices and 16 fixed WGCNA modules.
# Output: Preservation statistics, saved modulePreservation objects, and figures.
# Software: R
# Version: 4.5.2
# Random seed: 20260802; 500 permutations per external liver dataset
# Author: Study authors
########################################
#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 180)
suppressPackageStartupMessages({
  library(WGCNA)
  library(ggplot2)
  library(patchwork)
})

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
out <- file.path(root, "13_methodological_reanalysis_v9")
seed <- 20260802L
n_perm <- as.integer(Sys.getenv("PRESERVATION_PERM", "500"))
if (!is.finite(n_perm) || n_perm < 1L) stop("PRESERVATION_PERM must be a positive integer")
set.seed(seed)
allowWGCNAThreads(nThreads = 4)

modules_all <- readRDS(file.path(root, "04_intermediate", "gene_sets", "liver_modules_combined.rds"))
modules <- modules_all[grepl("^WGCNA_", names(modules_all))]
if (length(modules) != 16L) stop("Expected exactly 16 WGCNA modules; found ", length(modules))
modules <- lapply(modules, function(x) unique(toupper(trimws(as.character(x)))))

ref0 <- readRDS(file.path(out, "preprocessed_main_old_mapping", "liver_GSE135251.rds"))
tests <- list(
  GSE126848 = readRDS(file.path(out, "preprocessed_main_old_mapping", "liver_GSE126848.rds")),
  GSE167523 = readRDS(file.path(out, "preprocessed_main_old_mapping", "liver_GSE167523.rds"))
)
requested_test <- Sys.getenv("PRESERVATION_TEST", "")
if (nzchar(requested_test)) {
  if (!requested_test %in% names(tests)) stop("Unknown PRESERVATION_TEST: ", requested_test)
  tests <- tests[requested_test]
}
file_prefix <- if (nzchar(requested_test)) paste0("partial_", requested_test, "_") else ""

# Reproduce the reference gene eligibility rule used by the original network construction.
ref_var <- apply(ref0, 1, var, na.rm = TRUE)
var_threshold <- quantile(ref_var[ref_var > 0], 0.10, na.rm = TRUE)
ref_eligible <- names(ref_var)[is.finite(ref_var) & ref_var > var_threshold]
module_union <- unique(unlist(modules, use.names = FALSE))

extract_pair <- function(container, test_name) {
  # WGCNA returns reference-by-test nested lists; names vary slightly by release.
  ref_node <- container[[1L]]
  if (!is.list(ref_node)) stop("Unexpected modulePreservation result structure")
  hit <- grep(test_name, names(ref_node), fixed = TRUE)
  if (length(hit)) return(ref_node[[hit[1L]]])
  if (length(ref_node) >= 2L) return(ref_node[[2L]])
  stop("Could not locate test-network preservation results for ", test_name)
}

run_one <- function(test_name, test_expr) {
  common <- Reduce(intersect, list(ref_eligible, rownames(test_expr), module_union))
  color <- rep("grey", length(common)); names(color) <- common
  for (m in names(modules)) {
    color[intersect(common, modules[[m]])] <- sub("^WGCNA_", "", m)
  }
  keep <- color != "grey"
  common <- common[keep]
  color <- color[keep]
  if (anyDuplicated(common)) stop("Duplicated common genes")
  if (any(table(color) < 20L)) warning("At least one module has fewer than 20 common genes in ", test_name)

  multi_expr <- list(
    GSE135251 = list(data = t(ref0[common, , drop = FALSE])),
    test = list(data = t(test_expr[common, , drop = FALSE]))
  )
  names(multi_expr)[2L] <- test_name
  color_list <- list(GSE135251 = unname(color))

  message("Running modulePreservation: GSE135251 -> ", test_name,
          "; genes=", length(common), "; permutations=", n_perm)
  mp <- modulePreservation(
    multiData = multi_expr,
    multiColor = color_list,
    referenceNetworks = 1,
    nPermutations = n_perm,
    networkType = "unsigned",
    randomSeed = seed,
    quickCor = 0,
    maxGoldModuleSize = 1000,
    maxModuleSize = 10000,
    savePermutedStatistics = FALSE,
    parallelCalculation = TRUE,
    verbose = 3
  )
  saveRDS(mp, file.path(out, paste0("liver_module_preservation_", test_name, "_full.rds")))

  ztab <- as.data.frame(extract_pair(mp$preservation$Z, test_name))
  # In WGCNA 1.73, medianRank is an observed-statistic column rather than a top-level list.
  mtab <- as.data.frame(extract_pair(mp$preservation$observed, test_name))
  ztab$color <- rownames(ztab); mtab$color <- rownames(mtab)
  z_col <- grep("^Zsummary", names(ztab), value = TRUE)[1L]
  rank_col <- grep("^medianRank", names(mtab), value = TRUE)[1L]
  if (is.na(z_col) || is.na(rank_col)) stop("Expected Zsummary/medianRank columns were not returned")

  counts <- table(color)
  rows <- data.frame(
    module = paste0("WGCNA_", names(counts)),
    test_dataset = test_name,
    reference_dataset = "GSE135251",
    common_gene_n = as.integer(counts),
    reference_original_gene_n = vapply(paste0("WGCNA_", names(counts)), function(m) length(modules[[m]]), integer(1)),
    stringsAsFactors = FALSE
  )
  rows$coverage_fraction <- rows$common_gene_n / rows$reference_original_gene_n
  rows$Zsummary <- ztab[[z_col]][match(sub("^WGCNA_", "", rows$module), ztab$color)]
  rows$medianRank <- mtab[[rank_col]][match(sub("^WGCNA_", "", rows$module), mtab$color)]
  rows$common_gene_adequacy <- ifelse(rows$common_gene_n >= 20L, "adequate for interpretation", "too few; descriptive only")
  rows$preservation_class <- as.character(cut(rows$Zsummary, breaks = c(-Inf, 2, 10, Inf),
                                              labels = c("no clear evidence", "moderate", "strong"), right = FALSE))
  rows$preservation_class[rows$common_gene_n < 20L] <- "low common-gene count"
  rows$network_type <- "unsigned"
  rows$reference_gene_filter <- paste0("variance > 10th percentile of positive variance (", signif(var_threshold, 7), ")")
  rows$permutations <- n_perm
  rows$seed <- seed
  rows
}

res <- do.call(rbind, Map(run_one, names(tests), tests))
res <- res[order(res$test_dataset, res$medianRank), ]
write.csv(res, file.path(out, paste0(file_prefix, "liver_module_preservation_results.csv")), row.names = FALSE, na = "")

size_cor <- do.call(rbind, lapply(split(res, res$test_dataset), function(d) {
  ct <- cor.test(log10(d$common_gene_n), d$Zsummary, method = "spearman", exact = FALSE)
  data.frame(test_dataset = d$test_dataset[1], spearman_rho_log_size_vs_Zsummary = unname(ct$estimate), P = ct$p.value)
}))
write.csv(size_cor, file.path(out, paste0(file_prefix, "liver_module_preservation_size_effect.csv")), row.names = FALSE)

res$module_label <- sub("^WGCNA_", "", res$module)
p_z <- ggplot(res, aes(Zsummary, reorder(module_label, Zsummary))) +
  geom_vline(xintercept = 2, linetype = "dashed", color = "#767676") +
  geom_vline(xintercept = 10, linetype = "dotted", color = "#4D4D4D") +
  geom_point(aes(size = common_gene_n, color = preservation_class), alpha = 0.9) +
  facet_wrap(~test_dataset, scales = "free_y") +
  scale_color_manual(values = c("no clear evidence" = "#BDBDBD", "moderate" = "#3182BD", "strong" = "#D24B40",
                                "low common-gene count" = "#7B3294"), drop = FALSE) +
  scale_size_continuous(range = c(1.5, 5)) +
  theme_classic(base_size = 7, base_family = "Arial") +
  theme(legend.position = "bottom") +
  labs(x = "WGCNA preservation Zsummary", y = NULL, color = "Evidence", size = "Common genes",
       title = "External preservation of GSE135251 WGCNA modules")
p_rank <- ggplot(res, aes(medianRank, reorder(module_label, -medianRank))) +
  geom_point(aes(size = common_gene_n), color = "#4D4D4D", alpha = 0.85) +
  facet_wrap(~test_dataset, scales = "free_y") +
  scale_size_continuous(range = c(1.5, 5)) +
  theme_classic(base_size = 7, base_family = "Arial") +
  theme(legend.position = "bottom") +
  labs(x = "Preservation medianRank (lower is better)", y = NULL, size = "Common genes",
       title = "Size-aware relative preservation ranking")
p <- p_z / p_rank + plot_annotation(tag_levels = "A")
ggsave(file.path(out, paste0(file_prefix, "Figure_liver_module_preservation.pdf")), p, width = 180 / 25.4, height = 200 / 25.4,
       device = grDevices::cairo_pdf, family = "Arial")
ggsave(file.path(out, paste0(file_prefix, "Figure_liver_module_preservation.png")), p, width = 180 / 25.4, height = 200 / 25.4,
       dpi = 600, device = ragg::agg_png)

summary_line <- function(ds) {
  ok <- res$test_dataset == ds & res$common_gene_n >= 20L
  s <- sum(ok & res$Zsummary >= 10, na.rm = TRUE)
  m <- sum(ok & res$Zsummary >= 2 & res$Zsummary < 10, na.rm = TRUE)
  n <- sum(ok & res$Zsummary < 2, na.rm = TRUE)
  low <- sum(res$test_dataset == ds & res$common_gene_n < 20L)
  paste0("- ", ds, ": ", s, " strong, ", m, " moderate, ", n,
         " with no clear preservation evidence, and ", low, " modules marked descriptive because fewer than 20 genes were shared.")
}
writeLines(c(
  "# External liver module-preservation report", "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")), "",
  "## Design", "",
  "- Reference network: GSE135251; test networks: GSE126848 and GSE167523.",
  "- Only the 16 pre-defined GSE135251 WGCNA modules were evaluated.",
  paste0("- Formal WGCNA `modulePreservation` used an unsigned network, ", n_perm,
         " permutations, and fixed seed ", seed, "."),
  "- Each reference-test comparison used genes present in the original reference eligibility set, the named WGCNA modules, and the test expression matrix.",
  "- Zsummary was interpreted jointly with medianRank and common module size; color labels were not treated as biological functions.", "",
  "## Results", "", vapply(names(tests), summary_line, character(1)), "",
  "## Inferential boundary", "",
  "Preservation supports recurrence of within-module co-expression structure in external liver datasets. It does not establish that projected brain or blood signals originate in liver, circulate between tissues, mediate disease, or are causal."
), file.path(out, paste0(file_prefix, "liver_module_preservation_report.md")))

message("Completed liver module preservation with ", n_perm, " permutations per test.")
