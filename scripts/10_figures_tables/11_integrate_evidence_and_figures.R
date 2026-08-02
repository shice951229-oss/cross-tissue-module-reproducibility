########################################
# Script: 11_integrate_evidence_and_figures.R
# Purpose: Integrate locked result tables, assign locked evidence tiers, and generate main and supplementary figures.
# Input: Final model, negative-control, preservation, HPA, and preprocessing tables.
# Output: Evidence matrix, figure-source tables, and Figures 1-6.
# Software: R
# Version: 4.5.2
# Random seed: Not applicable (deterministic integration)
# Author: Study authors
########################################
#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 200)
suppressPackageStartupMessages({library(ggplot2); library(patchwork)})
root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
out <- file.path(root, "13_methodological_reanalysis_v9")
z <- function(x) as.numeric(scale(x))
readc <- function(name) read.csv(file.path(out, name), check.names = FALSE)
pretty_module <- function(x) {
  x <- sub("^hepatic_", "", x); x <- sub("^WGCNA_", "WGCNA ", x)
  x <- sub("_module$", "", x); gsub("_", " ", x)
}
save_plot <- function(p, stem, width_mm = 180, height_mm = 135) {
  ggsave(file.path(out, paste0(stem, ".pdf")), p, width = width_mm / 25.4, height = height_mm / 25.4,
         device = grDevices::cairo_pdf, family = "Arial")
  ggsave(file.path(out, paste0(stem, ".png")), p, width = width_mm / 25.4, height = height_mm / 25.4,
         dpi = 600, device = ragg::agg_png)
}
required <- c("matched_random_module_results_AD.csv", "matched_random_module_results_MDD.csv",
              "liver_module_preservation_results.csv", "module_liver_specificity_HPA_GTEx.csv")
if (any(!file.exists(file.path(out, required)))) stop("Required final analyses are missing: ", paste(required[!file.exists(file.path(out, required))], collapse = ", "))

mods <- readc("module_gene_lists_24_summary.csv")
names(mods)[names(mods) %in% c("gene_count", "n_genes")] <- "module_gene_n"
if (!"module" %in% names(mods)) names(mods)[1] <- "module"
mods$module_class <- ifelse(grepl("^WGCNA_", mods$module), "GSE135251 WGCNA", "curated liver-related")
mods$module_label <- pretty_module(mods$module)

# External liver disease/risk context, recomputed from audited TMM-logCPM ssGSEA scores.
s126 <- readRDS(file.path(out, "module_scores_main_old_mapping", "liver_GSE126848.rds"))
ph126 <- readRDS(file.path(root, "04_intermediate", "phenotype_tables", "liver_GSE126848", "pheno_harmonized.rds"))
ph126 <- ph126[match(colnames(s126), ph126$sample_id), ]
grp <- ifelse(tolower(ph126[["disease:ch1"]]) == "healthy", "healthy", "disease_risk")
ext_rows <- lapply(rownames(s126), function(m) {
  y1 <- as.numeric(s126[m, grp == "disease_risk"]); y0 <- as.numeric(s126[m, grp == "healthy"])
  n1 <- length(y1); n0 <- length(y0); sp <- sqrt(((n1 - 1) * var(y1) + (n0 - 1) * var(y0)) / (n1 + n0 - 2))
  d <- (mean(y1) - mean(y0)) / sp; se <- sqrt(1 / n1 + 1 / n0 + d^2 / (2 * (n1 + n0)))
  tt <- t.test(y1, y0)
  data.frame(module = m, liver_context_effect = d, liver_context_SE = se,
             liver_context_CI_low = d - 1.96 * se, liver_context_CI_high = d + 1.96 * se,
             liver_context_P = tt$p.value, liver_context_disease_risk_n = n1,
             liver_context_healthy_n = n0, liver_context_sample_n = n1 + n0,
             liver_context_donor_n = NA_integer_,
             liver_context_formula = "Welch two-sample t-test; Cohen d = disease/risk minus healthy")
})
ext <- merge(mods[, "module", drop = FALSE], do.call(rbind, ext_rows), by = "module", all.x = TRUE, sort = FALSE)
ext$liver_context_BH_FDR <- p.adjust(ext$liver_context_P, method = "BH", n = nrow(ext))
write.csv(ext, file.path(out, "external_liver_context_results_audited.csv"), row.names = FALSE, na = "")

ad <- readc("AD_dataset_specific_effects_strict_main.csv")
ad483 <- ad[ad$dataset == "GSE48350", ]; ad330 <- ad[ad$dataset == "GSE33000", ]
names(ad483)[names(ad483) %in% c("beta", "SE", "CI_low", "CI_high", "P", "BH_FDR", "sample_n", "donor_n", "formula")] <-
  paste0("GSE48350_", names(ad483)[names(ad483) %in% c("beta", "SE", "CI_low", "CI_high", "P", "BH_FDR", "sample_n", "donor_n", "formula")])
names(ad330)[names(ad330) %in% c("beta", "SE", "CI_low", "CI_high", "P", "BH_FDR", "sample_n", "donor_n", "formula")] <-
  paste0("GSE33000_", names(ad330)[names(ad330) %in% c("beta", "SE", "CI_low", "CI_high", "P", "BH_FDR", "sample_n", "donor_n", "formula")])
g483 <- readc("GSE48350_naive_vs_donoraware_full.csv")
reml <- readc("AD_two_dataset_REML_strict_main.csv")
hk <- readc("AD_two_dataset_Hartung_Knapp_strict_main.csv")
names(reml)[names(reml) != "module"] <- paste0("conventional_REML_", names(reml)[names(reml) != "module"])
names(hk)[names(hk) != "module"] <- paste0("Hartung_Knapp_", names(hk)[names(hk) != "module"])
mdd_blood <- readc("MDD_GSE98793_strict_main_results.csv")
names(mdd_blood)[names(mdd_blood) %in% c("effect", "SE", "CI_low", "CI_high", "P", "BH_FDR", "formula")] <-
  paste0("MDD_blood_", names(mdd_blood)[names(mdd_blood) %in% c("effect", "SE", "CI_low", "CI_high", "P", "BH_FDR", "formula")])
mdd_brain <- readc("GSE102556_overall_model_results.csv")
names(mdd_brain)[names(mdd_brain) != "module"] <- paste0("MDD_brain_", names(mdd_brain)[names(mdd_brain) != "module"])
null_ad <- readc("matched_random_module_results_AD.csv")
null_mdd <- readc("matched_random_module_results_MDD.csv")
hpa <- readc("module_liver_specificity_HPA_GTEx.csv")
pres <- readc("liver_module_preservation_results.csv")
if (!"permutations" %in% names(pres) || any(pres$permutations != 500L)) {
  stop("Formal liver preservation lock failed: every row must come from 500 permutations")
}

pres_wide <- reshape(pres[, c("module", "test_dataset", "common_gene_n", "coverage_fraction", "Zsummary", "medianRank", "preservation_class")],
                     idvar = "module", timevar = "test_dataset", direction = "wide")
names(pres_wide) <- gsub("\\.", "_", names(pres_wide))
pres$preservation_any_moderate <- as.integer(pres$common_gene_n >= 20 & pres$Zsummary >= 2)
pres$preservation_any_strong <- as.integer(pres$common_gene_n >= 20 & pres$Zsummary >= 10)
pres_summary <- aggregate(cbind(preservation_any_moderate, preservation_any_strong) ~ module, data = pres, FUN = max)
mdd_null_blood <- null_mdd[null_mdd$dataset == "GSE98793_blood", ]
mdd_null_brain <- null_mdd[null_mdd$dataset == "GSE102556_brain", ]
names(mdd_null_blood)[names(mdd_null_blood) != "module"] <- paste0("null_MDD_blood_", names(mdd_null_blood)[names(mdd_null_blood) != "module"])
names(mdd_null_brain)[names(mdd_null_brain) != "module"] <- paste0("null_MDD_brain_", names(mdd_null_brain)[names(mdd_null_brain) != "module"])

parts <- list(mods, ext, ad483, ad330, g483, reml, hk, mdd_blood, mdd_brain, null_ad, hpa, pres_wide, pres_summary,
              mdd_null_blood, mdd_null_brain)
for (i in seq_along(parts)) {
  hit <- names(parts[[i]]) == "mapping"
  names(parts[[i]])[hit] <- paste0("mapping_source_", i)
}
ev <- Reduce(function(a, b) merge(a, b, by = "module", all.x = TRUE, sort = FALSE), parts)
names(ev) <- make.unique(names(ev), sep = "_source")

# Evidence levels are algorithmic and deliberately conservative.
ev$liver_specificity_support <- !is.na(ev$BH_FDR_liver_enriched) & ev$BH_FDR_liver_enriched < .05
ev$external_liver_structure_support <- !is.na(ev$preservation_any_moderate) & ev$preservation_any_moderate == 1
ev$liver_source_support <- ev$liver_specificity_support | ev$external_liver_structure_support
ev$AD_direction_consistent <- sign(ev$GSE48350_beta) == sign(ev$GSE33000_beta)
ev$AD_both_dataset_FDR <- ev$GSE48350_BH_FDR < .05 & ev$GSE33000_BH_FDR < .05
ev$AD_any_dataset_FDR <- ev$GSE48350_BH_FDR < .05 | ev$GSE33000_BH_FDR < .05
ev$matched_random_support_any <- ev$empirical_BH_FDR_mean_abs_beta < .05 | ev$empirical_BH_FDR_cross_dataset_difference < .05 |
  ev$empirical_BH_FDR_abs_REML < .05
ev$final_evidence_tier <- "Tier D"
ev$final_evidence_tier[ev$AD_direction_consistent | ev$AD_any_dataset_FDR] <- "Tier C"
ev$final_evidence_tier[ev$liver_source_support & ev$AD_direction_consistent & ev$AD_any_dataset_FDR & ev$matched_random_support_any] <- "Tier B"
ev$final_evidence_tier[ev$liver_source_support & ev$AD_direction_consistent & ev$AD_both_dataset_FDR & ev$matched_random_support_any] <- "Tier A"
ev$allowed_conclusion <- ifelse(ev$final_evidence_tier == "Tier A",
  "Adjusted AD association in both datasets with independent liver support and matched-null excess; association only.",
  ifelse(ev$final_evidence_tier == "Tier B",
    "Limited directionally consistent AD projection with liver annotation and matched-null support; small-k uncertainty remains.",
    ifelse(ev$final_evidence_tier == "Tier C",
      "Dataset-specific or otherwise hypothesis-generating cross-tissue association.",
      "Dataset-dependent, directionally inconsistent, tissue-atypical, or technically vulnerable signal.")))
ev$forbidden_conclusion <- ifelse(ev$final_evidence_tier == "Tier A",
  "Hepatic origin, liver-to-brain transfer, mediation, mechanism, or causality.",
  "Independent replication, hepatic origin, liver-to-brain transfer, mediation, mechanism, or causality.")
ev$evidence_rule <- "A: liver support + direction consistency + both AD datasets FDR + matched-null FDR; B: liver support + direction consistency + >=1 AD dataset FDR + matched-null FDR; C: dataset-specific/limited; D: inconsistent or technically vulnerable."
ev <- ev[match(mods$module, ev$module), ]
write.csv(ev, file.path(out, "integrated_module_evidence_matrix.csv"), row.names = FALSE, na = "")

tier_counts <- table(factor(ev$final_evidence_tier, levels = c("Tier A", "Tier B", "Tier C", "Tier D")))
writeLines(c(
  "# Final evidence hierarchy", "", paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")), "",
  "## Prespecified conservative rules", "",
  "- Tier A: independent liver support, direction consistency, BH-FDR support in both AD datasets, and at least one matched-random empirical BH-FDR metric.",
  "- Tier B: independent liver support, direction consistency, BH-FDR support in at least one AD dataset, and matched-random empirical BH-FDR support; small-k instability or incomplete dataset support remains.",
  "- Tier C: single-dataset, non-random-exceeding, or otherwise hypothesis-generating evidence.",
  "- Tier D: directionally inconsistent, tissue-atypical, technically vulnerable, or unsupported evidence.", "",
  "## Counts", "",
  paste0("- Tier A: ", tier_counts["Tier A"]), paste0("- Tier B: ", tier_counts["Tier B"]),
  paste0("- Tier C: ", tier_counts["Tier C"]), paste0("- Tier D: ", tier_counts["Tier D"]), "",
  "## Universal inferential boundary", "",
  "No tier establishes hepatic tissue origin of a brain/blood score, circulation-mediated transfer, mediation, a liver-to-brain mechanism, or causality. Conventional two-dataset REML is secondary; Hartung-Knapp is retained as the small-k uncertainty analysis."
), file.path(out, "final_evidence_hierarchy.md"))

# Figure 1: design and inferential hierarchy.
nodes <- data.frame(x = c(1, 1, 1, 2.6, 4.2, 5.8), y = c(3.3, 2.2, 1.1, 2.2, 2.2, 2.2),
                    label = c("Liver reference\nGSE135251", "External liver\nGSE126848 / GSE167523", "HPA tissue\nconsensus",
                              "24 fixed modules\nssGSEA within dataset", "AD primary + MDD external\ndonor-aware / naive sensitivity",
                              "Matched random nulls\nEvidence tiers + limits"),
                    group = c("input", "input", "input", "score", "model", "inference"))
edges <- data.frame(x = c(1.45, 1.45, 1.45, 3.05, 4.65), y = c(3.3, 2.2, 1.1, 2.2, 2.2),
                    xend = c(2.15, 2.15, 2.15, 3.75, 5.35), yend = c(2.45, 2.2, 1.95, 2.2, 2.2))
p1 <- ggplot() + geom_segment(data = edges, aes(x, y, xend = xend, yend = yend), arrow = arrow(length = unit(2, "mm")), color = "#767676") +
  geom_label(data = nodes, aes(x, y, label = label, fill = group), size = 3, linewidth = .3, label.padding = unit(2.5, "mm"), lineheight = .95) +
  scale_fill_manual(values = c(input = "#EAF2F8", score = "#E8F5E9", model = "#FFF3E0", inference = "#FDEDEC")) +
  annotate("text", x = 3.4, y = 0.78, label = "Dataset-specific effects; no raw-score pooling", size = 3, color = "#4D4D4D") +
  annotate("text", x = 3.4, y = 0.48, label = "Association ≠ origin, mediation, or causality", size = 3, color = "#A63603") +
  coord_cartesian(xlim = c(.35, 6.45), ylim = c(.25, 3.75), clip = "off") + theme_void(base_family = "Arial") +
  theme(legend.position = "none", plot.title = element_text(face = "bold", size = 12)) +
  labs(title = "Audit-first cross-tissue module projection and evidence hierarchy")
save_plot(p1, "Figure1_study_design_evidence_hierarchy", 180, 100)

# Figure 2: preprocessing and donor overview.
reg <- readc("dataset_preprocessing_registry.csv")
reg$dataset <- factor(reg$GEO_accession, levels = rev(reg$GEO_accession))
donor_lookup <- c(GSE102556 = 48, GSE48350 = 84)
reg$donor_n <- unname(donor_lookup[reg$GEO_accession])
reg$sample_label <- ifelse(is.na(reg$donor_n), paste0("n=", reg$final_sample_n), paste0("n=", reg$final_sample_n, "; donors=", reg$donor_n))
p2a <- ggplot(reg, aes(final_sample_n, dataset, color = tissue, shape = technology)) +
  geom_point(size = 3) + geom_text(aes(label = sample_label), hjust = -.12, size = 2.5, color = "#333333") +
  scale_x_continuous(expand = expansion(mult = c(.03, .36))) + theme_classic(base_size = 7, base_family = "Arial") +
  theme(legend.position = "bottom") + labs(x = "Final analytical samples", y = NULL, color = "Tissue", shape = "Technology")
proc <- data.frame(dataset = reg$dataset, item = "Audited primary processing", value = reg$audited_main_processing)
p2b <- ggplot(proc, aes(1, dataset, label = value)) + geom_text(hjust = 0, size = 2.35, lineheight = .9) +
  coord_cartesian(xlim = c(1, 2.9), clip = "off") + theme_void(base_family = "Arial") +
  labs(title = "Dataset-specific decision")
p2 <- p2a + p2b + plot_layout(widths = c(1, 1.45)) + plot_annotation(title = "Preprocessing and independent-unit audit")
save_plot(p2, "Figure2_preprocessing_donor_overview", 180, 120)

# Figure 3: same-sample GSE48350 direct comparison.
g483$module_label <- pretty_module(g483$module)
g483$transition_label <- ifelse(g483$naive_BH_FDR < .05 & g483$donoraware_BH_FDR >= .05, "Naive FDR only", "Neither model FDR")
p3a <- ggplot(g483, aes(naive_beta, donoraware_beta, color = transition_label)) +
  geom_hline(yintercept = 0, color = "#D0D0D0") + geom_vline(xintercept = 0, color = "#D0D0D0") +
  geom_abline(slope = 1, linetype = "dashed", color = "#767676") + geom_point(size = 2.2) +
  scale_color_manual(values = c("Naive FDR only" = "#D24B40", "Neither model FDR" = "#767676")) +
  theme_classic(base_size = 7, base_family = "Arial") + theme(legend.position = "bottom") +
  labs(x = "Naive diagnosis coefficient (SD)", y = "Donor-aware coefficient (SD)", color = NULL,
       title = "A  Effect estimates")
g483$SE_ratio <- g483$donoraware_SE / g483$naive_SE
g483$module_label <- factor(g483$module_label, levels = g483$module_label[order(g483$SE_ratio)])
p3b <- ggplot(g483, aes(SE_ratio, module_label, color = transition_label)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#767676") + geom_point(size = 1.8) +
  scale_color_manual(values = c("Naive FDR only" = "#D24B40", "Neither model FDR" = "#767676")) +
  theme_classic(base_size = 6.5, base_family = "Arial") + theme(legend.position = "none") +
  labs(x = "Donor-aware SE / naive SE", y = NULL, title = "B  Uncertainty")
p3 <- p3a + p3b + plot_annotation(title = "GSE48350: identical samples and covariates, alternative correlation structures")
save_plot(p3, "Figure3_GSE48350_naive_vs_donoraware", 180, 125)

# Figure 4 and complete supplementary forest: dataset-specific AD effects.
ad$module_label <- pretty_module(ad$module)
ord <- unique(ad$module_label[order(ad$beta)])
ad$module_label <- factor(ad$module_label, levels = ord)
p4 <- ggplot(ad, aes(beta, module_label, color = dataset)) +
  geom_vline(xintercept = 0, color = "#BDBDBD") +
  geom_errorbar(aes(xmin = CI_low, xmax = CI_high), orientation = "y", width = .22, position = position_dodge(width = .55)) +
  geom_point(aes(shape = BH_FDR < .05), position = position_dodge(width = .55), size = 1.7) +
  scale_color_manual(values = c(GSE48350 = "#3182BD", GSE33000 = "#D24B40")) +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), labels = c("BH-FDR ≥0.05", "BH-FDR <0.05")) +
  theme_classic(base_size = 7, base_family = "Arial") + theme(legend.position = "bottom") +
  labs(x = "Adjusted AD diagnosis coefficient (within-dataset SD)", y = NULL, color = NULL, shape = NULL,
       title = "Dataset-specific AD effects and 95% confidence intervals")
save_plot(p4, "Figure4_AD_dataset_specific_effects", 180, 155)
save_plot(p4, "Supplementary_Figure_all_module_forest", 180, 155)

# Figure 5: matched random-module negative controls.
nad <- readRDS(file.path(out, "matched_random_module_null_AD.rds"))
write.csv(nad[, c("module", "replicate", "mean_abs_beta", "direction_consistent", "cross_dataset_difference",
                  "dual_nominal_same_direction", "reml_effect", "abs_beta_GSE48350", "abs_beta_GSE33000")],
          file.path(out, "figure_source_random_null_AD_long.csv"), row.names = FALSE, na = "")
nmdd <- readRDS(file.path(out, "matched_random_module_null_MDD.rds"))
write.csv(nmdd[, c("dataset", "module", "replicate", "effect", "P")],
          file.path(out, "figure_source_random_null_MDD_long.csv"), row.names = FALSE, na = "")
obs <- null_ad; obs$module_label <- pretty_module(obs$module)
ordn <- obs$module[order(obs$mean_abs_beta)]
nad$module_label <- pretty_module(nad$module); nad$module_label <- factor(nad$module_label, levels = pretty_module(ordn))
obs$module_label <- factor(obs$module_label, levels = pretty_module(ordn))
p5a <- ggplot(nad, aes(mean_abs_beta, module_label)) + geom_violin(fill = "#D9D9D9", color = NA, scale = "width") +
  geom_point(data = obs, aes(mean_abs_beta, module_label, color = empirical_BH_FDR_mean_abs_beta < .05), inherit.aes = FALSE, size = 1.5) +
  scale_color_manual(values = c(`FALSE` = "#767676", `TRUE` = "#D24B40"), labels = c("No", "Yes")) +
  theme_classic(base_size = 6.2, base_family = "Arial") + theme(legend.position = "bottom") +
  labs(x = "Mean absolute AD effect", y = NULL, color = "Empirical BH-FDR <0.05", title = "A  Effect magnitude")
p5b <- ggplot(nad, aes(cross_dataset_difference, module_label)) + geom_violin(fill = "#D9D9D9", color = NA, scale = "width") +
  geom_point(data = obs, aes(cross_dataset_difference, module_label, color = empirical_BH_FDR_cross_dataset_difference < .05), inherit.aes = FALSE, size = 1.5) +
  scale_color_manual(values = c(`FALSE` = "#767676", `TRUE` = "#3182BD"), labels = c("No", "Yes")) +
  theme_classic(base_size = 6.2, base_family = "Arial") + theme(legend.position = "bottom") +
  labs(x = "Absolute between-AD-dataset effect difference", y = NULL, color = "Empirical BH-FDR <0.05", title = "B  Reproducibility")
p5 <- p5a + p5b + plot_annotation(title = "Observed liver-defined modules against B=2000 matched random gene sets")
save_plot(p5, "Figure5_matched_random_negative_control", 180, 155)

# Figure 6: preservation plus final evidence tier.
ev$module_label <- pretty_module(ev$module)
tier_levels <- c("Tier D", "Tier C", "Tier B", "Tier A")
ev$final_evidence_tier <- factor(ev$final_evidence_tier, levels = tier_levels)
ev$module_label <- factor(ev$module_label, levels = ev$module_label[order(ev$final_evidence_tier)])
p6a <- ggplot(ev, aes(final_evidence_tier, module_label, fill = final_evidence_tier)) + geom_tile(color = "white") +
  scale_fill_manual(values = c("Tier A" = "#1B7837", "Tier B" = "#5AAE61", "Tier C" = "#FDB863", "Tier D" = "#BDBDBD")) +
  theme_minimal(base_size = 6.5, base_family = "Arial") + theme(legend.position = "none", panel.grid = element_blank()) +
  labs(x = NULL, y = NULL, title = "A  Final conservative tier")
pres$module_label <- pretty_module(pres$module)
p6b <- ggplot(pres, aes(Zsummary, module_label, color = preservation_class, size = common_gene_n)) +
  geom_vline(xintercept = 2, linetype = "dashed", color = "#767676") + geom_vline(xintercept = 10, linetype = "dotted", color = "#4D4D4D") +
  geom_point(alpha = .85) + facet_wrap(~test_dataset) +
  scale_color_manual(values = c("no clear evidence" = "#BDBDBD", moderate = "#3182BD", strong = "#D24B40", "low common-gene count" = "#7B3294")) +
  scale_size_continuous(range = c(1, 4)) +
  guides(size = guide_legend(order = 1, nrow = 1), color = guide_legend(order = 2, nrow = 2, byrow = TRUE)) +
  theme_classic(base_size = 6.5, base_family = "Arial") +
  theme(legend.position = "bottom", legend.box = "vertical", legend.text = element_text(size = 5.5)) +
  labs(x = "Preservation Zsummary", y = NULL, color = NULL, size = "Common genes", title = "B  External liver network preservation")
p6 <- p6a + p6b + plot_layout(widths = c(.55, 1.45)) + plot_annotation(title = "Evidence hierarchy and source-attribute checks")
save_plot(p6, "Figure6_final_evidence_and_liver_preservation", 180, 155)

# Supplementary donor-structure overview for GSE102556.
cw <- readc("GSE102556_sample_donor_crosswalk.csv")
per_donor <- aggregate(brain_region ~ donor_id + diagnosis, cw, function(x) length(unique(x)))
pS1a <- ggplot(per_donor, aes(brain_region, fill = diagnosis)) + geom_histogram(binwidth = 1, boundary = .5, position = "dodge") +
  scale_x_continuous(breaks = 1:6) + theme_classic(base_size = 7, base_family = "Arial") +
  labs(x = "Regions contributed per donor", y = "Independent donors", fill = NULL, title = "A  Repeated-region structure")
byreg <- aggregate(GSM ~ brain_region + diagnosis, cw, length)
pS1b <- ggplot(byreg, aes(GSM, brain_region, fill = diagnosis)) + geom_col(position = "dodge") +
  theme_classic(base_size = 7, base_family = "Arial") + labs(x = "Samples", y = NULL, fill = NULL, title = "B  Samples by region")
save_plot(pS1a + pS1b, "Supplementary_Figure_GSE102556_donor_structure", 180, 100)

message("Integrated evidence matrix and final R figures completed.")
