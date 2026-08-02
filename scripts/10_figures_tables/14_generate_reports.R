########################################
# Script: 14_generate_reports.R
# Purpose: Generate the locked methodological audit reports from final result tables.
# Input: Final CSV outputs and integrity records.
# Output: Markdown audit, evidence, preservation, and numeric consistency reports.
# Software: R
# Version: 4.5.2
# Random seed: Not applicable (deterministic)
# Author: Study authors
########################################
#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 200)
root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
out <- file.path(root, "13_methodological_reanalysis_v9")
rc <- function(x) read.csv(file.path(out, x), check.names = FALSE)
fmt <- function(x, d = 3) formatC(x, digits = d, format = "f")

g102 <- rc("GSE102556_naive_vs_donoraware.csv")
g483 <- rc("GSE48350_naive_vs_donoraware_full.csv")
ad <- rc("AD_dataset_specific_effects_strict_main.csv")
reml <- rc("AD_two_dataset_REML_strict_main.csv")
hk <- rc("AD_two_dataset_Hartung_Knapp_strict_main.csv")
adnull <- rc("matched_random_module_results_AD.csv")
adglob <- rc("global_negative_control_tests.csv")
mddnull <- rc("matched_random_module_results_MDD.csv")
mddglob <- rc("global_negative_control_tests_MDD.csv")
pres <- rc("liver_module_preservation_results.csv")
hpa <- rc("module_liver_specificity_HPA_GTEx.csv")
probe <- rc("probe_annotation_audit.csv")
psens <- rc("probe_mapping_sensitivity_summary.csv")
ev <- rc("integrated_module_evidence_matrix.csv")

g102_beta_cor <- cor(g102$naive_beta, g102$donoraware_beta)
g102_med_ratio <- median(g102$naive_SE / g102$donoraware_SE)
g483_beta_cor <- cor(g483$naive_beta, g483$donoraware_beta)
g483_med_ratio <- median(g483$donoraware_SE / g483$naive_SE)
ad330_sig <- sum(ad$dataset == "GSE33000" & ad$BH_FDR < .05)
ad483_sig <- sum(ad$dataset == "GSE48350" & ad$BH_FDR < .05)
reml_sig <- reml$module[reml$BH_FDR < .05]
hk_sig <- hk$module[hk$BH_FDR < .05]
ad_emp <- adnull$module[adnull$empirical_BH_FDR_mean_abs_beta < .05 | adnull$empirical_BH_FDR_abs_REML < .05]
mdd_blood_emp <- mddnull$module[mddnull$dataset == "GSE98793_blood" & mddnull$empirical_BH_FDR_abs_effect < .05]
mdd_brain_emp <- mddnull$module[mddnull$dataset == "GSE102556_brain" & mddnull$empirical_BH_FDR_abs_effect < .05]
pres_summary <- do.call(rbind, lapply(split(pres, pres$test_dataset), function(d) data.frame(
  dataset = d$test_dataset[1], strong = sum(d$common_gene_n >= 20 & d$Zsummary >= 10),
  moderate = sum(d$common_gene_n >= 20 & d$Zsummary >= 2 & d$Zsummary < 10),
  none = sum(d$common_gene_n >= 20 & d$Zsummary < 2), low_common = sum(d$common_gene_n < 20)
)))
tier_counts <- table(factor(ev$final_evidence_tier, levels = c("Tier A", "Tier B", "Tier C", "Tier D")))
hpa_sig <- hpa$module[hpa$BH_FDR_liver_enriched < .05]

audit <- c(
  "# Analysis audit report v9", "", paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")), "",
  "## Overall status", "", "**PASS WITH LIMITATIONS.** The requested primary audits, donor-aware comparisons, preprocessing/mapping sensitivities, B=2000 matched-random controls, formal 500-permutation liver module preservation, HPA specificity annotation, evidence integration and manuscript reframing were completed from local inputs. Limitations are explicit below.", "",
  "## Phase 0: inventory and raw-data lock", "",
  "- The project inventory, SHA-256 manifest, dependency map and preanalysis status report were created before reanalysis.",
  "- All 24 module definitions were complete (8 curated, 16 WGCNA).",
  "- Raw inputs were read-only. Final SHA-256 verification is recorded in `raw_data_integrity_check.md`; `raw_data_modified = FALSE` is required for closure.", "",
  "## GSE102556 donor audit", "",
  "- 263 brain samples mapped to 48 independently identified donors (26 MDD, 22 control) using explicit GEO title prefixes corroborated by project sample IDs and FPKM columns; no order- or diagnosis-based reconstruction was used.",
  "- Primary formula: `score_z ~ diagnosis + brain_region + age_z + sex + RIN_z + PMI_z + (1 | donor_id)`.",
  paste0("- Naive versus donor-aware beta correlation: ", fmt(g102_beta_cor, 4), "; median naive SE / donor-aware SE: ", fmt(g102_med_ratio, 4), "."),
  paste0("- BH-FDR modules: naive ", sum(g102$naive_BH_FDR < .05), "/24; donor-aware ", sum(g102$donoraware_BH_FDR < .05), "/24. Region analysis: 0/144 under global BH; sex interaction: 0/24."),
  paste0("- Singular/near-zero donor variance was flagged in ", sum(as.logical(g102$singular_fit)), " of 24 mixed models; models were retained and flagged rather than deleted."), "",
  "## GSE48350 same-sample model comparison", "",
  "- 253 samples from 84 verified donors; identical samples, score standardization and fixed covariates were used in both models.",
  "- Naive formula: `score_z ~ diagnosis + region + age_z + sex`; donor-aware formula adds `(1 | donor_id)`.",
  paste0("- Naive BH-FDR: ", sum(g483$naive_BH_FDR < .05), "/24; donor-aware BH-FDR: ", sum(g483$donoraware_BH_FDR < .05), "/24; beta correlation: ", fmt(g483_beta_cor, 4), "; median donor-aware SE / naive SE: ", fmt(g483_med_ratio, 4), "."),
  "- Accounting for within-donor repetition widened uncertainty and reduced naive significance in this dataset. It did not by itself explain all differences between GSE48350 and GSE33000.", "",
  "## Dataset-specific preprocessing", "",
  "- Count matrices (GSE135251, GSE126848, GSE167523) used edgeR TMM logCPM; GSE102556 used log2(FPKM+1); processed microarrays retained their reported processed values.",
  "- Four processed microarray datasets had received an additional quantile normalization in legacy code. Because ssGSEA was rank-based, module-score correlations, diagnosis-beta correlations and median SE ratios were all 1.000 and no FDR status changed. The redundant step was nevertheless removed from the final method.",
  "- Missingness was quantified for every matrix; only GSE5281 required nonzero median imputation, affecting 2.59e-5 of cells.", "",
  "## Probe mapping", "",
  paste0("- GPL570 contained ", probe$strict_ambiguous_probe_n[probe$platform == "GPL570"][1], " ambiguous probes and ", probe$strict_unmapped_probe_n[probe$platform == "GPL570"][1], " unmapped probes; strict IQR mapping retained ", probe$strict_IQR_final_gene_n[probe$platform == "GPL570"][1], " genes."),
  paste0("- GPL4372 retained ", probe$strict_IQR_final_gene_n[probe$platform == "GPL4372"][1], " genes through official Entrez-to-symbol mapping; no multi-symbol ambiguity was observed in that mapping."),
  paste0("- Minimum mean module-score correlation under strict-IQR mapping was ", fmt(min(psens$mean_score_correlation[psens$comparison_mapping == "strict_max_iqr"]), 3), "; GSE33000 had one FDR transition. The prespecified stability rule therefore promoted strict unambiguous/highest-IQR mapping to the final main analysis."), "",
  "## AD synthesis and matched-random control", "",
  paste0("- Strict main results: GSE33000 ", ad330_sig, "/24 BH-FDR; donor-aware GSE48350 ", ad483_sig, "/24."),
  paste0("- Conventional REML BH-FDR: ", length(reml_sig), "/24 (", paste(reml_sig, collapse = ", "), "); Hartung-Knapp BH-FDR: ", length(hk_sig), "/24."),
  paste0("- B=2000 empirical module-level effect/REML BH-FDR was met by: ", paste(ad_emp, collapse = ", "), ". No module passed empirical BH-FDR for smaller cross-dataset difference or dual same-direction support."),
  paste0("- Global mean absolute AD effect: empirical P=", signif(adglob$empirical_P[adglob$metric == "mean_abs_effect"], 4), ", BH-FDR=", signif(adglob$BH_FDR[adglob$metric == "mean_abs_effect"], 5), ". Direction consistency and reproducibility did not exceed random expectations."), "",
  "## MDD matched-random controls", "",
  paste0("- GSE98793 blood: global mean absolute effect empirical BH-FDR=", signif(mddglob$BH_FDR[mddglob$analysis == "GSE98793_blood" & mddglob$metric == "mean_abs_effect"], 4), "; module-level empirical BH-FDR: ", ifelse(length(mdd_blood_emp), paste(mdd_blood_emp, collapse = ", "), "none"), "."),
  paste0("- GSE102556 brain: global mean absolute effect empirical BH-FDR=", signif(mddglob$BH_FDR[mddglob$analysis == "GSE102556_brain" & mddglob$metric == "mean_abs_effect"], 4), "; module-level empirical BH-FDR: ", ifelse(length(mdd_brain_emp), paste(mdd_brain_emp, collapse = ", "), "none"), "."),
  "- Blood and brain were analyzed and reported separately; no cross-tissue MDD pooled estimate was calculated.", "",
  "## Liver source-attribute checks", "",
  unlist(lapply(seq_len(nrow(pres_summary)), function(i) paste0("- ", pres_summary$dataset[i], ": ", pres_summary$strong[i], " strong, ", pres_summary$moderate[i], " moderate, ", pres_summary$none[i], " no clear evidence, ", pres_summary$low_common[i], " low-common-gene descriptive modules."))),
  paste0("- HPA strict liver-enrichment BH-FDR: ", length(hpa_sig), "/24 (", paste(hpa_sig, collapse = ", "), "). This is descriptive specificity annotation, not tissue-origin evidence."), "",
  "## Final evidence tiers", "",
  paste0("- Tier A: ", tier_counts["Tier A"], "; Tier B: ", tier_counts["Tier B"], "; Tier C: ", tier_counts["Tier C"], "; Tier D: ", tier_counts["Tier D"], "."),
  "- No tier permits claims of hepatic origin, circulation-mediated transfer, mediation, a liver-to-brain mechanism or causality.", "",
  "## Remaining limitations", "",
  "1. GSE33000 donor IDs, RIN, PMI and technical batch were unavailable; its model remains sample-level age/sex adjusted.",
  "2. GSE5281 cross-region donor identities were unresolved; it remains descriptive and excluded from independent pooling.",
  "3. GSE144136 project expression is a 0-feature placeholder; only its published marker resource is retained.",
  "4. Some small WGCNA modules shared fewer than 20 genes with external liver tests and are not interpreted from Zsummary thresholds.",
  "5. Two AD datasets remain insufficient for stable random-effects inference; Hartung-Knapp results are prominently retained.",
  "6. Postmortem tissue composition, RNA quality, platform and cohort heterogeneity remain plausible explanations for cross-dataset differences."
)
writeLines(audit, file.path(out, "analysis_audit_report_v9.md"))

revision <- c(
  "# Revision log: v8 to v9", "", paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")), "",
  "## Reframing", "",
  "- Replaced the mechanism-seeking framing with an audit of reproducibility, dataset dependence, liver source attributes, negative controls and inferential boundaries.",
  "- New title: *Cross-tissue projection of liver-derived transcriptomic modules in Alzheimer disease and major depressive disorder: a reproducibility and negative-control analysis*.",
  "- Removed affirmative liver-brain-axis, hepatic-origin, liver-mediated, mediation and independent-replication claims.", "",
  "## Methods", "",
  "- Added project/code audit, source-aware preprocessing registry, GSE102556 donor recovery, same-sample GSE48350 naive comparison, strict probe mapping, B=2000 matched-random controls, formal 500-permutation module preservation, HPA specificity and evidence-tier rules.",
  "- Replaced the `max > 1000` classifier and automatic second quantile normalization with dataset-specific decisions documented from source metadata, platform and distributions.",
  "- Made every model formula, statistical unit and BH family explicit, including the 144-test global family for region analyses.", "",
  "## Locked numerical changes", "",
  paste0("- GSE33000 BH-FDR modules: 18/24 in v8 to ", ad330_sig, "/24 under final strict mapping."),
  paste0("- GSE48350 same-sample comparison added: naive ", sum(g483$naive_BH_FDR < .05), "/24 versus donor-aware ", sum(g483$donoraware_BH_FDR < .05), "/24."),
  paste0("- Conventional REML: two modules in v8 to ", length(reml_sig), "/24 (", paste(reml_sig, collapse = ", "), ") under strict mapping."),
  paste0("- Hartung-Knapp remains ", length(hk_sig), "/24 BH-FDR and is no longer visually or narratively subordinated to conventional REML."),
  "- GSE102556 now uses 263 samples from 48 verified donors in the primary mixed model; naive and donor-aware analyses both yielded 0/24 BH-FDR.",
  paste0("- Matched-random AD effect-magnitude/REML empirical BH-FDR modules: ", paste(ad_emp, collapse = ", "), "; reproducibility metrics: 0/24."),
  paste0("- Matched-random MDD blood empirical BH-FDR modules: ", paste(mdd_blood_emp, collapse = ", "), "; MDD brain: 0/24."), "",
  "## Interpretation and terminology", "",
  "- Replaced mediator tier/prioritization terminology with follow-up candidate tier and cross-context candidate genes.",
  "- The 32 maximum-score candidates are explicitly described as tied without within-tier ranking and as inheriting module-level evidence.",
  "- WGCNA green is discussed as translation/ribosome/mitochondrial bioenergetics with alternative explanations including cell composition, neuronal loss, RNA quality and general metabolic state.",
  "- WGCNA greenyellow is treated as keratinization/epidermal differentiation and a tissue-atypical technical/compositional warning.", "",
  "## Figures and tables", "",
  "- Replaced the previous figure sequence with six R-generated main figures covering design, preprocessing/donor structure, GSE48350 direct comparison, dataset-specific AD effects, matched-random controls and final evidence/preservation.",
  "- Added supplementary donor-structure, region, preprocessing, probe-mapping and full-forest figures; all are provided as PDF and 600-dpi PNG.",
  "- Added source workbooks for preprocessing decisions, integrated evidence and every figure.", "",
  "## Preserved administrative content", "",
  "- Author list, affiliations, corresponding-author email, ethics, consent, author contributions, funding, competing interests and acknowledgements were retained without factual change."
)
writeLines(revision, file.path(out, "revision_log_v8_to_v9.md"))

message("Audit report and revision log generated.")
