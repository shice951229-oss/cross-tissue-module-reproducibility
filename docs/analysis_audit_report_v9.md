# Analysis audit report v9

Generated: 2026-08-02 09:35:14 +0800

## Overall status

**PASS WITH LIMITATIONS.** The requested primary audits, donor-aware comparisons, preprocessing/mapping sensitivities, B=2000 matched-random controls, formal 500-permutation liver module preservation, HPA specificity annotation, evidence integration and manuscript reframing were completed from local inputs. Limitations are explicit below.

## Phase 0: inventory and raw-data lock

- The project inventory, SHA-256 manifest, dependency map and preanalysis status report were created before reanalysis.
- All 24 module definitions were complete (8 curated, 16 WGCNA).
- Raw inputs were read-only. Final SHA-256 verification is recorded in `raw_data_integrity_check.md`; `raw_data_modified = FALSE` is required for closure.

## GSE102556 donor audit

- 263 brain samples mapped to 48 independently identified donors (26 MDD, 22 control) using explicit GEO title prefixes corroborated by project sample IDs and FPKM columns; no order- or diagnosis-based reconstruction was used.
- Primary formula: `score_z ~ diagnosis + brain_region + age_z + sex + RIN_z + PMI_z + (1 | donor_id)`.
- Naive versus donor-aware beta correlation: 0.9949; median naive SE / donor-aware SE: 0.6988.
- BH-FDR modules: naive 0/24; donor-aware 0/24. Region analysis: 0/144 under global BH; sex interaction: 0/24.
- Singular/near-zero donor variance was flagged in 3 of 24 mixed models; models were retained and flagged rather than deleted.

## GSE48350 same-sample model comparison

- 253 samples from 84 verified donors; identical samples, score standardization and fixed covariates were used in both models.
- Naive formula: `score_z ~ diagnosis + region + age_z + sex`; donor-aware formula adds `(1 | donor_id)`.
- Naive BH-FDR: 4/24; donor-aware BH-FDR: 0/24; beta correlation: 0.9918; median donor-aware SE / naive SE: 1.3521.
- Accounting for within-donor repetition widened uncertainty and reduced naive significance in this dataset. It did not by itself explain all differences between GSE48350 and GSE33000.

## Dataset-specific preprocessing

- Count matrices (GSE135251, GSE126848, GSE167523) used edgeR TMM logCPM; GSE102556 used log2(FPKM+1); processed microarrays retained their reported processed values.
- Four processed microarray datasets had received an additional quantile normalization in legacy code. Because ssGSEA was rank-based, module-score correlations, diagnosis-beta correlations and median SE ratios were all 1.000 and no FDR status changed. The redundant step was nevertheless removed from the final method.
- Missingness was quantified for every matrix; only GSE5281 required nonzero median imputation, affecting 2.59e-5 of cells.

## Probe mapping

- GPL570 contained 1527 ambiguous probes and 10028 unmapped probes; strict IQR mapping retained 20823 genes.
- GPL4372 retained 19866 genes through official Entrez-to-symbol mapping; no multi-symbol ambiguity was observed in that mapping.
- Minimum mean module-score correlation under strict-IQR mapping was 0.929; GSE33000 had one FDR transition. The prespecified stability rule therefore promoted strict unambiguous/highest-IQR mapping to the final main analysis.

## AD synthesis and matched-random control

- Strict main results: GSE33000 20/24 BH-FDR; donor-aware GSE48350 0/24.
- Conventional REML BH-FDR: 1/24 (WGCNA_green); Hartung-Knapp BH-FDR: 0/24.
- B=2000 empirical module-level effect/REML BH-FDR was met by: hepatic_mitochondrial_dysfunction_module, WGCNA_green. No module passed empirical BH-FDR for smaller cross-dataset difference or dual same-direction support.
- Global mean absolute AD effect: empirical P=0.01249, BH-FDR=0.049975. Direction consistency and reproducibility did not exceed random expectations.

## MDD matched-random controls

- GSE98793 blood: global mean absolute effect empirical BH-FDR=0.0004998; module-level empirical BH-FDR: hepatic_complement_coagulation_module.
- GSE102556 brain: global mean absolute effect empirical BH-FDR=0.956; module-level empirical BH-FDR: none.
- Blood and brain were analyzed and reported separately; no cross-tissue MDD pooled estimate was calculated.

## Liver source-attribute checks

- GSE126848: 9 strong, 1 moderate, 1 no clear evidence, 5 low-common-gene descriptive modules.
- GSE167523: 7 strong, 2 moderate, 2 no clear evidence, 5 low-common-gene descriptive modules.
- HPA strict liver-enrichment BH-FDR: 5/24 (hepatic_complement_coagulation_module, hepatic_lipid_dysregulation_module, hepatic_acute_phase_response_module, hepatic_secretome_hepatokine_module, WGCNA_brown). This is descriptive specificity annotation, not tissue-origin evidence.

## Final evidence tiers

- Tier A: 0; Tier B: 1; Tier C: 21; Tier D: 2.
- No tier permits claims of hepatic origin, circulation-mediated transfer, mediation, a liver-to-brain mechanism or causality.

## Remaining limitations

1. GSE33000 donor IDs, RIN, PMI and technical batch were unavailable; its model remains sample-level age/sex adjusted.
2. GSE5281 cross-region donor identities were unresolved; it remains descriptive and excluded from independent pooling.
3. GSE144136 project expression is a 0-feature placeholder; only its published marker resource is retained.
4. Some small WGCNA modules shared fewer than 20 genes with external liver tests and are not interpreted from Zsummary thresholds.
5. Two AD datasets remain insufficient for stable random-effects inference; Hartung-Knapp results are prominently retained.
6. Postmortem tissue composition, RNA quality, platform and cohort heterogeneity remain plausible explanations for cross-dataset differences.
