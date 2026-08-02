# Preanalysis status report

## Reproducibility status at lock

- `raw_data_modified = FALSE`.
- All 1,115 files under `01_raw_data` were readable and received SHA-256 baselines.
- One non-raw supplementary DOCX was open in another process and could not be hashed during the first pass; this does not block any statistical analysis.
- The 24 module definitions are complete: eight curated modules and 16 WGCNA modules, 35-5,967 genes per module.
- R 4.5.2 and the required GSVA, limma, edgeR, lme4/lmerTest, metafor, clubSandwich, geepack and WGCNA packages are available.

## Analyses that can be rerun completely

- GSE102556 donor crosswalk, overall mixed model, identical-sample naive comparison, prespecified sex interaction and six-region descriptive analyses.
- GSE48350 identical-sample naive, donor-aware and donor-clustered robust comparisons for all 24 modules.
- Dataset-specific preprocessing and missingness audit for the eight bulk-expression datasets.
- Additional quantile-normalization sensitivity for all four microarray datasets.
- GPL570 and GPL4372 strict official-symbol mapping sensitivity.
- Matched random-module AD negative control using GSE135251, GSE48350 and GSE33000.
- External liver WGCNA module preservation in GSE126848 and GSE167523.
- HPA tissue-expression specificity because the HPA consensus archive is present locally.
- Conventional REML and Hartung-Knapp two-dataset synthesis.

## Verified sample information

- GSE102556: donor IDs are explicit in GEO sample titles and agree with both the project `subject_id`/sample name and the downloaded FPKM column names. No ordering- or diagnosis-based reconstruction is required.
- GSE48350: 253 samples map to 84 donors using explicit GEO individual/sample-title fields.
- GSE33000: the AD/control analytical set is 467 samples (310 AD, 157 controls); explicit donor identifiers remain unavailable.
- GSE5281: within-region samples are usable descriptively, but a complete cross-region donor map is not verified.
- GSE167523: all 98 samples belong to the disease spectrum under the project classifier; it cannot supply a binary healthy comparison but remains usable for unsupervised preservation.

## Information not yet confirmable from project files

- GSE33000 independent donor count, RIN, PMI and technical batch.
- GSE5281 cross-region donor identity.
- A full GSE144136 count matrix in the project; the saved expression object is a 0-feature placeholder, so the dataset is restricted to its published marker-resource role.
- The original free-text GPL570 multi-symbol annotation was not retained; strict sensitivity uses current official Bioconductor probe-to-symbol mappings and reports that provenance.

## Preanalysis code/manuscript discrepancies

- The old RNA-seq branch used only `max(expression)>1000` to choose `log2(x+1)` and otherwise applied quantile normalization. This is not source-aware and is replaced by dataset-specific decisions.
- The old microarray workflow reapplied quantile normalization to matrices whose GEO metadata report RMA, GC-RMA/PLIER or MAS5 processing.
- The old annotation code split multi-symbol annotations and kept the first symbol arbitrarily.
- v8 described GSE102556 without an explicit donor-aware overall model despite repeated donors.
- v8 stated that repeated-donor accounting materially changed AD inference without presenting a same-dataset GSE48350 naive comparison.
- v8 contains mediator terminology despite no formal mediation analysis.
- v8 has no matched-random gene-set negative control or formal external liver module-preservation analysis.

## Preanalysis conclusion

The core requested reanalysis is feasible from existing local inputs. Missing donor information for GSE33000 and GSE5281 limits only donor-level interpretation for those datasets and is not repaired by imputation or invented identifiers.

