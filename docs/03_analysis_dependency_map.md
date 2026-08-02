# Analysis dependency map

## Locked sources

| Analysis object | Locked input | Role |
|---|---|---|
| 24 module definitions | `04_intermediate/gene_sets/liver_modules_combined.rds` | Eight curated programs plus 16 GSE135251 WGCNA modules |
| GSE135251 | `04_intermediate/expression_matrices/liver_GSE135251/expression_matrix_raw.rds` | Liver reference counts and WGCNA-module source |
| GSE126848 | `04_intermediate/expression_matrices/liver_GSE126848/expression_matrix_raw.rds` | External liver test network / context |
| GSE167523 | `04_intermediate/expression_matrices/liver_GSE167523/expression_matrix_raw.rds` | External liver test network; no binary disease-control effect |
| GSE98793 | `04_intermediate/expression_matrices/MDD_GSE98793/expression_matrix_raw.rds` plus saved GPL570 mapping | MDD blood projection and microarray mapping sensitivity |
| GSE102556 | downloaded human FPKM table plus GEO phenotype RDS | MDD brain projection with repeated donors |
| GSE144136 | published supplementary marker table and project metadata | Marker-resource analysis only; 0-feature expression placeholder is not analyzed |
| GSE48350 | probe matrix, phenotype RDS and `08_donor_aware_reanalysis/02_GSE48350_donor_region_map.csv` | Primary donor-aware AD analysis |
| GSE5281 | probe matrix and phenotype RDS | Region-stratified descriptive AD analysis; excluded from independent pooling |
| GSE33000 | Rosetta-processed probe matrix, GPL4372 annotation and phenotype RDS | Single-region independent AD dataset |
| Target manuscript | `12_main_manuscript_final_closure/01_conventional_journal_manuscript_v8_revised.docx` | Retained style and content reference; never overwritten |

## Current matrix versions used before the v9 audit

All reported v8 module scores were generated dataset by dataset from `04_intermediate/normalized_data/<dataset>/expression_matrix_normalized.rds`; the raw module scores were never pooled across datasets. The v9 audit reconstructs the input scale from the corresponding `04_intermediate/expression_matrices` file and writes all alternatives only under `13_methodological_reanalysis_v9`.

| Dataset | v8 matrix | Audited v9 main candidate |
|---|---|---|
| GSE135251 | generic `log2(count+1)` output | TMM-normalized logCPM |
| GSE126848 | generic `log2(count+1)` output | TMM-normalized logCPM |
| GSE167523 | generic `log2(count+1)` output | TMM-normalized logCPM |
| GSE102556 | `log2(FPKM+1)` | `log2(FPKM+1)` retained |
| GSE98793 | additional quantile normalization of reported RMA values | reported RMA processed values; additional quantile as sensitivity |
| GSE48350 | additional quantile normalization of reported GC-RMA/PLIER values | reported processed values; additional quantile as sensitivity |
| GSE5281 | additional quantile normalization of MAS5 target-scaled values | reported processed values; additional quantile as sensitivity |
| GSE33000 | additional quantile normalization of Rosetta values | reported Rosetta processed values; additional quantile as sensitivity |
| GSE144136 | no full expression projection | published marker table only |

## Execution graph

1. `00_phase0_inventory.ps1` creates the file and SHA-256 baseline.
2. `01_discover_inputs.R` inventories expression objects, phenotype fields and all 24 gene lists.
3. `02_donor_models.R` verifies GSE102556 identities and performs same-sample naive versus donor-aware comparisons for GSE102556 and GSE48350.
4. `04_preprocessing_audit.R` applies source-justified preprocessing, reruns ssGSEA and quantifies additional-quantile sensitivity.
5. `05_probe_mapping_sensitivity.R` compares legacy first-symbol mapping with strict official one-to-one mapping and prespecified duplicate-probe collapse rules.
6. Matched random-module testing uses GSE135251 expression strata and the common GSE135251/GSE48350/GSE33000 gene universe, then fits the same primary AD models.
7. WGCNA module preservation uses GSE135251 as reference and GSE126848/GSE167523 as test networks.
8. The evidence matrix integrates dataset-specific effects, conventional REML, Hartung-Knapp, negative controls, liver preservation and MDD projections.
9. Only after the numeric lock, the v8 DOCX is copied and edited into the v9 manuscript, then rendered page by page.

## Statistical families

- GSE102556 overall diagnosis effects: 24-test BH family.
- GSE102556 sex interaction sensitivity: separate 24-test BH family.
- GSE102556 region results: primary descriptive BH family across 144 tests; within-region 24-test BH values are sensitivity only.
- GSE48350 naive, donor-aware and cluster-robust results: separate 24-test BH families.
- GSE33000 dataset-specific effects: 24-test BH family.
- Conventional REML and Hartung-Knapp: separate 24-test BH families.
- Matched-random empirical tests: 24 real-module empirical P values per stated metric, BH-adjusted within metric.

