# Cross-tissue transcriptomic module reproducibility analysis

[![DOI](https://zenodo.org/badge/1320331422.svg)](https://doi.org/10.5281/zenodo.21762838)

## Overview

This repository contains the analysis code and supporting materials for evaluating reproducibility and interpretational limits of cross-tissue transcriptomic module projections. It accompanies the manuscript:

> Reproducibility of liver-relevant transcriptomic module projections across Alzheimer disease and major depressive disorder

The analysis evaluates 8 predefined liver-relevant inflammatory-metabolic programs and 16 WGCNA modules derived in GSE135251 after projection into Alzheimer disease (AD) and major depressive disorder (MDD) transcriptomic datasets. The emphasis is reproducibility, dataset dependence, donor-aware uncertainty, matched random-module controls, external liver-network preservation, and conservative evidence grading. Cross-tissue association does not establish tissue origin, inter-organ transfer, mediation, or causality.

Scientific results are locked. This release does not redefine modules, alter preprocessing or models, rerun analyses to improve significance, remove negative findings, or change the Tier A/B/C/D hierarchy (0/1/21/2). `raw_data_modified = FALSE`; `modules = 24`; `Tier A/B/C/D = 0/1/21/2`.

## Release status

**PUBLIC RELEASE AVAILABLE**

The repository is publicly accessible at `https://github.com/shice951229-oss/cross-tissue-module-reproducibility`. The `v1.0.0` GitHub release archives validated commit `3eb6474`. Zenodo archived that release under version DOI `10.5281/zenodo.21762839`; the concept DOI for all versions is `10.5281/zenodo.21762838`. Author and institutional authorization for public release under the MIT License was confirmed before publication.

## Data availability

Raw transcriptomic datasets were obtained from the Gene Expression Omnibus (GEO) under the reported accession numbers:

- GSE135251, GSE167523, and GSE126848 (liver datasets)
- GSE98793 and GSE102556 (MDD blood and brain)
- GSE144136 (published marker resource)
- GSE48350, GSE5281, and GSE33000 (AD brain datasets)

GEO raw files, raw/processed expression matrices, large third-party source-data mirrors, and other third-party data are intentionally excluded. They must be downloaded from their original repositories. The local download staging directory is ignored by version control.

The Human Protein Atlas tissue-consensus file used for the optional liver-expression attribute analysis is also not redistributed; obtain the applicable release directly from the Human Protein Atlas and record its release information locally.

## Repository structure

- `environment/`: R/Python versions, recorded package versions, seed registry, and environment-lock status.
- `scripts/01_data_download/`: public-data download, input discovery, inventory, and raw-file integrity checks.
- `scripts/02_preprocessing/`: dataset-specific preprocessing and redundant-normalization audit.
- `scripts/03_probe_mapping/`: ambiguous-probe mapping sensitivity.
- `scripts/04_module_scoring/`: fixed-module GSVA/ssGSEA scoring.
- `scripts/05_donor_models/`: donor recovery and donor-aware versus naive model comparisons.
- `scripts/06_random_controls/`: matched random-module negative controls for AD and MDD.
- `scripts/07_meta_analysis/`: strict dataset effects, CR2 sensitivity, REML, and Hartung-Knapp synthesis.
- `scripts/08_module_preservation/`: WGCNA module-preservation analysis and postprocessing.
- `scripts/09_HPA_analysis/`: descriptive HPA liver-expression attribute analysis.
- `scripts/10_figures_tables/`: evidence integration, consistency audits, and final figure export.
- `metadata/`: dataset, sample, donor, module-definition, and probe-mapping registries.
- `results/`: locked tabular results, model outputs, figure source data, and evidence hierarchy.
- `figures/`: final six 600-dpi LZW TIFF figures.
- `docs/`: audit reports, reproducibility protocol, workflow diagrams, and public/private guidance.

## Requirements

The locked analysis was run with R 4.5.2. Directly used R packages and versions are recorded in `environment/package_versions.csv`; principal packages include GSVA, BiocParallel, edgeR, limma, WGCNA, lme4, lmerTest, metafor, clubSandwich, AnnotationDbi, hgu133plus2.db, org.Hs.eg.db, ggplot2, patchwork, matrixStats, ragg, scales, digest, yaml, GEOquery, and Biobase. `variancePartition` is not listed as a required package because the locked scripts do not use it.

Python 3.12.13 and Pillow 12.2.0 are used only for deterministic figure-format export; the consistency audit uses the Python standard library.

`environment/renv.lock` is a complete 188-package lock for R 4.5.2 and Bioconductor 3.22. It was restored in a newly initialized isolated project with zero missing packages and zero version mismatches; `renv::status()` reported a consistent state. The validation executed no scientific script body.

## Configuration and local data staging

`config.yaml` centralizes accessions, seeds, parameters, expected output locations, and the excluded local analysis staging root. The archived locked scripts retain the project-relative structure used in the analysis (`01_raw_data`, `04_intermediate`, and the locked result directory). Set `CTMR_ANALYSIS_ROOT` and `CTMR_RESULTS_DIR`, or stage those directories as described in `docs/reproducibility_protocol.pdf`, before executing them.

No repository script should contain a personal username or a Windows absolute path. Run the public-safety scan described in `code_audit_report.md` after any local edits.

## Reproduction workflow

The numbered folders represent the dependency order:

1. download or stage third-party inputs and verify source-file integrity;
2. apply the registered dataset-specific preprocessing and missingness audit;
3. create strict unambiguous probe mappings and sensitivity outputs;
4. score the fixed modules with GSVA/ssGSEA within each dataset;
5. recover verified donor structures and run donor-aware/naive comparisons;
6. run matched random-module negative controls with seed 20260802;
7. estimate dataset-specific effects, CR2 sensitivity, REML, and Hartung-Knapp results;
8. run 500-permutation external liver module-preservation analyses;
9. run the optional descriptive HPA liver-expression analysis;
10. integrate the locked evidence matrix, generate figures/tables, and run consistency checks.

The repository is not claimed to be one-command end-to-end reproducible because raw and processed third-party matrices are intentionally excluded. The dependency environment itself has passed clean-room restoration. Expected inputs and outputs for each retained script are documented in its header and in `code_manifest.csv`.

## Reproducing figures

- Fig1: `scripts/10_figures_tables/17_update_figures_v10.R` (text-only locked-table refinement)
- Fig2: `scripts/10_figures_tables/11_integrate_evidence_and_figures.R`
- Fig3: `scripts/10_figures_tables/11_integrate_evidence_and_figures.R`
- Fig4: `scripts/10_figures_tables/11_integrate_evidence_and_figures.R`
- Fig5: `scripts/10_figures_tables/11_integrate_evidence_and_figures.R`
- Fig6: `scripts/10_figures_tables/17_update_figures_v10.R` (text-only locked-table refinement)
- TIFF export: `scripts/10_figures_tables/10_export_tiff_figures.py`

The final plotted values are provided in `results/figure_source_data/figure_source_data.xlsx` and accompanying long-form CSV files.

## Reproducibility notes

- GSE102556 contains 263 samples from 48 verified donors; the donor is the independent statistical unit for the overall brain model.
- GSE48350 contains 253 samples from 84 verified donors; naive and donor-aware models use identical samples and covariates.
- GSE33000 lacks verified donor identifiers in the available project metadata. GSE5281 is descriptive and excluded from independent synthesis.
- Dataset-level scores are never merged across datasets. Diagnostic coefficients are standardized within dataset.
- BH families are documented explicitly, including the 24 overall-module families and the 144 region-module family.
- Conventional REML and Hartung-Knapp results remain separate. With two AD datasets, Hartung-Knapp is the required small-k uncertainty analysis.
- Matched random modules use seed 20260802 and B=2000 per observed module.
- Module preservation uses seed 20260802 and 500 permutations per external liver dataset.
- Evidence tiers remain Tier A/B/C/D = 0/1/21/2. No tier licenses tissue-origin, transfer, mediation, or causal claims.

## Citation

Use `CITATION.cff`. Cite the archived `v1.0.0` release with DOI `10.5281/zenodo.21762839`. The stable concept DOI for all repository versions is `10.5281/zenodo.21762838`.

## License

Author-generated code in this repository is released under the MIT License. This license does not relicense or redistribute GEO, Human Protein Atlas, or other third-party data.
