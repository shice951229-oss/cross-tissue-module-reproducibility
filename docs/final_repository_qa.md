# Final repository QA

Generated: 2026-08-02 (Asia/Shanghai)

## Overall release status

**NOT READY**

The local package passed the static code, public-safety, file-integrity and metadata checks described below. Release remains blocked because `environment/renv.lock` contains zero package records, the target GitHub remote has not been independently verified, no immutable release exists, and no Zenodo DOI has been issued.

## Verification summary

| Gate | Result | Evidence |
|---|---:|---|
| Scientific lock | PASS | `raw_data_modified = FALSE`; 1,115/1,115 raw files unchanged; 24 modules; Tier A/B/C/D = 0/1/21/2 |
| README and standard files | PASS | README, LICENSE, CITATION.cff, `.zenodo.json`, configuration, scripts, metadata, results, figure/table source data and environment records are present |
| Citation and Zenodo metadata | PASS | CFF/YAML and JSON parse; author order matches the manuscript; title, version 1.0.0, authors and MIT license are synchronized; no DOI, ORCID or email is asserted |
| Script headers and syntax | PASS | 20/20 scripts have the required standard header; 16/16 R scripts pass `parse()` |
| Dependency declarations | PASS | 20 directly detected non-base R dependencies are declared; 22/22 recorded package versions match the installed validation environment |
| Hidden manual steps | PASS | No interactive prompts, `readline`, file chooser, TODO/FIXME marker or undisclosed manual intervention was detected |
| Public/private separation | PASS | No GEO RAW archive, raw/processed expression matrix, RDS/RData object, download cache or third-party raw-data mirror is included |
| Privacy and credential scan | PASS | No local absolute user path, Windows username, password, token, API key, private key, author email or ORCID pattern was detected in public text or OOXML |
| Large-file gate | PASS | No file exceeds 25 MiB or 100 MiB |
| Spreadsheet integrity | PASS | 11/11 XLSX files reopen; expected sheets/ranges are present; zero formula-error hits |
| TIFF integrity | PASS | 6/6 figures are valid TIFF, 600 dpi and LZW-compressed |
| DOCX/PDF integrity | PASS | Package documents open structurally; the updated audit DOCX files and retained explanatory PDFs were rendered and visually inspected without clipping, overlap, broken tables or missing glyphs |
| Environment lock | **FAIL** | `renv.lock` is valid JSON but has zero package records; `renv::restore(prompt = FALSE)` was a no-op and `renv::status()` reported that the temporary project did not appear to use renv |
| GitHub remote | **FAIL** | Owner and intended URL are recorded, but repository creation/access and empty-remote status have not been independently verified |
| Immutable release and Zenodo DOI | **FAIL** | No `v1.0.0` release/tag and no Zenodo DOI exist |

## Environment validation

- R version: 4.5.2.
- Bioconductor version: 3.22.
- Temporary validation library: `renv` 1.2.3.
- Lock-file parse: PASS.
- Package records in lock file: 0.
- Restore process exit: success, but restoration result: no-op; this is not accepted as a reproducible restore.
- `renv::status()`: FAIL for release purposes because the temporary project was not recognized as an initialized renv project.
- Scientific analysis scripts executed: none.

## Final file and manifest reconciliation

- Final repository files, including the manifest itself: 110.
- Manifest rows: 109.
- Expected exclusions: one, the self-referential `docs/repository_file_manifest_sha256.csv`.
- SHA-256 or byte-size verification failures: 0.
- `repository_inventory.xlsx` rows: 110, including the manifest entry and documented release exclusions.

The manifest is authoritative for the frozen local package. No manifest-covered file may be edited after hash generation; any later metadata, dependency-lock, URL, DOI or release change requires regeneration of both the inventory and the manifest.

## Release decision and required human actions

1. Generate a real, complete `renv.lock` from the documented R 4.5.2/Bioconductor 3.22 environment and demonstrate a non-no-op clean-room restore followed by a clean `renv::status()`.
2. Create or independently verify the GitHub repository at the owner-supplied target, confirm the remote is empty before the first push, and push the frozen package. The intended URL is `https://github.com/shice951229-oss/cross-tissue-module-reproducibility`; this report does not claim that it is currently public or verified.
3. Obtain final author/institutional authorization for public release and use of the MIT license.
4. Only after all release gates pass, create the immutable `v1.0.0` tag/release, archive that release in Zenodo, obtain the real DOI, synchronize all citation and availability metadata, and rerun the final safety scan and hash freeze.

## Scientific lock statement

This engineering audit did not rerun or alter GSVA/ssGSEA, WGCNA, mixed-effects models, CR2, REML, Hartung-Knapp, matched random modules, module preservation, HPA analysis, figure source data or manuscript conclusions. Sample and donor counts, module definitions, effect estimates, P values, FDR values, figures and evidence tiers remain unchanged.
