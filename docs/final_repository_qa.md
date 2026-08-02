# Final repository QA

Updated: 2026-08-03 (Asia/Shanghai)

## Overall release status

**PASS WITH LIMITATIONS**

The repository is publicly accessible at `https://github.com/shice951229-oss/cross-tissue-module-reproducibility`. The immutable `v1.0.0` release remains pinned to validated commit `3eb6474`, and Zenodo has issued version DOI `10.5281/zenodo.21762839` and concept DOI `10.5281/zenodo.21762838`. The local package passed the static code, public-safety, file-integrity, metadata, and clean-project environment-restoration checks described below. The only remaining engineering action is to push this post-release metadata synchronization commit to the default branch; the existing release tag must not be moved.

## Verification summary

| Gate | Result | Evidence |
|---|---:|---|
| Scientific lock | PASS | `raw_data_modified = FALSE`; 1,115/1,115 raw files unchanged; 24 modules; Tier A/B/C/D = 0/1/21/2 |
| README and standard files | PASS | README, LICENSE, CITATION.cff, `.zenodo.json`, configuration, scripts, metadata, results, figure/table source data and environment records are present |
| Citation and Zenodo metadata | PASS | CFF/YAML and JSON parse; author order matches the manuscript; title, version 1.0.0, authors, MIT license, public repository URL and verified Zenodo version DOI are synchronized; unverified ORCID and email fields remain omitted |
| Script headers and syntax | PASS | 20/20 scripts have the required standard header; 16/16 R scripts pass `parse()` |
| Dependency declarations | PASS | 20 directly detected non-base R dependencies are declared; 22/22 recorded package versions match the installed validation environment |
| Hidden manual steps | PASS | No interactive prompts, `readline`, file chooser, TODO/FIXME marker or undisclosed manual intervention was detected |
| Public/private separation | PASS | No GEO RAW archive, raw/processed expression matrix, RDS/RData object, download cache or third-party raw-data mirror is included |
| Privacy and credential scan | PASS | No local absolute user path, Windows username, password, token, API key, private key, author email or ORCID pattern was detected in public text or OOXML |
| Large-file gate | PASS | No file exceeds 25 MiB or 100 MiB |
| Spreadsheet integrity | PASS | 11/11 XLSX files reopen; expected sheets/ranges are present; zero formula-error hits |
| TIFF integrity | PASS | 6/6 figures are valid TIFF, 600 dpi and LZW-compressed |
| DOCX/PDF integrity | PASS WITH LIMITATION | Package documents open structurally and current DOCX text/OOXML checks found no corruption or mojibake. LibreOffice/soffice was unavailable for a fresh visual rerender of the updated DOCX files, so no new clipping/layout claim is made |
| Environment lock | PASS | 188-package lock for R 4.5.2/Bioconductor 3.22; non-no-op sandboxed clean-project restore; zero missing packages; zero version mismatches; clean `renv::status()` |
| GitHub remote | PASS | Public repository verified at `https://github.com/shice951229-oss/cross-tissue-module-reproducibility` |
| Immutable release and Zenodo DOI | PASS | `v1.0.0` points to validated commit `3eb6474`; Zenodo version DOI is `10.5281/zenodo.21762839` and concept DOI is `10.5281/zenodo.21762838` |

## Environment validation

- R version: 4.5.2.
- Bioconductor version: 3.22.
- Temporary validation library: `renv` 1.2.3.
- Lock-file parse: PASS.
- Package records in lock file: 188.
- Restore process: PASS; non-no-op clean-project restoration completed in 395.94 seconds.
- Missing packages: 0; version mismatches: 0.
- `renv::status()`: PASS; `No issues found -- the project is in a consistent state.`
- Scientific analysis scripts executed: none.

## Final file and manifest reconciliation

- Final repository files, including the manifest itself: 110.
- Manifest rows: 109.
- Expected exclusions: one, the self-referential `docs/repository_file_manifest_sha256.csv`.
- SHA-256 or byte-size verification failures: 0.
- `repository_inventory.xlsx` rows: 110, including the manifest entry and documented release exclusions.

The manifest is authoritative for the frozen local package. No manifest-covered file may be edited after hash generation; any later metadata, dependency-lock, URL, DOI or release change requires regeneration of both the inventory and the manifest.

## Release decision and remaining human action

1. Push the post-release metadata synchronization commit to the public repository and confirm it is visible on the default branch.
2. Do not move or recreate the existing `v1.0.0` tag; it must continue to identify validated commit `3eb6474`.
3. In the Zenodo record editor, synchronize the descriptive Notes with `.zenodo.json` if the web record still displays the earlier pre-release wording. This does not change the DOI or archived files.

## Scientific lock statement

This engineering audit did not rerun or alter GSVA/ssGSEA, WGCNA, mixed-effects models, CR2, REML, Hartung-Knapp, matched random modules, module preservation, HPA analysis, figure source data or manuscript conclusions. Sample and donor counts, module definitions, effect estimates, P values, FDR values, figures and evidence tiers remain unchanged.
