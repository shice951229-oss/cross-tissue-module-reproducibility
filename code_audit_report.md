# Code audit report

## Scope and non-execution boundary

This pre-release audit examined repository structure, file references, script headers, static package declarations, R syntax, path safety, credential patterns, obvious external-object dependencies, hidden interactive steps, and file-size thresholds. It did not execute GSVA, WGCNA, mixed models, CR2, REML, Hartung-Knapp, matched random modules, module preservation, HPA analysis, figure-generation scripts, or any other scientific workflow.

Scientific locks remain unchanged: `raw_data_modified = FALSE`; module count = 24; Tier A/B/C/D = 0/1/21/2.

## Final static findings

| Check | Status | Evidence and implication |
|---|---|---|
| README repository references | PASS | All repository files and directories named as execution or documentation resources exist. References to `01_raw_data`, `04_intermediate`, and other local staging inputs are explicitly described as excluded external inputs. |
| Standard script headers | PASS | All 20 scripts (16 R, 2 Python, 2 PowerShell) include Script, Purpose, Input, Output, Software, Version, Random seed, and Author fields. |
| R syntax | PASS | All 16 R scripts passed `parse()` without executing their bodies. |
| Direct R dependencies | PASS | All 20 statically detected non-base packages are recorded in `environment/package_versions.csv`. `R.utils` and `yaml` are retained as locked-environment or repository-infrastructure dependencies and are labeled accordingly. |
| Installed package check | PASS WITH LIMITATION | All 22 recorded packages are installed in the audit environment at the recorded canonical versions; hyphen/dot display differences for lme4, lmerTest, and metafor are equivalent package-version renderings. This does not replace a lockfile restore test. |
| Obvious undefined variables | PASS WITH LIMITATION | A `codetools::findGlobals()` heuristic produced only base calls, data-column names used through non-standard evaluation, and documented upstream objects. No obvious misspelled standalone symbol was identified. Semantic execution was intentionally not attempted because it would require locked scientific inputs. |
| Hidden manual or interactive steps | PASS | No `interactive()`, `readline()`, `file.choose()`, `choose.files()`, `menu()`, TODO, FIXME, or explicit copy-by-hand instruction occurs in the retained R scripts. Command-line arguments and environment variables are explicit. |
| Personal/local paths | PASS | No Windows drive-letter path, `/Users/` path, local username, or temporary personal directory is embedded in public scripts. Project-relative staging paths remain documented in `config.yaml`. |
| Credentials and identifiers | PASS | No password, token, API key, private key, email address, ORCID, or credential-like string was detected in public text or OOXML content. |
| Third-party raw data | PASS | No RDS/RData, CEL, RAW, compressed GEO archive, HDF5, FASTQ, BAM, or other raw expression object is present. Hash-only raw-data audit records do not contain the underlying data. |
| File-size thresholds | PASS | No repository file exceeds 25 MiB; consequently none exceeds 100 MiB. |
| Required directories | PASS | `figures`, `metadata`, `results`, `docs`, `environment`, and all numbered script directories are present. |
| Deterministic figure export | PASS | The TIFF conversion utility is explicit; six final figures are 600-dpi LZW TIFF files. |
| Complete environment lock | PASS | `environment/renv.lock` records 188 packages for R 4.5.2/Bioconductor 3.22. A non-no-op sandboxed clean-project restore completed with zero missing packages and zero version mismatches; `renv::status()` reported a consistent state. |
| Remote repository | PASS WITH LIMITATION | Local `origin` is configured to `https://github.com/shice951229-oss/cross-tissue-module-reproducibility`, and the owner reports that the private repository was created and the initial commit was pushed. A read-only automated remote check was interrupted by a network connection reset, so public accessibility is not independently asserted. |

## Direct dependency record

The locked analysis used R 4.5.2 and Bioconductor 3.22. Directly referenced non-base packages are AnnotationDbi, Biobase, BiocParallel, clubSandwich, digest, edgeR, GEOquery, ggplot2, GSVA, hgu133plus2.db, limma, lme4, lmerTest, matrixStats, metafor, org.Hs.eg.db, patchwork, ragg, scales, and WGCNA. Recorded supporting dependencies include R.utils and yaml. Python is limited to deterministic format and consistency utilities.

## Explicit execution boundary

The retained scientific scripts require GEO/HPA inputs and outputs from prior numbered steps. Their required objects are documented in script headers, `config.yaml`, `code_manifest.csv`, and the reproducibility protocol. Because those third-party matrices are intentionally excluded, static checks cannot prove end-to-end semantic execution. This is a documented input boundary, not a hidden manual step.

## Release assessment

The code package is syntactically valid, public-safe, and reproducibly environment-locked. Release status remains **NOT READY** because no immutable `v1.0.0` tag/release or Zenodo DOI exists and final public-release authorization is still required. No scientific result was rerun or changed during environment validation.
