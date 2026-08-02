# Environment validation report

Generated: 2026-08-02 (Asia/Shanghai)

## Scope

This validation was restricted to environment metadata, dependency availability, lockfile parsing, an isolated temporary `renv` restore attempt, `renv::status()`, and syntax-only parsing of R scripts. No scientific script body was executed and no scientific output was created or changed.

## Runtime and dependency record

| Item | Result |
|---|---|
| R runtime | R 4.5.2 (2025-10-31, UCRT) |
| Bioconductor | 3.22 |
| Temporary validation renv | 1.2.3, installed only in an isolated temporary validation library outside the repository |
| Recorded packages | 22 |
| Recorded packages installed in audit environment | 22/22 |
| Canonical recorded-version matches | 22/22 |
| R scripts parsed | 16/16 |
| R parse failures | 0 |

The displayed installed versions `1.1.38`, `3.1.3`, and `5.0.1` correspond to the canonical package versions `1.1-38`, `3.1-3`, and `5.0-1` recorded for lme4, lmerTest, and metafor.

## Lockfile inspection

`environment/renv.lock` is valid JSON and can be read by renv. It records R 4.5.2 but contains an empty `Packages` object and nonstandard audit-placeholder fields. It is therefore not a complete dependency lock.

| Test | Technical result | Release interpretation |
|---|---|---|
| JSON parsing | PASS | Structural parsing only |
| renv lockfile parsing | PASS | Structural parsing only |
| Recorded R version | 4.5.2 | Matches the locked analysis runtime |
| Locked package count | 0 | FAIL: major and transitive dependencies are absent |
| `renv::restore(prompt = FALSE)` in isolated temporary project | Returned success with "library is already synchronized" | FAIL for reproducibility: restore was a no-op because the lock contains zero packages |
| `renv::status()` | Reported that the temporary project does not appear to be using renv | FAIL |

## R syntax and static dependency checks

All 16 retained R scripts passed `parse()` without evaluation. Twenty non-base packages were detected through `library()`, `require()`, or namespace operators; every detected package is present in `environment/package_versions.csv`. No interactive input, file chooser, TODO/FIXME marker, or explicit hidden manual step was detected.

Static analysis cannot establish that every dataset-dependent object exists at runtime. Candidate external symbols identified by `codetools` were reviewed and consisted of base calls, ggplot/data-column names used through non-standard evaluation, or documented outputs from upstream numbered steps. No obvious standalone undefined variable was identified.

## Conclusion

**Environment validation: FAIL - NOT READY FOR RELEASE.**

The local audit environment contains the recorded packages, and every R script is syntactically valid. However, the checked-in lockfile is an empty audit placeholder. A complete lockfile must be generated in the intended R 4.5.2/Bioconductor 3.22 environment, restored in a genuinely clean project, and followed by a synchronized `renv::status()` result. Until then, the repository must not claim turnkey reproducibility or publish a formal `1.0.0` release.

Scientific lock: `raw_data_modified = FALSE`; module count = 24; Tier A/B/C/D = 0/1/21/2.
