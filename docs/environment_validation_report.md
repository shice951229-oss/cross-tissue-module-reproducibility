# Environment validation report

Generated: 2026-08-03 (Asia/Shanghai)

## Scope

Validation was restricted to dependency-lock construction and parsing, an isolated clean-project `renv::restore(prompt = FALSE)`, `renv::status()`, package/version reconciliation, and syntax-only parsing of R scripts. No scientific script body was executed. GSVA/ssGSEA, WGCNA, donor models, CR2, REML, Hartung-Knapp, matched random modules, module preservation, HPA analysis, figure generation, and result-table generation were not run.

## Locked environment

| Item | Result |
|---|---|
| R runtime | R 4.5.2 (UCRT) |
| Bioconductor | 3.22 |
| Validation `renv` | 1.2.3 |
| Declared direct/supporting packages | 22 |
| Packages recorded in `environment/renv.lock` | 188 |
| Lock-file size | 418,588 bytes |
| Lock-file SHA-256 | `0fe3a7ab5c14117750b9e43f25e71729d9e72df4b0db436692ff1b3f8aac2a49` |
| Local-path records in lock file | 0 |

The 188 records contain the declared dependencies and their transitive closure. R-recommended packages supplied by the isolated R sandbox are also reconciled when evaluating the lock; no package was counted as restored merely because an unrelated user library contained it.

## Clean-project restoration

A newly initialized temporary project outside the repository was created. Its temporary `DESCRIPTION` was generated from the 22 entries in `environment/package_versions.csv` solely so `renv::status()` could evaluate the repository's declared dependency closure. Sandbox mode hid third-party packages from the system R library. Windows binaries were preferred where available; the exact locked `magick` 2.9.0 source package was also successfully built against its version-pinned ImageMagick Windows bundle without changing the package version or installing a global system dependency.

| Test | Result |
|---|---|
| `renv::restore(prompt = FALSE)` | PASS; non-no-op restoration completed |
| Locked packages reconciled | 188/188 |
| Missing locked packages | 0 |
| Version mismatches | 0 |
| Restore error | None |
| Restore elapsed time | 395.94 seconds |
| `renv::status()` | PASS: `No issues found -- the project is in a consistent state.` |
| R scripts parsed without evaluation | 16/16 |
| R parse failures | 0 |

The temporary project library contained 177 package installations; remaining locked R-recommended packages were supplied by the isolated R sandbox. Across the active isolated libraries, every one of the 188 lock records was present at the recorded canonical version.

## Interpretation and limitations

The environment lock is complete, parseable, restorable, and synchronized for R 4.5.2/Bioconductor 3.22 on the validated Windows UCRT platform. This removes the previous empty-lock release blocker.

This test does not establish end-to-end scientific reproducibility because GEO/HPA expression inputs and other third-party matrices are intentionally excluded from the repository. It also does not validate package restoration on macOS or Linux. Those are documented portability and input-availability boundaries, not evidence that scientific analyses were rerun.

## Conclusion

**Environment validation: PASS.**

Repository release status remains **NOT READY** only because no immutable `v1.0.0` release/tag or Zenodo DOI exists and final public-release authorization is still required. The repository owner reports that the private GitHub repository was created and the initial commit was pushed; automated read-only remote verification during this audit was interrupted by a network connection reset, so no independent public-access claim is made.

Scientific lock: `raw_data_modified = FALSE`; 1,115/1,115 raw files unchanged; module count = 24; Tier A/B/C/D = 0/1/21/2.
