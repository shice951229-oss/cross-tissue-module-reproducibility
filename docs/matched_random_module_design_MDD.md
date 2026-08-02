# Matched random-module negative-control design: MDD

- Fixed seed: 20260802; B=2000 random sets per observed module and tissue.
- Matching was conducted separately for MDD blood and MDD brain using module size, GSE135251 mean-expression decile, GSE135251 IQR decile, and target detectability.
- GSE98793 used the final strict unambiguous/IQR probe mapping and the prespecified Welch/Cohen-d analysis.
- GSE102556 used the verified 48-donor structure and the donor-aware mixed model.
- Blood and brain were tested and reported separately; no cross-tissue pooled MDD estimate was calculated.
- Empirical P values use (1 + at-least-as-extreme null count)/(B + 1); BH correction spans 24 modules within each tissue.
