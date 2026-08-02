# Matched random-module negative-control design

Generated: 2026-08-02 02:31:35 +0800

- Fixed seed: 20260802
- Random sets per real module: B=2000
- Common universe: 17000 genes measured in GSE135251, strict-mapped GSE48350 and strict-mapped GSE33000.
- Exact matching: common-universe module size and joint GSE135251 mean-expression/IQR decile counts.
- GSVA/ssGSEA was applied jointly to all real and random sets; scores were standardized separately within each target dataset.
- GSE48350: `score_z ~ diagnosis + region + age_z + sex + (1 | donor_id)`; GSE33000: `score_z ~ diagnosis + age_z + sex`.
- Empirical P=(1+at-least-as-extreme null count)/(B+1); lower difference and higher magnitude/consistency/support are the stated one-sided tails.
- BH correction spans the 24 real modules separately for every empirical metric.
