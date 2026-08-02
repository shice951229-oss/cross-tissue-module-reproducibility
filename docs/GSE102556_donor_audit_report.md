# GSE102556 donor audit report

Generated: 2026-08-02 00:53:54 +0800

## Identity evidence

The donor identifier is the explicit numeric prefix before the colon in each GEO sample title. It was cross-checked against the project `subject_id`, the prefix of the expression-table sample name, and presence of the complete sample name among the downloaded FPKM columns. No donor was inferred from ordering or diagnosis.

- Total samples: 263
- Independent donors: 48
- Missing donor identifiers: 0
- Samples with any donor conflict: 0
- Donor-aware analytic samples: 263

## Donors by diagnosis

 diagnosis donors
   Control     22
       MDD     26

## Regions

                   brain_region samples donors
      anterior_cingulate_cortex      28     28
                anterior_insula      48     48
 dorsolateral_prefrontal_cortex      48     48
              nucleus_accumbens      48     48
           orbitofrontal_cortex      48     48
                      subiculum      43     43

## Regions contributed per donor


 4  5  6 
 4 17 27 

## Covariate audit

 covariate missing_n missing_fraction unique_nonmissing                                                                                                       values_by_diagnosis
       age         0                0                33                       Control                           MDD  "n=122 mean=47.639 sd=16.482" "n=141 mean=45.688 sd=13.156" 
       sex         0                0                 2                                                                                         Control     MDD  "72;50" "67;74" 
       RIN         0                0                45                             Control                         MDD  "n=122 mean=6.996 sd=0.916" "n=141 mean=7.200 sd=0.974" 
       PMI         0                0                40                       Control                           MDD  "n=122 mean=26.707 sd=21.733" "n=141 mean=27.330 sd=16.824" 

## Primary model

`score_z ~ diagnosis + brain_region + age_z + sex + RIN_z + PMI_z + (1 | donor_id)`

The 24 overall diagnosis effects form one BH family. The naive comparison uses the identical analytic sample set and fixed-effect design but omits the donor random intercept. Region-specific analyses are descriptive; their primary correction spans all 144 tests, with within-region 24-test BH values supplied only as sensitivity results.

Effect correlation (naive vs donor-aware): 0.9949
Naive FDR<0.05 modules: 0
Donor-aware FDR<0.05 modules: 0
Median naive SE / donor-aware SE: 0.6988
