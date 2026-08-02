# Final numeric consistency check

Overall automated status: **PASS**

This audit cross-checked the locked CSV outputs, integrated 24-module evidence matrix, final DOCX text, figure/table references, citation continuity, administrative sections, forbidden-term context and closure files.

| Check | Status | Evidence |
|---|---:|---|
| GSE102556 samples/donors | PASS | Expected 263 samples and 48 donors. |
| GSE102556 group donors | PASS | Expected 26 MDD and 22 control donors. |
| GSE102556 beta correlation | PASS | Expected 0.9949. |
| GSE102556 singular-fit count | PASS | Expected 3/24 singular or near-zero donor-variance fits. |
| GSE48350 samples/donors | PASS | Expected 253 samples and 84 donors. |
| GSE48350 beta correlation | PASS | Expected 0.9918. |
| GSE48350 FDR counts | PASS | CSV expects naive 4/24 and donor-aware 0/24. |
| AD dataset FDR counts | PASS | CSV expects GSE33000 20/24 and GSE48350 0/24. |
| REML/Hartung-Knapp separation | PASS | CSV expects REML 1/24 and Hartung-Knapp 0/24. |
| Evidence matrix complete | PASS | Expected 24 classified modules. |
| GSE5281 not independently pooled | PASS | GSE5281 must be descriptive/excluded from pooling. |
| GSE167523 not a negative case-control effect | PASS | GSE167523 must be described as an unsupervised preservation test. |
| Raw-data lock stated | PASS | Raw-data immutability statement required. |
| Main figure references | PASS | Found main figure references: [1, 2, 3, 4, 5, 6]. |
| Supplementary table numbering | PASS | Numbering begins at S1 when present. |
| Reference list continuity | PASS | Found 37 sequential reference entries. |
| Citation range valid | PASS | Maximum citation 37; references 37. |
| Forbidden-term context: proves | PASS | No occurrence or every occurrence is explicitly negated/forbidden. |
| Forbidden-term context: demonstrates a liver-to-brain mechanism | PASS | No occurrence or every occurrence is explicitly negated/forbidden. |
| Forbidden-term context: hepatic origin | PASS | No occurrence or every occurrence is explicitly negated/forbidden. |
| Forbidden-term context: liver-mediated | PASS | No occurrence or every occurrence is explicitly negated/forbidden. |
| Forbidden-term context: causal pathway | PASS | No occurrence or every occurrence is explicitly negated/forbidden. |
| Forbidden-term context: mediation effect | PASS | No occurrence or every occurrence is explicitly negated/forbidden. |
| Forbidden-term context: independent replication | PASS | No occurrence or every occurrence is explicitly negated/forbidden. |
| Administrative section: author contributions | PASS | Required administrative material retained. |
| Administrative section: funding | PASS | Required administrative material retained. |
| Administrative section: competing interests | PASS | Required administrative material retained. |
| Administrative section: acknowledgements | PASS | Required administrative material retained. |
| Administrative section: ethics approval | PASS | Required administrative material retained. |
| Administrative section: data availability | PASS | Required administrative material retained. |
| Required closure files | PASS | Missing: none |

## Manual review record

- Effect sizes, standard errors, confidence intervals, P values, BH-FDR values, sample counts, donor counts and model formulas were visually checked in the rendered manuscript tables against their source CSVs.
- Abstract, Results, tables and figure captions were checked for consistent sample/donor terminology and multiple-testing families.
- Conventional REML and Hartung-Knapp were kept distinct; single-dataset significance was not described as replication.
- Cross-tissue association was not interpreted as hepatic origin, transfer, mediation, mechanism or causality.
- Figure references, supplementary numbering, author/administrative sections and reference order were checked after final rendering.
