SUPPLEMENTARY / REVISION VALIDATION DEPOSIT
Human Islet Transcriptomics in Type 2 Diabetes Reveals Conserved Dysfunction Axes and a Candidate Gene Panel

CONTENTS
- Supplementary_Material.docx: Supplementary Figures S1-S2 and Tables S1-S7. Supplementary Table S4 is a single combined group-level and donor-level clinical/demographic table; no S4B split table remains.
- analysis_scripts/study_level_crossfit_modules.py: rotating three-cohort module definition / held-out fourth-cohort validation plus matched gene-set null.
- analysis_scripts/rds_miniparser.py: minimal reader used to access stored RDS expression matrices in the cross-fitting audit runtime.
- analysis_scripts/nested_5fold_ml.py: fully nested stratified five-fold CV with DEG filtering and feature selection inside outer training folds.
- analysis_scripts/make_figure7_nested.py: revised Figure 7 generation.
- analysis_scripts/delong_auc_ci.py: DeLong confidence interval for GSE50244 external AUC.
- results/: donor-level GSE221156 data, cross-fitted module outputs, matched-null outputs, and nested-CV predictions/metrics/stability.

PRIMARY REVISION VALIDATION RESULTS
- Cross-fitted IsletDysfunctionScore d: GSE25724 1.012; GSE20966 1.359; GSE38642 0.390; GSE164416 1.129.
- Pooled held-out Hedges g: 0.905 (95% CI 0.501-1.309), I2=5.5%.
- Matched-null empirical p for composite: 0.3973, 0.0245, 0.4188, 0.0020, respectively.
- Fully nested 5-fold pooled OOF AUC: 0.960 (bootstrap 95% CI 0.900-0.997).
- External GSE50244 AUC: 0.907 (DeLong 95% CI 0.818-0.996).

The public source datasets are GEO GSE25724, GSE20966, GSE38642, GSE164416, GSE50244, and GSE221156. Core pipeline code is maintained in the public repositories cited in the manuscript. This deposit captures the exact supplementary and revision-validation analyses used in the resubmitted manuscript.
