# las-crc-prediction

**Knowledge-enhanced machine learning to predict lactylation-associated transcriptional states (LAS) in colorectal cancer (TCGA-COAD/READ).**

![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)
![R](https://img.shields.io/badge/R-%3E%3D4.3-276DC3?logo=r)
![Python](https://img.shields.io/badge/Python-%3E%3D3.10-3776AB?logo=python)
![Status](https://img.shields.io/badge/status-work%20in%20progress-orange)

Companion open-source repository for a PhD research proposal: we derive a
lactylation-associated transcriptional score from TCGA colorectal cancer
bulk RNA-seq, model it with interpretable machine learning, and map it to
drug-sensitivity resources.

## Pipeline overview

```
 TCGA-COAD/READ            QC + inclusion flow          LAS scoring
 (GDC, STAR-Counts)  --->  (edgeR log2-CPM, n=600)  ---> (ssGSEA on LRG set)
        R/01                    R/02                        R/03
                                                              |
                                                              v
 Drug mapping          <---   SHAP interpretation  <---  XGBoost model
 (GDSC/CTRP, TODO)            (python/05)                 (python/04)
        python/06
```

| Step | Script | Input | Output |
|------|--------|-------|--------|
| 1 | `R/01_download_tcga.R` | GDC API | `data/raw/tcga_coadread_se.rds`, `tcga_coadread_clinical.rds` |
| 2 | `R/02_preprocess_qc.R` | raw RDS | `data/processed/expression_log2cpm.rds`, `clinical_cohort.csv`, `cohort_flow.csv` |
| 3 | `R/03_las_ssgsea.R` | expression matrix + `data/gene_sets/LRG_geneset.gmt` | `data/processed/las_scores.csv`, `results/tables/construct_validity.csv` |
| 4 | `python/04_train_xgboost.py` | `las_scores.csv` + `clinical_cohort.csv` | `results/tables/model_metrics.csv`, `results/models/xgb_model.joblib` |
| 5 | `python/05_shap_interpretation.py` | fitted model + features | `results/figures/shap_summary.png`, `shap_dependence_*.png` |
| 6 | `python/06_drug_sensitivity.py` | GDSC/CTRP downloads (manual) | `results/tables/drug_sensitivity_correlations.csv` |

The join key across all files is `patient_id` (first 12 characters of the
TCGA barcode, e.g. `TCGA-AA-1234`).

## Quickstart

```bash
# 1. R environment (Bioconductor): see environment-r.txt
Rscript R/01_download_tcga.R     # download TCGA-COAD/READ from GDC
Rscript R/02_preprocess_qc.R     # QC, inclusion flow, log2-CPM
Rscript R/03_las_ssgsea.R        # LAS scoring via ssGSEA

# 2. Python environment
pip install -r requirements.txt
python python/04_train_xgboost.py --input-dir data/processed --output-dir results
python python/05_shap_interpretation.py
python python/06_drug_sensitivity.py   # skeleton; requires manual GDSC/CTRP downloads
```

## Data availability

**No TCGA data are included in this repository.** TCGA raw and processed data
must be downloaded individually from the
[NCI Genomic Data Commons](https://portal.gdc.cancer.gov/) and are subject to
the GDC Data Use Agreement; redistribution is not permitted. This repository
contains code only. `data/raw/`, `data/processed/`, and `results/` are
git-ignored by design. GDSC and CTRP data likewise must be obtained from
their respective portals.

## Methodological note

The **LAS score is a proxy of lactylation-associated transcriptional activity**
inferred from mRNA co-expression of lactylation-related genes (LRGs). It is
**not** a direct measurement of protein lactylation levels (e.g. pan-Kla
immunoblotting or lactylome mass spectrometry). This interpretation follows
the lactylation-related-gene framework popularised in recent literature
(e.g. *Cell Bioscience*, 2025). The gene set shipped in
`data/gene_sets/LRG_geneset.gmt` is a ~40-gene starter scaffold and must be
expanded to the full literature-derived set (e.g. the 332-gene compendium)
before formal analysis.

## Preliminary results (GDC download: October 2026; exact freeze in `data/raw/download_log.txt`)

Cohort inclusion flow (TCGA-COAD/READ, STAR-Counts):

| Step | n |
|------|---|
| RNA-seq samples (tumour + adjacent normal) | 701 |
| Primary tumour only (barcode positions 14-15 == "01") | 647 |
| Deduplicated patients (one aliquot per patient) | 624 |
| Complete TNM stage + survival data (modeling cohort) | **600** |

Construct validity of LAS (Spearman correlation with LAS, N = 600):

| Gene | rho | p value | Interpretation |
|------|-----|---------|----------------|
| ALDOA | 0.66 | 3.3e-75 | Glycolytic flux (lactate substrate supply) |
| ENO1 | 0.47 | 1.8e-34 | Glycolytic flux |
| PFKP | 0.47 | 1.1e-33 | Glycolytic flux |
| PFKM | 0.18 | 1.1e-05 | Glycolytic flux |
| EP300 | -0.31 | 6.8e-15 | Lactylation writer (mRNA level) |

LAS correlates strongly with the glycolytic / lactate-supply transcriptional
programme but not with EP300 mRNA abundance. This is consistent with LAS
capturing substrate-driven lactylation-associated transcription rather than
writer expression — enzymatic lactyltransferase activity is not reflected by
EP300 mRNA levels — and reinforces that LAS must be interpreted as a
transcriptional proxy (see Methodological note).

## Status

Work in progress. Completed: TCGA-COAD/READ download, preprocessing and QC
pipeline, ssGSEA LAS scoring with construct-validity checks (n = 600). In
progress: model tuning, SHAP interpretation, GDSC/CTRP drug-sensitivity
mapping.

## Contact

Xingxing Tao — prospective PhD applicant (data science / biomedical
informatics). Issues and pull requests are welcome.

## License

MIT License — see [LICENSE](LICENSE).
