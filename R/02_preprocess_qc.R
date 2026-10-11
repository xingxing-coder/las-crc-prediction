# ==============================================================================
# 02_preprocess_qc.R
# ------------------------------------------------------------------------------
# Purpose : Preprocess and quality-control the TCGA-COAD/READ cohort:
#           apply the sample inclusion flow (primary tumour only, one sample
#           per patient, complete stage/MSI/survival information), normalise
#           raw counts to log2-CPM with edgeR, and export the analysis-ready
#           expression matrix plus a tidy clinical table.
# Input   : data/raw/tcga_coadread_se.rds        (from R/01_download_tcga.R)
#           data/raw/tcga_coadread_clinical.rds  (from R/01_download_tcga.R)
# Output  : data/processed/expression_log2cpm.rds (gene x sample matrix,
#                                                 rows = gene symbol)
#           data/processed/clinical_cohort.csv    (columns include patient_id,
#                                                  age, gender, tnm_stage,
#                                                  msi_status, os_time, os_event)
#           data/processed/cohort_flow.csv        (step, n_remaining, note)
# Author  : Xingxing Tao
# Date    : 2026-09
# Run     : From the repository root:  Rscript R/02_preprocess_qc.R
# Notes   : The final cohort size (approx. 620 patients expected) varies
#           slightly between GDC data freezes. The download date and data
#           version are recorded in data/raw/download_log.txt; always cite
#           them when reporting sample counts.
#           MSI status is not part of the standard GDC clinical XML/indexed
#           record; this script looks for an MSI annotation column and, if
#           absent, inserts NA with a TODO to merge a curated pan-cancer MSI
#           call set (see inline comment).
# ==============================================================================

## ---- 0. Dependencies --------------------------------------------------------
library(SummarizedExperiment)
library(edgeR)
library(data.table)

## ---- 1. Configuration -------------------------------------------------------
se_in     <- file.path("data", "raw", "tcga_coadread_se.rds")
clin_in   <- file.path("data", "raw", "tcga_coadread_clinical.rds")

proc_dir  <- file.path("data", "processed")
dir.create(proc_dir, showWarnings = FALSE, recursive = TRUE)

expr_out  <- file.path(proc_dir, "expression_log2cpm.rds")
cohort_out <- file.path(proc_dir, "clinical_cohort.csv")
flow_out  <- file.path(proc_dir, "cohort_flow.csv")

# Cohort inclusion flow tracker: step / n_remaining / note
flow <- data.frame(
  step        = character(),
  n_remaining = integer(),
  note        = character(),
  stringsAsFactors = FALSE
)
record_step <- function(step, n, note) {
  message(sprintf("[02] %-55s n = %d", step, n))
  flow[nrow(flow) + 1, ] <<- list(step, as.integer(n), note)
}

## ---- 2. Load raw data --------------------------------------------------------
se   <- readRDS(se_in)
clin <- as.data.table(readRDS(clin_in))

counts <- assay(se, "unstranded")          # STAR - Counts raw counts
barcodes <- colnames(counts)               # full TCGA aliquot barcodes

## ---- 3. Inclusion flow ------------------------------------------------------
# Step 1. All RNA-seq samples (tumour + adjacent normal)
record_step("1_raw_rnaseq_samples", ncol(counts),
            "All TCGA-COAD/READ STAR-Counts samples (tumour + normal)")

# Step 2. Keep primary tumour only.
# TCGA barcode positions 14-15 encode the sample type; "01" = Primary Tumor.
# Equivalent check: colData(se)$sample_type == "Primary Tumor".
sample_type_code <- substr(barcodes, 14, 15)
keep_tumor <- sample_type_code == "01"
stopifnot(all(colData(se)$sample_type[keep_tumor] %in%
                c("Primary Tumor", "Primary Blood Derived Cancer - Bone Marrow")))
counts <- counts[, keep_tumor, drop = FALSE]
record_step("2_primary_tumor_only", ncol(counts),
            "Removed normal tissue (kept barcode positions 14-15 == '01')")

# Step 3. De-duplicate patients: patient_id = first 12 chars of the barcode
# (e.g. TCGA-AA-1234); keep the first aliquot per patient.
patient_id <- substr(colnames(counts), 1, 12)
keep_first <- !duplicated(patient_id)
counts <- counts[, keep_first, drop = FALSE]
patient_id <- patient_id[keep_first]
colnames(counts) <- patient_id
record_step("3_deduplicated_patients", ncol(counts),
            "One aliquot per patient (first aliquot kept; key = barcode[1:12])")

## ---- 4. Build tidy clinical table ------------------------------------------
# Map GDC clinical fields to analysis columns.
clin[, patient_id := submitter_id]

# Survival: days_to_death for deceased, days_to_last_follow_up otherwise.
clin[, os_event := as.integer(vital_status == "Dead")]
os_time <- ifelse(!is.na(clin$days_to_death),
                  clin$days_to_death,
                  clin$days_to_last_follow_up)
clin[, os_time := as.numeric(os_time)]

# TNM stage / gender / age harmonisation across TCGAbiolinks versions.
# Column names in GDCquery_clinic() output vary between TCGAbiolinks
# releases, so resolve columns from a candidate list and fall back to NA
# with a warning (instead of crashing) when a column is absent.
pick_col <- function(dt, candidates, label) {
  hit <- intersect(candidates, names(dt))[1]
  if (is.na(hit)) {
    hit <- grep(paste0("^(", paste(candidates, collapse = "|"), ")$"),
                names(dt), ignore.case = TRUE, value = TRUE)[1]
  }
  if (is.na(hit)) {
    message("[02] WARNING: no ", label, " column found; set to NA. ",
            "Available columns: ", paste(names(dt), collapse = ", "))
    return(NA_character_)
  }
  hit
}

stage_col  <- pick_col(clin, c("ajcc_pathologic_stage", "figo_stage"),
                       "TNM stage")
gender_col <- pick_col(clin, c("gender", "sex_at_birth", "sex"), "gender")
age_col    <- pick_col(clin, c("age_at_diagnosis", "age_at_index"), "age")

clin[, tnm_stage := if (is.na(stage_col))  NA_character_ else get(stage_col)]
clin[, gender    := if (is.na(gender_col)) NA_character_ else get(gender_col)]
clin[, age       := if (is.na(age_col))    NA_real_      else as.numeric(get(age_col))]
if (!all(is.na(clin$age)) && max(clin$age, na.rm = TRUE) > 150)
  clin[, age := age / 365.25]  # days -> years

# MSI status ---------------------------------------------------------------
# TODO(msi): The default GDCquery_clinic(type = "clinical") table has no MSI
# column. Merge a curated MSI call set here before final analysis, e.g.:
#   Bonneville et al., JCO Precis Oncol 2017 (pan-cancer MANTIS MSI calls), or
#   cBioPortal COAD/READ "MSI Status" clinical attribute, or the
#   TCGA PanCanAtlas MSI subtype annotation.
# Until merged, msi_status is NA and step 4 of the flow uses a placeholder
# completeness rule; DO NOT treat the resulting n as the final cohort size.
msi_col <- grep("msi", names(clin), ignore.case = TRUE, value = TRUE)
if (length(msi_col) > 0) {
  clin[, msi_status := get(msi_col[1])]
} else {
  clin[, msi_status := NA_character_]
  message("[02] NOTE: no MSI column found in clinical table; msi_status set to NA.",
          " Merge external MSI calls before final analysis (see TODO above).")
}

clinical_cohort <- clin[, .(patient_id, age, gender, tnm_stage,
                            msi_status, os_time, os_event)]

# Restrict clinical table to patients present in the expression matrix.
clinical_cohort <- clinical_cohort[patient_id %in% colnames(counts)]

## ---- 5. Step 4 of the flow: complete stage / MSI / survival ----------------
complete_flag <- complete.cases(
  clinical_cohort[, .(tnm_stage, os_time, os_event)]
)
# MSI completeness is only enforced once real MSI calls are merged (TODO above).
if (!all(is.na(clinical_cohort$msi_status))) {
  complete_flag <- complete_flag & !is.na(clinical_cohort$msi_status)
}

clinical_cohort <- clinical_cohort[complete_flag]
counts <- counts[, clinical_cohort$patient_id, drop = FALSE]
record_step("4_complete_clinical", ncol(counts),
            paste0("Removed patients with missing TNM stage / MSI / survival;",
                   " approx. 620 expected (varies by GDC data freeze)"))

## ---- 6. Normalisation: log2-CPM ---------------------------------------------
# Gene symbols as row names (SummarizedExperiment carries them in rowData).
gene_symbols <- rowData(se)$gene_name
rownames(counts) <- make.unique(gene_symbols)

# Filter very lowly expressed genes, then edgeR log2-CPM.
keep_genes <- rowSums(counts >= 10) >= ceiling(0.1 * ncol(counts))
counts <- counts[keep_genes, , drop = FALSE]
log2cpm <- edgeR::cpm(counts, log = TRUE, prior.count = 1)

saveRDS(log2cpm, expr_out)
fwrite(clinical_cohort, cohort_out)
fwrite(flow, flow_out)

message("[02] Saved expression matrix  -> ", expr_out,
        " (", nrow(log2cpm), " genes x ", ncol(log2cpm), " patients)")
message("[02] Saved clinical cohort    -> ", cohort_out)
message("[02] Saved cohort flow table  -> ", flow_out)
message("[02] Done. Reminder: cohort size depends on the GDC data freeze;",
        " cite data/raw/download_log.txt when reporting n.")
