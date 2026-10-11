# ==============================================================================
# 01_download_tcga.R
# ------------------------------------------------------------------------------
# Purpose : Download TCGA-COAD and TCGA-READ RNA-seq (STAR - Counts) gene
#           expression data and clinical data from the NCI Genomic Data
#           Commons (GDC) via the TCGAbiolinks package, and save them as
#           local RDS objects for downstream processing.
# Input   : None (data are fetched from the GDC API; internet connection and
#           ~4 GB of free disk space required).
# Output  : data/raw/tcga_coadread_se.rds        (SummarizedExperiment object)
#           data/raw/tcga_coadread_clinical.rds  (clinical data.frame)
#           data/raw/download_log.txt            (download date + GDC version)
# Author  : Xingxing Tao
# Date    : 2026-09
# Run     : From the repository root:
#             Rscript R/01_download_tcga.R
#           or inside an interactive R session:
#             source("R/01_download_tcga.R")
# Notes   : TCGA data are subject to the GDC Data Use Agreement. Raw data must
#           NOT be redistributed or committed to this repository (see
#           .gitignore). Record the download date and GDC data release version
#           because sample counts vary slightly between data freezes.
# ==============================================================================

## ---- 0. Dependencies --------------------------------------------------------
# Install TCGAbiolinks once via Bioconductor if needed:
#   if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
#   BiocManager::install("TCGAbiolinks")

library(TCGAbiolinks)
library(SummarizedExperiment)

## ---- 1. Configuration -------------------------------------------------------
raw_dir    <- file.path("data", "raw")
gdc_cache  <- file.path(raw_dir, "GDCdata")   # local GDC download cache
dir.create(raw_dir, showWarnings = FALSE, recursive = TRUE)

projects   <- c("TCGA-COAD", "TCGA-READ")

se_out     <- file.path(raw_dir, "tcga_coadread_se.rds")
clin_out   <- file.path(raw_dir, "tcga_coadread_clinical.rds")
log_out    <- file.path(raw_dir, "download_log.txt")

## ---- 2. Query, download and prepare transcriptome data ----------------------
message("[01] Querying GDC for TCGA-COAD/READ STAR - Counts expression data ...")

query_exp <- GDCquery(
  project       = projects,
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

message("[01] Downloading files (this can take a while on first run) ...")
GDCdownload(
  query     = query_exp,
  directory = gdc_cache
)

message("[01] Preparing SummarizedExperiment object ...")
se_coadread <- GDCprepare(
  query     = query_exp,
  directory = gdc_cache
)

saveRDS(se_coadread, se_out)
message("[01] Saved expression SummarizedExperiment -> ", se_out)
message("[01] Dimensions: ", nrow(se_coadread), " genes x ",
        ncol(se_coadread), " samples")

## ---- 3. Download clinical data ----------------------------------------------
message("[01] Downloading clinical data via GDCquery_clinic ...")

clinical_list <- lapply(projects, function(p) {
  clin <- GDCquery_clinic(project = p, type = "clinical")
  clin$project_id <- p
  clin
})

# NOTE: GDCquery_clinic returns tables whose columns differ between projects
# (e.g. COAD vs READ), so a plain rbind() fails with "variables' column
# number mismatch". rbindlist(fill = TRUE) binds rows by column *name* and
# pads missing columns with NA. data.table is a TCGAbiolinks dependency.
clinical <- data.table::rbindlist(clinical_list, fill = TRUE)
clinical <- as.data.frame(clinical)

saveRDS(clinical, clin_out)
message("[01] Saved clinical data -> ", clin_out,
        " (", nrow(clinical), " patients)")

## ---- 4. Record provenance ----------------------------------------------------
writeLines(
  c(
    paste0("download_date: ", Sys.Date()),
    paste0("projects: ", paste(projects, collapse = ", ")),
    paste0("workflow_type: STAR - Counts"),
    paste0("TCGAbiolinks_version: ", as.character(packageVersion("TCGAbiolinks"))),
    paste0("R_version: ", R.version.string),
    "note: Sample counts vary slightly between GDC data freezes; record this",
    "      date when reporting cohort sizes in the manuscript/application."
  ),
  con = log_out
)

message("[01] Done. Provenance written to ", log_out)
