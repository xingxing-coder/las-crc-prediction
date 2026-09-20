# ==============================================================================
# 03_las_ssgsea.R
# ------------------------------------------------------------------------------
# Purpose : Score each TCGA-COAD/READ patient for the Lactylation-Associated
#           transcriptional State (LAS) using single-sample GSEA (ssGSEA) on
#           the lactylation-related gene (LRG) set, and run construct-validity
#           checks against EP300 (a writer of protein lactylation) and core
#           glycolysis genes (PFKM / PFKP / ALDOA / ENO1), which supply the
#           lactate substrate of lactylation.
# Input   : data/processed/expression_log2cpm.rds (from R/02_preprocess_qc.R)
#           data/gene_sets/LRG_geneset.gmt        (lactylation-related genes)
# Output  : data/processed/las_scores.csv          (columns: patient_id, LAS)
#           results/tables/construct_validity.csv  (Spearman correlations of
#                                                   LAS vs EP300 / glycolysis
#                                                   gene expression)
# Author  : Xingxing Tao
# Date    : 2026-09
# Run     : From the repository root:  Rscript R/03_las_ssgsea.R
# Notes   : IMPORTANT METHODOLOGICAL CAVEAT --------------------------------
#           The LAS score is a PROXY of lactylation-associated transcriptional
#           activity inferred from mRNA co-expression of lactylation-related
#           genes. It is NOT a direct measurement of protein lactylation
#           levels (e.g. pan-Kla immunoblotting or lactylome mass
#           spectrometry). All downstream interpretations must reflect this.
#           ------------------------------------------------------------------
#           Gene-set size: the starter set in data/gene_sets/LRG_geneset.gmt
#           (~40 genes) is a minimal scaffold. For the formal analysis it must
#           be expanded to the literature-derived set cited in the research
#           proposal (e.g. the 332-gene lactylation-related compendium,
#           Cell Biosci 2025).
# ==============================================================================

## ---- 0. Dependencies --------------------------------------------------------
library(GSVA)        # >= 1.50 for the *Param() API used below
library(data.table)

## ---- 1. Configuration -------------------------------------------------------
expr_in  <- file.path("data", "processed", "expression_log2cpm.rds")
gmt_in   <- file.path("data", "gene_sets", "LRG_geneset.gmt")

las_out  <- file.path("data", "processed", "las_scores.csv")
valid_dir <- file.path("results", "tables")
dir.create(valid_dir, showWarnings = FALSE, recursive = TRUE)
valid_out <- file.path(valid_dir, "construct_validity.csv")

## ---- 2. Load data ------------------------------------------------------------
mat <- readRDS(expr_in)                 # genes x patients (log2-CPM)
stopifnot(is.matrix(mat))

gene_sets <- GSEABase::getGmt(gmt_in)   # GeneSetCollection from the GMT file
message("[03] Loaded gene set(s): ", paste(names(gene_sets), collapse = ", "))
message("[03] Expression matrix: ", nrow(mat), " genes x ",
        ncol(mat), " patients")

## ---- 3. ssGSEA scoring (GSVA >= 1.50 new API) --------------------------------
# New parameter-object API (recommended):
param <- ssgseaParam(exprData = mat, geneSets = gene_sets)
las_mat <- gsva(param)

# --- Legacy API (GSVA < 1.50), kept as a fallback: ---------------------------
# las_mat <- gsva(expr = mat, gset.idx.list = gene_sets, method = "ssgsea",
#                 kcdf = "Gaussian", verbose = TRUE)
# -----------------------------------------------------------------------------

las_scores <- data.table(
  patient_id = colnames(las_mat),
  LAS        = as.numeric(las_mat["LACTYLATION_RELATED_GENES", ])
)
fwrite(las_scores, las_out)
message("[03] Saved LAS scores -> ", las_out,
        " (", nrow(las_scores), " patients)")

## ---- 4. Construct-validity checks --------------------------------------------
# If LAS genuinely captures lactylation-associated transcriptional activity,
# it should correlate positively with the lactylation writer EP300 and with
# glycolytic flux genes (PFKM/PFKP/ALDOA/ENO1) that drive lactate production.
check_genes <- c("EP300", "PFKM", "PFKP", "ALDOA", "ENO1")
check_genes <- intersect(check_genes, rownames(mat))

validity <- rbindlist(lapply(check_genes, function(g) {
  ct <- suppressWarnings(
    cor.test(las_scores$LAS, mat[g, las_scores$patient_id],
             method = "spearman", exact = FALSE)
  )
  data.table(
    gene      = g,
    rho       = unname(ct$estimate),
    p_value   = ct$p.value,
    rationale = switch(g,
      "EP300" = "Lactylation writer (histone lysine lactyltransferase)",
      "PFKM"  = "Glycolytic flux -> lactate substrate for lactylation",
      "PFKP"  = "Glycolytic flux -> lactate substrate for lactylation",
      "ALDOA" = "Glycolytic flux -> lactate substrate for lactylation",
      "ENO1"  = "Glycolytic flux -> lactate substrate for lactylation",
      "Construct validity check"
    )
  )
}))

fwrite(validity, valid_out)
message("[03] Construct-validity correlations -> ", valid_out)
print(validity)

message("[03] Done. Reminder: LAS is a proxy of lactylation-associated ",
        "transcriptional activity, not a direct measure of protein lactylation.")
