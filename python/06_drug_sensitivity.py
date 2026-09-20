#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
06_drug_sensitivity.py
----------------------
Purpose : SKELETON (work in progress). Map the Lactylation-Associated
          transcriptional State (LAS) onto pharmacogenomic resources
          (GDSC / CTRP): score cancer cell lines for an LRG-based ssGSEA-like
          signature and correlate the score with drug IC50 values (Spearman),
          to nominate compounds preferentially active in high-LAS tumours.
Input   : (to be downloaded manually by the user - NOT redistributable here)
          - GDSC cell line expression:
              https://www.cancerrxgene.org/downloads  (e.g. "Cell line RMA
              normalised expression" or the GDSC1/2 bulk RNA-seq release)
          - GDSC drug response:
              https://www.cancerrxgene.org/downloads  ("GDSC1_fitted_dose_response")
          - CTRP v2 (Broad DepMap portal):
              https://depmap.org/portal/download/     ("CTRP v2 viability /
              AUC" + matching expression from DepMap "OmicsExpressionProteinCodingGenesTPMLogp1")
          - data/gene_sets/LRG_geneset.gmt (already in this repository)
Output  : results/tables/drug_sensitivity_correlations.csv
          (planned columns: source, drug, n_cell_lines, spearman_rho, p_value)
Author  : Xingxing Tao
Date    : 2026-09
Run     : python python/06_drug_sensitivity.py \
              --gdsc-expr <path> --gdsc-ic50 <path> [--ctrp-dir <path>] \
              --gene-set data/gene_sets/LRG_geneset.gmt \
              --output-dir results/tables
Status  : IN PROGRESS - this is a scaffold. Major TODOs are flagged inline.
"""

import argparse
import os

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Correlate LAS-like ssGSEA scores with drug sensitivity "
                    "(GDSC/CTRP). SKELETON - see TODOs.")
    parser.add_argument("--gdsc-expr", default=None,
                        help="Path to GDSC cell-line expression matrix "
                             "(download: https://www.cancerrxgene.org/downloads)")
    parser.add_argument("--gdsc-ic50", default=None,
                        help="Path to GDSC fitted dose-response table "
                             "(download: https://www.cancerrxgene.org/downloads)")
    parser.add_argument("--ctrp-dir", default=None,
                        help="Directory with CTRP v2 expression + AUC tables "
                             "(download: https://depmap.org/portal/download/)")
    parser.add_argument("--gene-set", default=os.path.join(
                            "data", "gene_sets", "LRG_geneset.gmt"),
                        help="GMT file with the lactylation-related gene set")
    parser.add_argument("--output-dir", default=os.path.join("results", "tables"),
                        help="Directory for the correlation table")
    return parser.parse_args()


def read_gmt(path: str) -> dict:
    """Parse a GMT file into {set_name: [genes]}."""
    sets = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            sets[parts[0]] = parts[2:]
    return sets


def score_cell_lines(expr: pd.DataFrame, genes: list) -> pd.Series:
    """
    Placeholder LRG signature score for cell lines.

    TODO: Replace this mean-z-score approximation with a proper ssGSEA
    implementation (e.g. call the `gsva` R package via rpy2, or the
    `gseapy.ssgsea` Python implementation) so scores are comparable to the
    patient-level LAS computed in R/03_las_ssgsea.R.
    """
    present = [g for g in genes if g in expr.index]
    if not present:
        raise ValueError("None of the LRG genes found in the expression matrix; "
                         "check that rows are gene symbols.")
    z = expr.loc[present].apply(lambda col: (col - col.mean()) / col.std(ddof=0))
    return z.mean(axis=0)  # mean z-score across LRG genes, per cell line


def main() -> None:
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    gene_sets = read_gmt(args.gene_set)
    lrg_genes = gene_sets["LACTYLATION_RELATED_GENES"]
    print(f"[06] Loaded {len(lrg_genes)} LRG genes from {args.gene_set}")

    results = []

    # ------------------------------------------------------------------ GDSC
    # TODO: implement once the user has downloaded the GDSC tables.
    #   1. expr = pd.read_csv(args.gdsc_expr, index_col=0)   # genes x cell lines
    #   2. ic50 = pd.read_csv(args.gdsc_ic50)                # fitted dose response
    #   3. scores = score_cell_lines(expr, lrg_genes)
    #   4. for each drug: merge scores with -log10(IC50) by cell line, then
    #      rho, p = spearmanr(...); append to results.
    if args.gdsc_expr and args.gdsc_ic50:
        print("[06] GDSC inputs detected, but the GDSC branch is not yet "
              "implemented (see TODO block in source).")
    else:
        print("[06] GDSC inputs not provided; skipping. "
              "Download from https://www.cancerrxgene.org/downloads")

    # ------------------------------------------------------------------ CTRP
    # TODO: same pattern for CTRP v2 (DepMap portal). Note CTRP reports AUC of
    # the dose-response curve rather than IC50; use Spearman correlation of
    # LRG score vs AUC (lower AUC = more sensitive).
    if args.ctrp_dir:
        print("[06] CTRP directory detected, but the CTRP branch is not yet "
              "implemented (see TODO block in source).")
    else:
        print("[06] CTRP directory not provided; skipping. "
              "Download from https://depmap.org/portal/download/")

    # ----------------------------------------------------------------- output
    out_path = os.path.join(args.output_dir, "drug_sensitivity_correlations.csv")
    if results:
        pd.DataFrame(results).to_csv(out_path, index=False)
        print(f"[06] Saved -> {out_path}")
    else:
        # Write an empty table with the planned schema so downstream tooling
        # has a stable contract while this module is under development.
        pd.DataFrame(columns=["source", "drug", "n_cell_lines",
                              "spearman_rho", "p_value"]).to_csv(out_path, index=False)
        print(f"[06] No analyses run yet; wrote empty schema table -> {out_path}")

    print("[06] Reminder: cell-line LRG scores are an approximation of the "
          "patient-level LAS proxy; interpret correlations accordingly.")


if __name__ == "__main__":
    main()
