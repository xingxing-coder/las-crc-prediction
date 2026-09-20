#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
05_shap_interpretation.py
-------------------------
Purpose : Interpret the trained XGBoost LAS model with SHAP: produce a global
          summary (beeswarm) plot and dependence plots for the top-3 most
          important features.
Input   : --model   results/models/xgb_model.joblib  (from 04_train_xgboost.py)
          --input-dir data/processed
            - las_scores.csv, clinical_cohort.csv (same merge as step 04)
Output  : results/figures/shap_summary.png
          results/figures/shap_dependence_<feature>.png (top 3 features)
Author  : Xingxing Tao
Date    : 2026-09
Run     : python python/05_shap_interpretation.py \
              --input-dir data/processed \
              --model results/models/xgb_model.joblib \
              --output-dir results/figures
"""

import argparse
import os

import joblib
import matplotlib
matplotlib.use("Agg")  # headless rendering
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import shap


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="SHAP interpretation of the XGBoost LAS model.")
    parser.add_argument("--input-dir", default=os.path.join("data", "processed"),
                        help="Directory with las_scores.csv and clinical_cohort.csv")
    parser.add_argument("--model", default=os.path.join("results", "models",
                                                        "xgb_model.joblib"),
                        help="Path to the fitted model artefact from step 04")
    parser.add_argument("--output-dir", default=os.path.join("results", "figures"),
                        help="Directory for SHAP figures")
    parser.add_argument("--top-n", type=int, default=3,
                        help="Number of top features for dependence plots")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    # --- Load model artefact and rebuild the feature matrix ------------------
    artefact = joblib.load(args.model)
    reg = artefact["regressor"]
    features = artefact["features"]

    las = pd.read_csv(os.path.join(args.input_dir, "las_scores.csv"))
    clin = pd.read_csv(os.path.join(args.input_dir, "clinical_cohort.csv"))
    df = las.merge(clin, on="patient_id", how="inner")
    X = df[features]
    print(f"[05] Feature matrix: {X.shape[0]} patients x {X.shape[1]} features")

    # Transform with the fitted preprocessing pipeline so SHAP sees the exact
    # feature space the booster was trained on.
    prep = reg.named_steps["prep"]
    booster = reg.named_steps["model"]
    X_enc = prep.transform(X)
    if hasattr(X_enc, "toarray"):  # sparse one-hot matrix -> dense for plotting
        X_enc = X_enc.toarray()
    try:
        feature_names = list(prep.get_feature_names_out())
    except AttributeError:  # older scikit-learn
        feature_names = [f"f{i}" for i in range(X_enc.shape[1])]
    X_enc = pd.DataFrame(X_enc, columns=feature_names)

    # --- SHAP values ----------------------------------------------------------
    explainer = shap.Explainer(booster, X_enc)
    shap_values = explainer(X_enc)

    # Summary (beeswarm) plot.
    plt.figure()
    shap.summary_plot(shap_values, X_enc, show=False)
    summary_path = os.path.join(args.output_dir, "shap_summary.png")
    plt.tight_layout()
    plt.savefig(summary_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"[05] Saved -> {summary_path}")

    # Dependence plots for the top-N features by mean |SHAP|.
    importance = np.abs(shap_values.values).mean(axis=0)
    top_features = [feature_names[i] for i in np.argsort(importance)[::-1][:args.top_n]]
    for feat in top_features:
        plt.figure()
        shap.dependence_plot(feat, shap_values.values, X_enc, show=False)
        dep_path = os.path.join(args.output_dir, f"shap_dependence_{feat}.png")
        plt.tight_layout()
        plt.savefig(dep_path, dpi=300, bbox_inches="tight")
        plt.close()
        print(f"[05] Saved -> {dep_path}")

    print("[05] Done.")


if __name__ == "__main__":
    main()
