#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
04_train_xgboost.py
-------------------
Purpose : Train a knowledge-enhanced XGBoost regression model to predict the
          Lactylation-Associated transcriptional State (LAS) score from
          clinical features, and evaluate it with regression and binary-
          classification metrics.
Input   : --input-dir (default: data/processed)
            - las_scores.csv       (patient_id, LAS; from R/03_las_ssgsea.R)
            - clinical_cohort.csv  (patient_id, age, gender, tnm_stage,
                                    msi_status, os_time, os_event;
                                    from R/02_preprocess_qc.R)
Output  : results/tables/model_metrics.csv  (R2/MAE/RMSE + AUC/sensitivity/
                                             specificity for the median-split
                                             binary LAS task)
          results/models/xgb_model.joblib   (fitted best estimator)
Author  : Xingxing Tao
Date    : 2026-09
Run     : python python/04_train_xgboost.py \
              --input-dir data/processed --output-dir results
Notes   : LAS is a proxy of lactylation-associated transcriptional activity,
          not a direct measure of protein lactylation levels.
"""

import argparse
import os

import joblib
import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.metrics import (accuracy_score, confusion_matrix, mean_absolute_error,
                             mean_squared_error, r2_score, roc_auc_score)
from sklearn.model_selection import GridSearchCV, train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from xgboost import XGBClassifier, XGBRegressor

RANDOM_STATE = 42

# ---------------------------------------------------------------------------
# Feature schema.
# TODO(haematology): the research proposal specifies "clinical + haematological
# indicators" as predictors. TCGA clinical data do NOT contain routine blood
# counts (WBC, neutrophils, lymphocytes, NLR, platelets, ...). Options:
#   (1) drop haematological features for the TCGA analysis and state this
#       limitation explicitly;
#   (2) validate the haematology-enhanced model on an independent in-house /
#       public cohort where blood tests are available.
# Until then, the model is demonstrated with the available clinical features.
# ---------------------------------------------------------------------------
NUMERIC_FEATURES = ["age", "os_time"]
CATEGORICAL_FEATURES = ["gender", "tnm_stage", "msi_status"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Train XGBoost model predicting the LAS score in TCGA-COAD/READ."
    )
    parser.add_argument("--input-dir", default=os.path.join("data", "processed"),
                        help="Directory with las_scores.csv and clinical_cohort.csv")
    parser.add_argument("--output-dir", default="results",
                        help="Directory for metrics table and fitted model")
    return parser.parse_args()


def load_data(input_dir: str) -> pd.DataFrame:
    las = pd.read_csv(os.path.join(input_dir, "las_scores.csv"))
    clin = pd.read_csv(os.path.join(input_dir, "clinical_cohort.csv"))
    df = las.merge(clin, on="patient_id", how="inner")
    if df.empty:
        raise ValueError("Merged dataset is empty - check that both input files "
                         "share the same patient_id values (TCGA barcode[1:12]).")
    print(f"[04] Merged dataset: {df.shape[0]} patients x {df.shape[1]} columns")
    return df


def build_preprocessor() -> ColumnTransformer:
    numeric = Pipeline(steps=[
        ("impute", SimpleImputer(strategy="median")),
        ("scale", StandardScaler()),
    ])
    categorical = Pipeline(steps=[
        ("impute", SimpleImputer(strategy="most_frequent")),
        ("onehot", OneHotEncoder(handle_unknown="ignore")),
    ])
    return ColumnTransformer(transformers=[
        ("num", numeric, NUMERIC_FEATURES),
        ("cat", categorical, CATEGORICAL_FEATURES),
    ])


def tune_regressor(X_train: pd.DataFrame, y_train: pd.Series) -> Pipeline:
    pipe = Pipeline(steps=[
        ("prep", build_preprocessor()),
        ("model", XGBRegressor(random_state=RANDOM_STATE, n_jobs=-1,
                               objective="reg:squarederror")),
    ])
    grid = {
        "model__n_estimators": [200, 500],
        "model__max_depth": [3, 5],
        "model__learning_rate": [0.01, 0.1],
    }
    search = GridSearchCV(pipe, grid, cv=5, scoring="r2", n_jobs=-1)
    search.fit(X_train, y_train)
    print(f"[04] Best hyper-parameters: {search.best_params_} "
          f"(CV R2 = {search.best_score_:.3f})")
    return search.best_estimator_


def tune_classifier(X_train: pd.DataFrame, y_train: pd.Series) -> Pipeline:
    pipe = Pipeline(steps=[
        ("prep", build_preprocessor()),
        ("model", XGBClassifier(random_state=RANDOM_STATE, n_jobs=-1,
                                eval_metric="logloss")),
    ])
    grid = {
        "model__n_estimators": [200, 500],
        "model__max_depth": [3, 5],
        "model__learning_rate": [0.01, 0.1],
    }
    search = GridSearchCV(pipe, grid, cv=5, scoring="roc_auc", n_jobs=-1)
    search.fit(X_train, y_train)
    print(f"[04] Best classifier parameters: {search.best_params_} "
          f"(CV AUC = {search.best_score_:.3f})")
    return search.best_estimator_


def main() -> None:
    args = parse_args()
    tables_dir = os.path.join(args.output_dir, "tables")
    models_dir = os.path.join(args.output_dir, "models")
    os.makedirs(tables_dir, exist_ok=True)
    os.makedirs(models_dir, exist_ok=True)

    df = load_data(args.input_dir)
    X = df[NUMERIC_FEATURES + CATEGORICAL_FEATURES]
    y_reg = df["LAS"]

    # 7:3 split (70% train / 30% test)
    X_train, X_test, y_train, y_test = train_test_split(
        X, y_reg, test_size=0.3, random_state=RANDOM_STATE
    )
    print(f"[04] Train/test split: {len(X_train)} / {len(X_test)} patients")

    # --- Regression task: predict continuous LAS ----------------------------
    reg = tune_regressor(X_train, y_train)
    y_pred = reg.predict(X_test)
    r2 = r2_score(y_test, y_pred)
    mae = mean_absolute_error(y_test, y_pred)
    rmse = float(np.sqrt(mean_squared_error(y_test, y_pred)))
    print(f"[04] Regression  R2={r2:.3f}  MAE={mae:.4f}  RMSE={rmse:.4f}")

    # --- Binary task: high vs low LAS (median split) ------------------------
    median_las = y_reg.median()
    y_bin = (y_reg >= median_las).astype(int)
    y_bin_train = (y_train >= median_las).astype(int)
    y_bin_test = (y_test >= median_las).astype(int)

    clf = tune_classifier(X_train, y_bin_train)
    prob = clf.predict_proba(X_test)[:, 1]
    pred_bin = (prob >= 0.5).astype(int)
    auc = roc_auc_score(y_bin_test, prob)
    tn, fp, fn, tp = confusion_matrix(y_bin_test, pred_bin).ravel()
    sensitivity = tp / (tp + fn) if (tp + fn) else np.nan
    specificity = tn / (tn + fp) if (tn + fp) else np.nan
    acc = accuracy_score(y_bin_test, pred_bin)
    print(f"[04] Binary (median-split LAS)  AUC={auc:.3f}  "
          f"sensitivity={sensitivity:.3f}  specificity={specificity:.3f}")

    # --- Persist metrics and model ------------------------------------------
    metrics = pd.DataFrame([{
        "task": "regression", "metric": "R2", "value": r2}, {
        "task": "regression", "metric": "MAE", "value": mae}, {
        "task": "regression", "metric": "RMSE", "value": rmse}, {
        "task": "binary_median_split", "metric": "AUC", "value": auc}, {
        "task": "binary_median_split", "metric": "sensitivity", "value": sensitivity}, {
        "task": "binary_median_split", "metric": "specificity", "value": specificity}, {
        "task": "binary_median_split", "metric": "accuracy", "value": acc},
    ])
    metrics_path = os.path.join(tables_dir, "model_metrics.csv")
    metrics.to_csv(metrics_path, index=False)

    model_path = os.path.join(models_dir, "xgb_model.joblib")
    joblib.dump({"regressor": reg, "classifier": clf,
                 "features": NUMERIC_FEATURES + CATEGORICAL_FEATURES,
                 "las_median": median_las}, model_path)

    print(f"[04] Saved metrics -> {metrics_path}")
    print(f"[04] Saved model   -> {model_path}")
    print("[04] Reminder: LAS is a proxy of lactylation-associated "
          "transcriptional activity, not a direct protein-lactylation readout.")


if __name__ == "__main__":
    main()
