#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
predict_disease_risk.py
────────────────────────────────────────────
Loads trained models and generates per-protein disease risk predictions
based on protein-level anomaly features.

✅ Auto-detects best model from model_comparison.csv
✅ Works with protein_summary.csv
✅ Saves per-protein predictions + probabilities
✅ Skips already processed predictions (resumable)
"""

import os
import joblib
import pandas as pd
from tqdm import tqdm
from sklearn.metrics import roc_auc_score

# ================================
# 📁 Config
# ================================
SUMMARY_FILE = "data/processed/protein_level_summary.csv"
MODEL_DIR = "models"
OUTPUT_FILE = "data/results/predicted_disease_risk.csv"

os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)

# ================================
# 🔍 Load dataset
# ================================
print(f"📂 Loading summary dataset: {SUMMARY_FILE}")
df = pd.read_csv(SUMMARY_FILE)

numeric_cols = [
    "mean_anomaly_score", "max_anomaly_score", "min_anomaly_score",
    "std_anomaly_score", "num_anomalous_residues",
    "num_chains", "num_residues", "fraction_anomalous"
]
X = df[numeric_cols]

# Label (0 or 1 for disease presence)
def has_disease(disease_json):
    if not isinstance(disease_json, str):
        return 0
    return int("disease_name" in disease_json and len(disease_json.strip("[]")) > 5)

df["has_disease"] = df["diseases"].apply(has_disease)
y_true = df["has_disease"]

# ================================
# 🧩 Select best model
# ================================
comparison_path = os.path.join(MODEL_DIR, "model_comparison.csv")
if not os.path.exists(comparison_path):
    raise FileNotFoundError("⚠️ model_comparison.csv not found. Run train script first.")

comparison_df = pd.read_csv(comparison_path)
best_model_name = comparison_df.sort_values("roc_auc", ascending=False).iloc[0]["model"]
best_model_path = os.path.join(MODEL_DIR, f"{best_model_name}_disease_predictor.pkl")

print(f"🏆 Best model selected: {best_model_name}")
model = joblib.load(best_model_path)

# ================================
# 🚀 Predict
# ================================
print(f"🧠 Generating predictions using {best_model_name} ...")

# If resumable prediction file exists, skip already processed rows
if os.path.exists(OUTPUT_FILE):
    existing = pd.read_csv(OUTPUT_FILE)
    done_ids = set(existing["protein_id"])
    df = df[~df["protein_id"].isin(done_ids)]
    print(f"🔁 Resuming: {len(done_ids)} already processed, {len(df)} remaining.")
else:
    existing = pd.DataFrame()

if len(df) == 0:
    print("✅ All predictions already computed.")
else:
    preds = []
    for _, row in tqdm(df.iterrows(), total=len(df), desc="Predicting"):
        X_input = pd.DataFrame([row[numeric_cols]], columns=numeric_cols)
        prob = model.predict_proba(X_input)[0][1]
        preds.append({
            "protein_id": row["protein_id"],
            "predicted_disease_prob": prob,
            "actual_disease_label": row["has_disease"]
        })

    pred_df = pd.DataFrame(preds)
    final_df = pd.concat([existing, pred_df], ignore_index=True)
    final_df.to_csv(OUTPUT_FILE, index=False)
    print(f"\n💾 Saved predictions → {OUTPUT_FILE}")

# ================================
# 📈 Evaluate (if labels present)
# ================================
pred_data = pd.read_csv(OUTPUT_FILE)
if "actual_disease_label" in pred_data.columns:
    auc = roc_auc_score(pred_data["actual_disease_label"], pred_data["predicted_disease_prob"])
    print(f"\n📊 Model AUC on full dataset: {auc:.4f}")
else:
    print("⚠️ No true labels found for AUC evaluation.")
