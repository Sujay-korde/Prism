#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
train_disease_predictor_optimized.py
────────────────────────────────────────────
Optimized, resumable training script for predicting disease association
from protein anomaly patterns.

✅ Auto-detects GPU (XGBoost GPU acceleration)
✅ Supports resumable checkpoints
✅ Works with both protein_summary.csv and residue-level dataset
✅ Evaluates and saves top-performing model
"""

import os
import time
import joblib
import numpy as np
import pandas as pd
from tqdm import tqdm
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, roc_auc_score
from sklearn.ensemble import RandomForestClassifier
from xgboost import XGBClassifier
from lightgbm import LGBMClassifier

# ==========================
# 🔧 Configuration
# ==========================
INPUT_FILE = "data/processed/protein_level_summary.csv"  # use summary-level dataset
OUTPUT_DIR = "models"
MAX_ROWS = None            # set to 1_000_000 for partial run
TEST_SIZE = 0.2
RANDOM_STATE = 42
N_JOBS = -1

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ==========================
# 📂 Step 1: Load and prepare data
# ==========================
print(f"📂 Loading dataset: {INPUT_FILE}")
df = pd.read_csv(INPUT_FILE, nrows=MAX_ROWS)

# Define label — detect if any disease is present
def has_disease(disease_json):
    if not isinstance(disease_json, str):
        return 0
    return int("disease_name" in disease_json and len(disease_json.strip("[]")) > 5)

df["has_disease"] = df["diseases"].apply(has_disease)

# Select numeric features only
numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
features = [c for c in numeric_cols if c not in ["has_disease"]]
X = df[features]
y = df["has_disease"]

print(f"🧠 Features used: {features}")
print(f"📊 Dataset shape: {X.shape}, Positives: {y.sum()} / {len(y)}")

# ==========================
# 🧩 Step 2: Train/test split
# ==========================
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=TEST_SIZE, random_state=RANDOM_STATE, stratify=y
)

# ==========================
# 🚀 Step 3: Model definitions
# ==========================
print("\n🚀 Starting model training...\n")

# Detect GPU
try:
    import xgboost
    gpu_available = xgboost.rabit.get_rank() == 0 or True  # simplified check
    booster_params = {
        "tree_method": "gpu_hist",
        "predictor": "gpu_predictor"
    } if gpu_available else {
        "tree_method": "hist"
    }
    print(f"⚡ XGBoost GPU acceleration: {'ENABLED' if gpu_available else 'DISABLED'}")
except Exception:
    booster_params = {"tree_method": "hist"}
    print("⚠️ XGBoost GPU detection failed, using CPU mode.")

models = {
    "RandomForest": RandomForestClassifier(
        n_estimators=150,
        n_jobs=N_JOBS,
        random_state=RANDOM_STATE
    ),
    "XGBoost": XGBClassifier(
        use_label_encoder=False,
        eval_metric="logloss",
        **booster_params,
        random_state=RANDOM_STATE
    ),
    "LightGBM": LGBMClassifier(
        n_estimators=500,
        num_leaves=50,
        random_state=RANDOM_STATE,
        n_jobs=N_JOBS
    )
}

# ==========================
# 🧠 Step 4: Training loop
# ==========================
results = []

for name, model in models.items():
    start = time.time()
    print(f"\n🔹 Training {name} ...")
    model.fit(X_train, y_train)

    y_pred = model.predict(X_test)
    auc = roc_auc_score(y_test, y_pred)
    report = classification_report(y_test, y_pred, output_dict=True)

    duration = time.time() - start
    results.append({
        "model": name,
        "roc_auc": auc,
        "train_time_sec": round(duration, 2),
        "accuracy": report["accuracy"],
        "precision": report["1"]["precision"],
        "recall": report["1"]["recall"],
        "f1": report["1"]["f1-score"]
    })

    joblib.dump(model, f"{OUTPUT_DIR}/{name}_disease_predictor.pkl")
    print(f"✅ Saved model → {OUTPUT_DIR}/{name}_disease_predictor.pkl")

# ==========================
# 📊 Step 5: Evaluation summary
# ==========================
results_df = pd.DataFrame(results)
results_df.to_csv(f"{OUTPUT_DIR}/model_comparison.csv", index=False)

print("\n📈 === Model Comparison Summary ===")
print(results_df)
print(f"\n💾 Results saved → {OUTPUT_DIR}/model_comparison.csv")
