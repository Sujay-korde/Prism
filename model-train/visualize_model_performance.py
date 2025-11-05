#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
visualize_model_performance.py
──────────────────────────────
Generates key evaluation visualizations for the disease predictor model.
"""

import os
import joblib
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, roc_auc_score, precision_recall_curve, confusion_matrix, ConfusionMatrixDisplay
import seaborn as sns

# ========================
# ⚙️ Config
# ========================
MODEL_PATH = "models/RandomForest_disease_predictor.pkl"
DATA_PATH = "data/risk_results/predicted_disease_risk.csv"
OUTPUT_DIR = "data/visualizations"

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ========================
# 📂 Load data
# ========================
print("📂 Loading prediction results...")
df = pd.read_csv(DATA_PATH)
y_true = df["actual_disease_label"]
y_score = df["predicted_disease_prob"]

model = joblib.load(MODEL_PATH)
print("✅ Model loaded successfully.")

# ========================
# 📈 ROC Curve
# ========================
fpr, tpr, _ = roc_curve(y_true, y_score)
auc = roc_auc_score(y_true, y_score)

plt.figure()
plt.plot(fpr, tpr, label=f"AUC = {auc:.3f}")
plt.plot([0, 1], [0, 1], "k--")
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("ROC Curve — Disease Prediction")
plt.legend()
plt.savefig(f"{OUTPUT_DIR}/roc_curve.png", dpi=300)
plt.close()

# ========================
# 🔍 Precision-Recall Curve
# ========================
precision, recall, _ = precision_recall_curve(y_true, y_score)
plt.figure()
plt.plot(recall, precision, color="blue")
plt.xlabel("Recall")
plt.ylabel("Precision")
plt.title("Precision–Recall Curve")
plt.savefig(f"{OUTPUT_DIR}/precision_recall_curve.png", dpi=300)
plt.close()

# ========================
# 🧠 Feature Importances
# ========================
try:
    importances = model.feature_importances_
    features = model.feature_names_in_
    sns.barplot(x=importances, y=features)
    plt.title("Feature Importances")
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR}/feature_importances.png", dpi=300)
    plt.close()
except Exception as e:
    print(f"⚠️  Skipped feature importances: {e}")

# ========================
# 🧾 Confusion Matrix
# ========================
cm = confusion_matrix(y_true, (y_score > 0.5).astype(int))
disp = ConfusionMatrixDisplay(confusion_matrix=cm)
disp.plot(cmap="Blues")
plt.title("Confusion Matrix")
plt.savefig(f"{OUTPUT_DIR}/confusion_matrix.png", dpi=300)
plt.close()

print(f"\n✅ Visualizations saved → {OUTPUT_DIR}")
