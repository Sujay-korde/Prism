#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
aggregate_protein_level.py
──────────────────────────────────────────────
Creates a protein-level summary dataset for supervised learning.

Each protein gets one record summarizing:
- Mean, max, min, std of anomaly scores
- Fraction of anomalous residues
- Chain and residue counts
- Disease info (if available)
✅ Designed for Stage 1 learning (global disease prediction)
✅ Keeps residue-level data intact for later Stage 2 & 3 modeling
"""

import pandas as pd
import os
from tqdm import tqdm

# === Paths ===
INPUT_FILE = "data/processed/protein_anomalies_with_diseases.csv"
OUTPUT_FILE = "data/processed/protein_level_summary.csv"

print(f"📂 Loading {INPUT_FILE} ...")
df = pd.read_csv(INPUT_FILE)

# --- Sanity checks ---
required_cols = ["protein_id", "chain", "residue_id", "anomaly_score", "anomaly_label"]
missing = [c for c in required_cols if c not in df.columns]
if missing:
    raise ValueError(f"Missing columns: {missing}")

# --- Numeric conversions ---
df["anomaly_score"] = pd.to_numeric(df["anomaly_score"], errors="coerce")
df = df.dropna(subset=["anomaly_score"])

# --- Define aggregation logic ---
agg_dict = {
    "anomaly_score": ["mean", "max", "min", "std"],
    "anomaly_label": "sum",
    "chain": pd.Series.nunique,
    "residue_id": "count",
}

# --- Aggregate per protein ---
print("🔄 Aggregating protein-level statistics ...")
summary = (
    df.groupby("protein_id")
      .agg(agg_dict)
      .reset_index()
)

# --- Rename columns ---
summary.columns = [
    "protein_id",
    "mean_anomaly_score",
    "max_anomaly_score",
    "min_anomaly_score",
    "std_anomaly_score",
    "num_anomalous_residues",
    "num_chains",
    "num_residues"
]

# --- Compute derived features ---
summary["fraction_anomalous"] = summary["num_anomalous_residues"] / summary["num_residues"]

# --- Merge with disease info (if present) ---
disease_cols = ["protein_name", "organism", "diseases"]
if all(col in df.columns for col in disease_cols):
    disease_info = df.groupby("protein_id")[disease_cols].first().reset_index()
    summary = summary.merge(disease_info, on="protein_id", how="left")

# --- Save output ---
os.makedirs("data/processed", exist_ok=True)
summary.to_csv(OUTPUT_FILE, index=False)

print(f"\n✅ Saved aggregated protein-level summary → {OUTPUT_FILE}")
print(f"🧩 Total proteins summarized: {len(summary)}")
print(f"📊 Columns: {list(summary.columns)}")
