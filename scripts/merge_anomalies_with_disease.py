#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
merge_anomalies_with_diseases.py
────────────────────────────────────────────
Merges:
  (1) data/processed/protein_anomalies_detected.csv   ← anomaly features
  (2) data/disease_labels/pdb_to_disease.csv          ← UniProt + disease info

✅ Keeps all anomaly + coordinate + encoding columns
✅ Adds columns: uniprot_id, protein_name, organism, diseases
✅ Handles missing mappings safely
✅ Produces disease-supervised dataset ready for ML
"""

import os
import pandas as pd

# === Paths ===
ANOMALY_FILE = "data/processed/protein_anomalies_detected.csv"
DISEASE_FILE = "data/disease_labels/pdb_uniprot_disease_final.csv"
OUTPUT_FILE = "data/processed/protein_anomalies_with_diseases.csv"

# === Load data ===
print(f"📂 Loading anomaly dataset: {ANOMALY_FILE}")
anom_df = pd.read_csv(ANOMALY_FILE)

print(f"📂 Loading disease mapping: {DISEASE_FILE}")
disease_df = pd.read_csv(DISEASE_FILE)

# Normalize column names
anom_df["protein_id"] = anom_df["protein_id"].astype(str).str.upper()
disease_df["pdb_id"] = disease_df["pdb_id"].astype(str).str.upper()

# === Merge datasets on PDB ID ===
merged = pd.merge(
    anom_df,
    disease_df,
    how="left",
    left_on="protein_id",
    right_on="pdb_id"
)

# Drop duplicate PDB ID column if needed
if "pdb_id" in merged.columns:
    merged.drop(columns=["pdb_id"], inplace=True)

# === Save final merged dataset ===
os.makedirs("data/processed", exist_ok=True)
merged.to_csv(OUTPUT_FILE, index=False)

# === Summary ===
total = len(merged)
with_disease = merged["diseases"].notna().sum()
coverage = (with_disease / total) * 100

print("\n📊 === Merge Summary ===")
print(f"Total records: {total}")
print(f"With disease info: {with_disease}")
print(f"Coverage: {coverage:.2f}%")

print(f"\n💾 Saved merged dataset → {OUTPUT_FILE}")
