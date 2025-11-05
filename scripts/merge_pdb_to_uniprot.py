#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
merge_pdb_uniprot_disease.py
───────────────────────────────
Merges PDB→UniProt mapping with UniProt→Disease labels into a final dataset.

Inputs:
  - data/disease_labels/pdb_to_uniprot.json
  - data/disease_labels/uniprot_disease_labels_final.csv

Output:
  - data/disease_labels/pdb_uniprot_disease_final.csv
"""

import json
import pandas as pd
import os

# === Paths ===
PDB_TO_UNIPROT = "data/disease_labels/pdb_to_uniprot.json"
UNIPROT_DISEASE = "data/disease_labels/uniprot_disease_labels_final.csv"
OUTPUT_FILE = "data/disease_labels/pdb_uniprot_disease_final.csv"

# === Ensure folders exist ===
os.makedirs("data/disease_labels", exist_ok=True)

print("📂 Loading files...")
with open(PDB_TO_UNIPROT, "r") as f:
    pdb_to_uniprot = json.load(f)

df_disease = pd.read_csv(UNIPROT_DISEASE)

# Normalize columns
df_disease["uniprot_id"] = df_disease["uniprot_id"].astype(str).str.strip().str.upper()

# Reverse map: UniProt → PDB list
uniprot_to_pdbs = {}
for pdb_id, uni_id in pdb_to_uniprot.items():
    if not uni_id or pd.isna(uni_id):
        continue
    uni_id = str(uni_id).strip().upper()
    uniprot_to_pdbs.setdefault(uni_id, []).append(pdb_id)

records = []

print(f"🔗 Merging {len(df_disease)} UniProt entries with {len(pdb_to_uniprot)} PDB mappings...")

for _, row in df_disease.iterrows():
    uni_id = row["uniprot_id"].strip().upper()
    pdb_list = uniprot_to_pdbs.get(uni_id, [])

    if not pdb_list:
        # Skip UniProt IDs not linked to any PDB
        continue

    for pdb_id in pdb_list:
        records.append({
            "pdb_id": pdb_id,
            "uniprot_id": uni_id,
            "protein_name": row.get("protein_name", ""),
            "organism": row.get("organism", ""),
            "diseases": row.get("diseases", "[]")
        })

# Create final dataframe
df_final = pd.DataFrame(records)

print(f"✅ Created {len(df_final)} merged PDB–UniProt–Disease records")

# Save output
df_final.to_csv(OUTPUT_FILE, index=False)
print(f"💾 Saved merged dataset → {OUTPUT_FILE}")
