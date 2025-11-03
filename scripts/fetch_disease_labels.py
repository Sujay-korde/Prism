#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
fetch_disease_labels.py
────────────────────────────
Fetches disease annotations for PDB proteins using the UniProt REST API.

INPUT  : protein_anomalies_detected.csv  (contains column 'protein_id')
OUTPUT : data/disease_labels/uniprot_disease_labels.csv
"""

import os
import time
import json
import requests
import pandas as pd
from tqdm import tqdm

# =============================
# Configuration
# =============================
INPUT_FILE  = "data/protein_anomalies_detected.csv"
OUTPUT_DIR  = "data/disease_labels"
OUTPUT_FILE = f"{OUTPUT_DIR}/uniprot_disease_labels.csv"

# API Endpoints
RCSB_MAPPING_API = "https://data.rcsb.org/rest/v1/core/entry/"
UNIPROT_API      = "https://rest.uniprot.org/uniprotkb/"

os.makedirs(OUTPUT_DIR, exist_ok=True)

# =============================
# Helper Functions
# =============================

def pdb_to_uniprot(pdb_id: str):
    """
    Maps a PDB ID (e.g. '1A08') to its UniProt accession (e.g. 'P69905').
    Uses the RCSB PDB REST API.
    """
    try:
        r = requests.get(RCSB_MAPPING_API + pdb_id.lower(), timeout=10)
        if r.status_code != 200:
            return None
        data = r.json()
        # Extract UniProt accession if present
        refs = data.get("struct", {}).get("rcsb_external_references", [])
        for ref in refs:
            if ref.get("resource_name") == "UniProt":
                return ref.get("accession")
    except Exception:
        return None
    return None


def fetch_uniprot_disease_info(uniprot_id: str):
    """
    Fetches disease annotations for a given UniProt ID.
    """
    try:
        url = f"{UNIPROT_API}{uniprot_id}.json"
        r = requests.get(url, timeout=10)
        if r.status_code != 200:
            return None

        data = r.json()
        protein_name = (
            data.get("proteinDescription", {})
            .get("recommendedName", {})
            .get("fullName", {})
            .get("value", "")
        )

        organism = data.get("organism", {}).get("scientificName", "")
        diseases = []

        for comment in data.get("comments", []):
            if comment.get("commentType") == "DISEASE":
                disease_obj = comment.get("disease", {})
                diseases.append({
                    "disease_id": disease_obj.get("diseaseId", ""),
                    "disease_name": disease_obj.get("diseaseId", "").replace("DI-", ""),
                    "description": disease_obj.get("diseaseDescription", "")
                })

        return {
            "uniprot_id": uniprot_id,
            "protein_name": protein_name,
            "organism": organism,
            "diseases": json.dumps(diseases, ensure_ascii=False)
        }
    except Exception:
        return None


# =============================
# Main Workflow
# =============================

def main():
    print(f"📂 Loading {INPUT_FILE} ...")
    df = pd.read_csv(INPUT_FILE)
    pdb_ids = sorted(df["protein_id"].dropna().unique())
    print(f"🔍 Found {len(pdb_ids)} unique PDB IDs")

    records = []
    for pdb_id in tqdm(pdb_ids, desc="Fetching UniProt mappings"):
        uni_id = pdb_to_uniprot(pdb_id)
        if not uni_id:
            print(f"⚠️  No UniProt mapping found for {pdb_id}")
            continue

        disease_info = fetch_uniprot_disease_info(uni_id)
        if not disease_info:
            print(f"⚠️  No disease info for UniProt {uni_id}")
            continue

        disease_info["pdb_id"] = pdb_id
        records.append(disease_info)
        time.sleep(0.3)  # avoid API throttling

    if not records:
        print("❌ No records fetched. Check your network or API limits.")
        return

    df_out = pd.DataFrame(records)
    df_out.to_csv(OUTPUT_FILE, index=False)
    print(f"✅ Saved {len(df_out)} mappings → {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
