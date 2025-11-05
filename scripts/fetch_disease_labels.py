#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
fetch_disease_labels_resumable.py
────────────────────────────────────────────
Resumable version that safely fetches disease info from UniProt using PDB IDs.

✅ Auto-saves progress (cache + partial results)
✅ Skips already processed PDBs
✅ Recovers automatically after crashes or network errors
"""

import os
import json
import time
import pandas as pd
import requests
from tqdm import tqdm

# === Paths ===
INPUT_FILE = "data/processed/protein_anomalies_detected.csv"
CACHE_FILE = "data/disease_labels/uniprot_cache.json"
PARTIAL_FILE = "data/disease_labels/fetched_records.csv"
OUTPUT_FILE = "data/disease_labels/uniprot_disease_labels_final.csv"

# === Create folders ===
os.makedirs("data/disease_labels", exist_ok=True)

GRAPHQL_URL = "https://data.rcsb.org/graphql"
UNIPROT_API = "https://rest.uniprot.org/uniprotkb/"

# ================================
# 🔹 Helper Functions
# ================================

def pdb_to_uniprot(pdb_id):
    """Try RCSB GraphQL API to map PDB → UniProt."""
    query = {
        "query": f"""
        {{
          entry(entry_id: "{pdb_id.upper()}") {{
            polymer_entities {{
              rcsb_polymer_entity_container_identifiers {{
                reference_sequence_identifiers {{
                  database_name
                  database_accession
                }}
              }}
            }}
          }}
        }}
        """
    }

    try:
        r = requests.post(GRAPHQL_URL, json=query, timeout=10)
        if r.status_code == 200:
            result = r.json()
            entities = result.get("data", {}).get("entry", {}).get("polymer_entities", [])
            for e in entities:
                refs = e.get("rcsb_polymer_entity_container_identifiers", {}).get("reference_sequence_identifiers", [])
                for ref in refs:
                    if ref.get("database_name") == "UniProt":
                        return ref.get("database_accession")
    except Exception:
        pass
    return None


def fetch_uniprot_disease(uniprot_id):
    """Fetch disease information for a UniProt ID."""
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
                dis = comment.get("disease", {})
                diseases.append({
                    "disease_id": dis.get("diseaseId", ""),
                    "disease_name": dis.get("diseaseId", "").replace("DI-", ""),
                    "description": dis.get("diseaseDescription", "")
                })

        return {
            "uniprot_id": uniprot_id,
            "protein_name": protein_name,
            "organism": organism,
            "diseases": json.dumps(diseases, ensure_ascii=False)
        }

    except Exception:
        return None


# ================================
# 🔹 Main Resumable Script
# ================================

def main():
    print(f"📂 Loading {INPUT_FILE} ...")
    df = pd.read_csv(INPUT_FILE)
    pdb_ids = sorted(df["protein_id"].dropna().unique())
    print(f"🔍 Found {len(pdb_ids)} unique PDB IDs")

    # Load existing cache
    if os.path.exists(CACHE_FILE):
        with open(CACHE_FILE, "r") as f:
            cache = json.load(f)
    else:
        cache = {}

    # Load partial results
    if os.path.exists(PARTIAL_FILE):
        done_df = pd.read_csv(PARTIAL_FILE)
        done_pdbs = set(done_df["pdb_id"])
        results = done_df.to_dict(orient="records")
        print(f"🔁 Resuming: {len(done_pdbs)} already processed")
    else:
        results = []
        done_pdbs = set()

    for pdb_id in tqdm(pdb_ids, desc="Fetching UniProt mappings"):
        if pdb_id in done_pdbs:
            continue

        # --- Step 1: Map PDB → UniProt ---
        if pdb_id in cache:
            uni_id = cache[pdb_id]
        else:
            uni_id = pdb_to_uniprot(pdb_id)
            cache[pdb_id] = uni_id
            with open(CACHE_FILE, "w") as f:
                json.dump(cache, f, indent=2)

        if not uni_id:
            print(f"⚠️  No UniProt mapping found for {pdb_id}")
            continue

        # --- Step 2: Fetch disease data ---
        disease_info = fetch_uniprot_disease(uni_id)
        if not disease_info:
            print(f"⚠️  No disease info found for UniProt {uni_id}")
            continue

        disease_info["pdb_id"] = pdb_id
        results.append(disease_info)

        # --- Step 3: Auto-save every record ---
        pd.DataFrame(results).to_csv(PARTIAL_FILE, index=False)
        time.sleep(0.25)

    # Save final merged output
    pd.DataFrame(results).to_csv(OUTPUT_FILE, index=False)
    print(f"\n✅ Done! Saved all {len(results)} mappings → {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
