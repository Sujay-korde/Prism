#!/usr/bin/env python3
"""
map_pdb_to_uniprot_parallel.py
------------------------------------------------
High-performance PDB → UniProt mapper with resume capability.

✅ Features:
- Parallel API calls (10 threads)
- Resumable (loads existing JSON)
- Fetches only missing/null entries
- Dual-source mapping: UniProt + RCSB
- Progress autosave & log support
"""

import os
import json
import time
import requests
import pandas as pd
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed

# ==========================
# ⚙️ CONFIGURATION
# ==========================
INPUT_FILE = "data/processed/protein_anomalies_detected.csv"
OUTPUT_JSON = "data/disease_labels/pdb_to_uniprot.json"
UNMAPPED_LOG = "data/disease_labels/unmapped_pdbs.txt"
MAX_THREADS = 10  # 🧵 Adjust for speed vs. stability
SAVE_EVERY = 100  # Save progress every 100 results

os.makedirs(os.path.dirname(OUTPUT_JSON), exist_ok=True)

# ==========================
# 🧩 API HELPERS
# ==========================

def fetch_uniprot_from_uniprot_api(pdb_id):
    """Query UniProt mapping API."""
    try:
        url = f"https://rest.uniprot.org/uniprotkb/stream?fields=accession&query=database:(type:pdb {pdb_id})"
        r = requests.get(url, timeout=10)
        if r.status_code == 200 and r.text.strip():
            return r.text.strip().split("\n")[0]
    except Exception:
        pass
    return None


def fetch_uniprot_from_rcsb(pdb_id):
    """Fallback using RCSB polymer entity API."""
    for entity_id in range(1, 5):
        try:
            url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_id}"
            r = requests.get(url, timeout=10)
            if r.status_code == 200:
                data = r.json()
                refs = data.get("rcsb_polymer_entity_container_identifiers", {}).get("reference_sequence_identifiers", [])
                for ref in refs:
                    if ref.get("database_name") == "UniProt":
                        return ref.get("database_accession")
        except Exception:
            continue
    return None


def fetch_uniprot(pdb_id):
    """Try both APIs with retry."""
    for attempt in range(2):
        uid = fetch_uniprot_from_uniprot_api(pdb_id)
        if uid:
            return uid
        uid = fetch_uniprot_from_rcsb(pdb_id)
        if uid:
            return uid
        time.sleep(0.3)
    return None


# ==========================
# 🚀 MAIN EXECUTION
# ==========================
if __name__ == "__main__":
    print(f"📂 Loading dataset: {INPUT_FILE}")
    df = pd.read_csv(INPUT_FILE, usecols=["protein_id"])
    pdb_ids = sorted(set(df["protein_id"].dropna().astype(str)))
    print(f"🔍 Found {len(pdb_ids):,} unique PDB IDs")

    # Load existing progress
    if os.path.exists(OUTPUT_JSON):
        with open(OUTPUT_JSON, "r") as f:
            pdb_to_uniprot = json.load(f)
        print(f"♻️ Resuming from {len(pdb_to_uniprot)} mappings")
    else:
        pdb_to_uniprot = {}

    # Determine unmapped IDs (null or missing)
    unmapped_ids = [p for p in pdb_ids if p not in pdb_to_uniprot or not pdb_to_uniprot[p]]
    print(f"🧩 Remaining unmapped IDs: {len(unmapped_ids):,}")

    if not unmapped_ids:
        print("✅ All PDBs are already mapped!")
        exit()

    # Fetch in parallel
    results = {}
    with ThreadPoolExecutor(max_workers=MAX_THREADS) as executor:
        futures = {executor.submit(fetch_uniprot, pdb_id): pdb_id for pdb_id in unmapped_ids}

        for i, future in enumerate(tqdm(as_completed(futures), total=len(futures), desc="Fetching mappings")):
            pdb_id = futures[future]
            try:
                uniprot_id = future.result()
                results[pdb_id] = uniprot_id
                pdb_to_uniprot[pdb_id] = uniprot_id
                if not uniprot_id:
                    with open(UNMAPPED_LOG, "a") as log:
                        log.write(f"{pdb_id}\n")
            except Exception as e:
                print(f"⚠️ Error fetching {pdb_id}: {e}")

            # Auto-save every N results
            if (i + 1) % SAVE_EVERY == 0:
                json.dump(pdb_to_uniprot, open(OUTPUT_JSON, "w"), indent=2)

    # Final save
    json.dump(pdb_to_uniprot, open(OUTPUT_JSON, "w"), indent=2)

    # Summary
    mapped = sum(1 for v in pdb_to_uniprot.values() if v)
    missing = sum(1 for v in pdb_to_uniprot.values() if not v)
    print("\n📊 === Final Mapping Summary ===")
    print(f"✅ Mapped: {mapped:,}")
    print(f"❌ Missing: {missing:,}")
    print(f"📈 Coverage: {mapped / (mapped + missing):.2%}")
    print(f"💾 Saved mappings → {OUTPUT_JSON}")
    print(f"🗒️  Unmapped list → {UNMAPPED_LOG}")
