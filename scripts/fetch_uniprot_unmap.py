import requests
import json
import time
from tqdm import tqdm
import pandas as pd
import os

# Input mapping
MAPPING_FILE = "data/disease_labels/pdb_to_uniprot_supplement.json"
OUTPUT_FILE = "data/disease_labels/uniprot_disease_labels_final.csv"

with open(MAPPING_FILE, "r") as f:
    pdb_to_uniprot = json.load(f)

# Unique UniProt IDs
uniprot_ids = list({v for v in pdb_to_uniprot.values() if v})

print(f"🔍 Found {len(uniprot_ids)} unique UniProt IDs to fetch.")

records = []
if os.path.exists(OUTPUT_FILE):
    records = pd.read_csv(OUTPUT_FILE).to_dict(orient="records")

fetched_ids = {r["uniprot_id"] for r in records}

for uid in tqdm(uniprot_ids):
    if uid in fetched_ids:
        continue

    url = f"https://rest.uniprot.org/uniprotkb/{uid}.json"
    try:
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            protein_name = data.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", "")
            organism = data.get("organism", {}).get("scientificName", "")
            diseases = []

            for c in data.get("comments", []):
                if c.get("commentType") == "DISEASE":
                    disease = c.get("disease", {})
                    diseases.append({
                        "disease_id": disease.get("diseaseId", ""),
                        "disease_name": disease.get("acronym", ""),
                        "description": c.get("text", [{}])[0].get("value", "")
                    })

            records.append({
                "uniprot_id": uid,
                "protein_name": protein_name,
                "organism": organism,
                "diseases": json.dumps(diseases)
            })

            # Save progress every 20 requests
            if len(records) % 20 == 0:
                pd.DataFrame(records).to_csv(OUTPUT_FILE, index=False)
        else:
            print(f"⚠️ Failed for {uid} ({response.status_code})")

    except Exception as e:
        print(f"❌ Error fetching {uid}: {e}")

    time.sleep(0.5)  # avoid hitting API rate limits

# Final save
pd.DataFrame(records).to_csv(OUTPUT_FILE, index=False)
print(f"\n✅ Saved disease mappings → {OUTPUT_FILE}")
