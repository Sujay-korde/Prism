import requests, json, time, os
from tqdm import tqdm

UNMAPPED_FILE = "data/disease_labels/unmapped_pdbs.txt"
OUTPUT_FILE = "data/disease_labels/pdb_to_uniprot_supplement.json"

# Load previous progress if file exists
if os.path.exists(OUTPUT_FILE):
    with open(OUTPUT_FILE) as f:
        results = json.load(f)
else:
    results = {}

with open(UNMAPPED_FILE) as f:
    unmapped = [line.strip() for line in f if line.strip()]

to_process = [p for p in unmapped if p not in results]

print(f"🧩 Remaining to fetch: {len(to_process)} / {len(unmapped)}")

for pdb_id in tqdm(to_process, desc="Fetching supplementary mappings from RCSB"):
    try:
        url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            uniprot_ids = []
            for entity in data.get("rcsb_entry_container_identifiers", {}).get("polymer_entity_ids", []):
                eurl = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity}"
                r = requests.get(eurl, timeout=10)
                if r.status_code == 200:
                    edata = r.json()
                    refseqs = edata.get("rcsb_polymer_entity_container_identifiers", {}).get("reference_sequence_identifiers", [])
                    for ref in refseqs:
                        if ref.get("database_name") == "UniProt":
                            uniprot_ids.append(ref["database_accession"])
            if uniprot_ids:
                results[pdb_id] = uniprot_ids[0]
        time.sleep(0.1)  # reduce wait a bit
    except Exception as e:
        print(f"⚠️ {pdb_id}: {e}")

    # save progress every 50 entries
    if len(results) % 50 == 0:
        with open(OUTPUT_FILE, "w") as f:
            json.dump(results, f, indent=2)

# final save
with open(OUTPUT_FILE, "w") as f:
    json.dump(results, f, indent=2)

print(f"\n✅ Saved progress → {OUTPUT_FILE}")
print(f"📊 Total mapped so far: {len(results)}")
