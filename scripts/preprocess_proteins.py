import os
import pandas as pd
from Bio.PDB import PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionException

# === Configuration ===
PDB_DIR = "pdb_files"
OUTPUT = "data/processed_protein_dataset.csv"
os.makedirs("data", exist_ok=True)

# === Safe PDB Parser ===
def safe_parse(pdb_path):
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("protein", pdb_path)
        return structure
    except PDBConstructionException as e:
        print(f"⚠️ Warning: Skipping invalid coordinates in {os.path.basename(pdb_path)} ({e})")
        return None
    except Exception as e:
        print(f"❌ Error parsing {os.path.basename(pdb_path)}: {e}")
        return None

# === Extract CA Atom Coordinates ===
def extract_features(structure):
    data = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.has_id("CA"):
                    atom = residue["CA"]
                    data.append({
                        "chain": chain.id,
                        "residue_name": residue.resname,
                        "residue_id": residue.id[1],
                        "x": atom.coord[0],
                        "y": atom.coord[1],
                        "z": atom.coord[2],
                    })
    return pd.DataFrame(data)

# === Main Preprocessing ===
merged_data = []
pdb_files = [f for f in os.listdir(PDB_DIR) if f.endswith(".pdb")]
print(f"🔍 Found {len(pdb_files)} PDB files to process.\n")

for i, pdb_file in enumerate(pdb_files, 1):
    pdb_path = os.path.join(PDB_DIR, pdb_file)
    structure = safe_parse(pdb_path)
    if structure is None:
        continue

    df = extract_features(structure)
    if df.empty:
        print(f"⚠️ No residues found in {pdb_file}, skipping.")
        continue

    df["protein_id"] = pdb_file.replace(".pdb", "")
    merged_data.append(df)

    print(f"[{i}/{len(pdb_files)}] ✅ Processed: {pdb_file} ({len(df)} residues)")

# === Combine & Save ===
if merged_data:
    final_df = pd.concat(merged_data, ignore_index=True)
    final_df.to_csv(OUTPUT, index=False)
    print(f"\n🎯 All done! Saved combined dataset to: {OUTPUT}")
else:
    print("\n⚠️ No valid PDB files were processed.")
