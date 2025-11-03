import os
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
from Bio.PDB import PDBParser, DSSP
from Bio.PDB.PDBExceptions import PDBConstructionException

# === CONFIGURATION ===
PDB_DIR = "pdb_files"
DSSP_DIR = "dssp_outputs"
OUTPUT_DIR = "data"
FINAL_OUTPUT = os.path.join(OUTPUT_DIR, "processed_protein_dataset_with_dssp.csv")
TEMP_FILE = os.path.join(OUTPUT_DIR, "processed_partial.csv")
os.makedirs(OUTPUT_DIR, exist_ok=True)

MAX_WORKERS = 6  # adjust based on your CPU
BATCH_SAVE_INTERVAL = 50  # save every 50 proteins processed

# === SAFE PARSER ===
def safe_parse(pdb_path):
    parser = PDBParser(QUIET=True)
    try:
        return parser.get_structure("protein", pdb_path)
    except PDBConstructionException as e:
        return None
    except Exception:
        return None

# === PROCESS ONE PDB ===
def process_single(pdb_file):
    pdb_path = os.path.join(PDB_DIR, pdb_file)
    dssp_path = os.path.join(DSSP_DIR, pdb_file.replace(".pdb", ".dssp"))

    structure = safe_parse(pdb_path)
    if structure is None:
        return None

    # DSSP parsing (if available)
    dssp_data = {}
    if os.path.exists(dssp_path):
        try:
            model = next(structure.get_models())
            dssp = DSSP(model, pdb_path, dssp=dssp_path)
            for key, value in dssp.property_dict.items():
                chain, res_id = key
                resnum = res_id[1]
                dssp_data[(chain, resnum)] = {
                    "sec_structure": value["SS"],
                    "accessibility": value["ASA"]
                }
        except Exception:
            pass

    # Extract coordinates
    data = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.has_id("CA"):
                    atom = residue["CA"]
                    resnum = residue.id[1]
                    dssp_info = dssp_data.get((chain.id, resnum), {})
                    data.append({
                        "protein_id": pdb_file.replace(".pdb", ""),
                        "chain": chain.id,
                        "residue_name": residue.resname,
                        "residue_id": resnum,
                        "x": atom.coord[0],
                        "y": atom.coord[1],
                        "z": atom.coord[2],
                        "sec_structure": dssp_info.get("sec_structure", "NA"),
                        "accessibility": dssp_info.get("accessibility", None)
                    })
    return pd.DataFrame(data)

# === MAIN FUNCTION ===
def main():
    all_pdbs = [f for f in os.listdir(PDB_DIR) if f.endswith(".pdb")]
    processed = set()

    # Resume support
    if os.path.exists(TEMP_FILE):
        existing = pd.read_csv(TEMP_FILE)
        processed = set(existing["protein_id"].unique())
        print(f"🔁 Resuming: {len(processed)} already processed, skipping...")

    remaining_pdbs = [f for f in all_pdbs if f.replace(".pdb", "") not in processed]
    print(f"🧩 Total to process now: {len(remaining_pdbs)} PDBs\n")

    partial_results = []

    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {executor.submit(process_single, pdb): pdb for pdb in remaining_pdbs}
        for i, future in enumerate(tqdm(as_completed(futures), total=len(futures), desc="Processing Proteins")):
            df = future.result()
            if df is not None and not df.empty:
                partial_results.append(df)

            # Auto-save every few proteins
            if len(partial_results) >= BATCH_SAVE_INTERVAL:
                pd.concat(partial_results).to_csv(TEMP_FILE, mode='a', index=False, header=not os.path.exists(TEMP_FILE))
                partial_results.clear()

    # Final save
    if partial_results:
        pd.concat(partial_results).to_csv(TEMP_FILE, mode='a', index=False, header=not os.path.exists(TEMP_FILE))

    print("✅ All batches completed. Merging results...")

    # Merge all partial data
    df_final = pd.read_csv(TEMP_FILE)
    df_final.drop_duplicates(subset=["protein_id", "residue_id"], inplace=True)
    df_final.to_csv(FINAL_OUTPUT, index=False)
    print(f"🎯 Final dataset saved to: {FINAL_OUTPUT}")

if __name__ == "__main__":
    main()
