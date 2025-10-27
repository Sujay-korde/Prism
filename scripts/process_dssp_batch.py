import os
import glob
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
import time

# === Configuration ===
PDB_DIR = "pdb_files"          # Folder with PDB files
OUTPUT_DIR = "dssp_outputs"    # Folder for DSSP results
os.makedirs(OUTPUT_DIR, exist_ok=True)

BATCH_SIZE = 5                 # Files per batch
MAX_WORKERS = 1                # Parallel processes
TIMEOUT = 60                   # Max seconds per file

# === DSSP runner ===
def run_dssp(pdb_file):
    pdb_name = os.path.basename(pdb_file)
    output_file = os.path.join(OUTPUT_DIR, pdb_name.replace(".pdb", ".dssp"))

    # Skip existing output
    if os.path.exists(output_file):
        return pdb_file, "skipped"

    if os.path.getsize(pdb_file) == 0:
        return pdb_file, "empty_file"

    try:
        result = subprocess.run(
            ["dssp", pdb_file.replace("\\", "/"), output_file.replace("\\", "/")],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=TIMEOUT
        )

        if result.stderr.strip():
            with open("dssp_warnings.log", "a") as f:
                f.write(f"File: {pdb_file}\n{result.stderr}\n\n")

        return pdb_file, "processed"

    except subprocess.TimeoutExpired:
        with open("dssp_errors.log", "a") as f:
            f.write(f"File: {pdb_file} Error: DSSP timed out after {TIMEOUT}s\n")
        return pdb_file, "timeout"

    except Exception as e:
        with open("dssp_errors.log", "a") as f:
            f.write(f"File: {pdb_file} Error: {e}\n")
        return pdb_file, "error"


# === File discovery ===
all_pdb_files = sorted(glob.glob(os.path.join(PDB_DIR, "*.pdb")))
processed_files = set(os.path.basename(f).replace(".dssp", ".pdb")
                      for f in glob.glob(os.path.join(OUTPUT_DIR, "*.dssp")))
remaining_files = [f for f in all_pdb_files if os.path.basename(f) not in processed_files]

print(f"Found {len(all_pdb_files)} PDB files.")
print(f"{len(remaining_files)} files remaining to process.")

# === Batch processing ===
processed_count = skipped_count = error_count = 0

for batch_start in range(0, len(remaining_files), BATCH_SIZE):
    batch_files = remaining_files[batch_start:batch_start + BATCH_SIZE]
    batch_num = batch_start // BATCH_SIZE + 1
    print(f"\n🚀 Processing batch {batch_num}/{(len(remaining_files)//BATCH_SIZE)+1} ({len(batch_files)} files)...")

    batch_results = []
    start_time = time.time()

    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = [executor.submit(run_dssp, pdb_file) for pdb_file in batch_files]
        for future in tqdm(as_completed(futures), total=len(futures), desc=f"Batch {batch_num} progress"):
            pdb_file, status = future.result()
            batch_results.append((pdb_file, status))

    for _, status in batch_results:
        if status == "processed":
            processed_count += 1
        elif status in ["skipped", "empty_file", "timeout"]:
            skipped_count += 1
        else:
            error_count += 1

    print(f"\n✅ Processed: {processed_count} | ⚠️ Skipped: {skipped_count} | ❌ Errors: {error_count}")
    print(f"Remaining: {len(remaining_files) - (batch_start + len(batch_files))}")
    print(f"⏱ Batch {batch_num} completed in {time.time() - start_time:.2f}s\n", flush=True)

print("\n🎯 All batches completed!")
print(f"DSSP outputs: {OUTPUT_DIR}")
print("Warnings logged in: dssp_warnings.log")
print("Errors logged in: dssp_errors.log")
