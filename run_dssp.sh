#!/bin/bash

# Paths
PDB_DIR="./pdb_files"
DSSP_DIR="./dssp_outputs"
mkdir -p "$DSSP_DIR"

# Number of parallel jobs (adjust to your CPU cores, e.g., 4, 6, 8)
PARALLEL_JOBS=4

# Batch size (number of files per batch)
BATCH_SIZE=200

# List all PDB files
PDB_FILES=($PDB_DIR/*.pdb)
TOTAL=${#PDB_FILES[@]}
echo "Found $TOTAL PDB files."

# Function to process one file
process_file() {
    pdb_file=$1
    filename=$(basename "$pdb_file" .pdb)
    output_file="$DSSP_DIR/$filename.dssp"

    # Skip if already processed
    if [[ -f "$output_file" ]]; then
        echo "Skipping $filename (already processed)"
        return
    fi

    # Run DSSP
    if ! dssp "$pdb_file" "$output_file" --quiet; then
        echo "⚠️ Skipped $filename: Failed DSSP assignment"
        rm -f "$output_file"
    else
        echo "Processed $filename"
    fi
}

export -f process_file
export DSSP_DIR

# Process in batches
for ((i=0; i<TOTAL; i+=BATCH_SIZE)); do
    echo "Processing batch $((i/BATCH_SIZE + 1))..."
    batch_files=("${PDB_FILES[@]:i:BATCH_SIZE}")

    # Run batch in parallel
    printf "%s\n" "${batch_files[@]}" | parallel -j $PARALLEL_JOBS process_file {}
done

echo "All batches completed."
