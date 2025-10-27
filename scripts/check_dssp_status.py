#!/usr/bin/env python3
import glob, os
from pathlib import Path

PDB_DIR = Path("pdb_files")
OUTPUT_DIR = Path("dssp_outputs")
SAMPLE_PDB_DIR = Path("data/sample_pdb")
SAMPLE_DSSP_DIR = Path("data/sample_dssp")

def check_folder_counts(pdb_dir, dssp_dir):
    pdb_count = len(list(pdb_dir.glob("*.pdb"))) if pdb_dir.exists() else 0
    dssp_count = len(list(dssp_dir.glob("*.dssp"))) if dssp_dir.exists() else 0
    remaining = max(0, pdb_count - dssp_count)
    pct = (dssp_count / pdb_count * 100) if pdb_count else 0
    return pdb_count, dssp_count, remaining, pct

def main():
    total_pdb, total_dssp, rem, pct = check_folder_counts(PDB_DIR, OUTPUT_DIR)
    print("=== Project (local full dataset) ===")
    print(f"Total PDB files   : {total_pdb}")
    print(f"Total DSSP files  : {total_dssp}")
    print(f"Remaining         : {rem}")
    print(f"Progress          : {pct:.2f}%\n")

    # Sample dataset summary (for repo)
    sp_pdb, sp_dssp, sp_rem, sp_pct = check_folder_counts(SAMPLE_PDB_DIR, SAMPLE_DSSP_DIR)
    print("=== Repo sample dataset (data/) ===")
    print(f"Sample PDB files  : {sp_pdb}")
    print(f"Sample DSSP files : {sp_dssp}")
    print(f"Sample remaining  : {sp_rem}")
    print(f"Sample progress   : {sp_pct:.2f}%")

if __name__ == '__main__':
    main()
