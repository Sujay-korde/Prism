import os
import pandas as pd
import numpy as np
from sklearn.ensemble import IsolationForest
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from tqdm import tqdm

# === CONFIG ===
INPUT = "data/protein_anomalies_detected.csv"
OUTPUT_DIR = "results/local_anomalies"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# === LOAD DATA ===
print(f"📂 Loading dataset: {INPUT}")
df = pd.read_csv(INPUT, low_memory=False)

required_cols = ["protein_id", "x", "y", "z"]
missing = [c for c in required_cols if c not in df.columns]
if missing:
    raise ValueError(f"❌ Missing columns: {missing}")

df = df.dropna(subset=["x", "y", "z"])
print(f"✅ Loaded {len(df):,} residues across {df['protein_id'].nunique()} proteins")

# === DETECT LOCAL ANOMALIES PER PROTEIN ===
proteins = df["protein_id"].unique()

for pid in tqdm(proteins, desc="🔍 Processing proteins"):
    sub_df = df[df["protein_id"] == pid].copy()

    if len(sub_df) < 10:
        continue  # too small to model reliably

    features = sub_df[["x", "y", "z"]].to_numpy()

    # Train IsolationForest
    model = IsolationForest(contamination=0.05, random_state=42)
    model.fit(features)

    sub_df["local_anomaly_label"] = model.predict(features)
    sub_df["local_anomaly_score"] = model.decision_function(features)

    # Convert sklearn outputs: (-1 → 1 for anomaly, 1 → 0 for normal)
    sub_df["local_anomaly_label"] = sub_df["local_anomaly_label"].map({1: 0, -1: 1})

    # Save per-protein results
    out_csv = os.path.join(OUTPUT_DIR, f"{pid}_local_anomalies.csv")
    sub_df.to_csv(out_csv, index=False)

    # === Optional Visualization ===
    fig = plt.figure(figsize=(6, 5))
    ax = fig.add_subplot(111, projection="3d")

    normal = sub_df[sub_df["local_anomaly_label"] == 0]
    anomaly = sub_df[sub_df["local_anomaly_label"] == 1]

    ax.scatter(normal["x"], normal["y"], normal["z"], c="blue", s=10, alpha=0.5, label="Normal")
    ax.scatter(anomaly["x"], anomaly["y"], anomaly["z"], c="red", s=20, alpha=0.8, label="Anomaly")

    ax.set_title(f"Local Anomalies in Protein {pid}")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.legend()

    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, f"{pid}_3d.png"))
    plt.close(fig)

print(f"\n✅ Done! Saved per-protein results and 3D plots → {OUTPUT_DIR}/")
