import pandas as pd
from sklearn.ensemble import IsolationForest
import os

INPUT = "data/protein_anomalies.csv"
OUTPUT = "data/protein_anomalies_detected.csv"

# === Define expected columns ===
expected_cols = [
    "protein_id", "chain", "residue_name", "residue_id",
    "x", "y", "z", "sec_structure", "accessibility",
    "residue_name_encoded", "sec_structure_encoded",
    "anomaly_label", "anomaly_score"
]

print(f"📂 Loading dataset: {INPUT}")
df = pd.read_csv(
    INPUT,
    names=expected_cols,
    header=0,
    on_bad_lines="skip",
    engine="python"
)
print(f"✅ Loaded {len(df)} protein records with columns: {df.columns.tolist()}")

# === Keep only numeric coordinates ===
df = df[pd.to_numeric(df["x"], errors="coerce").notna()]
df = df[pd.to_numeric(df["y"], errors="coerce").notna()]
df = df[pd.to_numeric(df["z"], errors="coerce").notna()]
df[["x", "y", "z"]] = df[["x", "y", "z"]].astype(float)

print(f"📉 Cleaned numeric coordinates: {len(df)} valid rows")

# === Fit Isolation Forest for anomaly detection ===
features = ["x", "y", "z"]
model = IsolationForest(contamination=0.05, random_state=42)

print("🧠 Training IsolationForest model (this may take a few minutes)...")
model.fit(df[features])

# === Compute anomaly scores and labels ===
df["anomaly_score"] = model.decision_function(df[features])
df["anomaly_label"] = model.predict(df[features])
df["anomaly_label"] = df["anomaly_label"].map({1: 0, -1: 1})  # 1=normal, -1=anomaly

num_anomalies = df["anomaly_label"].sum()
print(f"🚨 Detected {num_anomalies} potential anomalies out of {len(df)} residues")

# === Save output ===
os.makedirs(os.path.dirname(OUTPUT), exist_ok=True)
df.to_csv(OUTPUT, index=False)
print(f"💾 Saved processed results → {OUTPUT}")
  