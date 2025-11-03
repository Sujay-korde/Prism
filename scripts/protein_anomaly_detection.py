# ============================================================
# 🧬 PRISM: Protein Structure Anomaly Detection Pipeline
# ============================================================

import pandas as pd
import numpy as np
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.ensemble import IsolationForest
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os

# === CONFIGURATION ===
INPUT_FILE = "data/processed_protein_dataset_with_dssp.csv"
OUTPUT_FILE = "data/protein_anomalies.csv"

# === STEP 1: Load & Validate Data ===
print("🔍 Loading dataset...")
if not os.path.exists(INPUT_FILE):
    raise FileNotFoundError(f"Input file not found: {INPUT_FILE}")

df = pd.read_csv(INPUT_FILE)
print(f"✅ Loaded {len(df)} residues across {df['protein_id'].nunique()} proteins.")

# Basic checks
print("\n📊 Dataset Summary:")
print(df.info())
print(df.head(5))

# Handle missing DSSP data
df["sec_structure"] = df["sec_structure"].fillna("NA")
df["accessibility"] = df["accessibility"].fillna(df["accessibility"].median())

# === STEP 2: Encode and Normalize Features ===
print("\n⚙️ Encoding and scaling features...")

# Encode categorical features
le_residue = LabelEncoder()
df["residue_name_encoded"] = le_residue.fit_transform(df["residue_name"])

le_ss = LabelEncoder()
df["sec_structure_encoded"] = le_ss.fit_transform(df["sec_structure"])

# Scale coordinates and accessibility
scaler = StandardScaler()
df[["x", "y", "z", "accessibility"]] = scaler.fit_transform(df[["x", "y", "z", "accessibility"]])

# Prepare feature set
features = df[["x", "y", "z", "residue_name_encoded", "sec_structure_encoded", "accessibility"]]

# === STEP 3: Train Isolation Forest ===
print("\n🧠 Training Isolation Forest for anomaly detection...")

model = IsolationForest(
    n_estimators=200,
    contamination=0.05,   # 5% of data treated as anomalies
    random_state=42,
    n_jobs=-1
)
df["anomaly_label"] = model.fit_predict(features)
df["anomaly_score"] = model.decision_function(features)

# Label interpretation: -1 = anomaly, 1 = normal
num_anomalies = (df["anomaly_label"] == -1).sum()
print(f"✅ Model trained. Detected {num_anomalies} potential anomalies out of {len(df)} residues.")

# === STEP 4: Save Results ===
df.to_csv(OUTPUT_FILE, index=False)
print(f"💾 Results saved to: {OUTPUT_FILE}")

# === STEP 5: Visualize Example Protein ===
print("\n🎨 Generating a 3D scatter visualization...")

# Pick one protein to visualize
sample_id = df["protein_id"].unique()[0]
sample = df[df["protein_id"] == sample_id]

fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111, projection="3d")

normal = sample[sample["anomaly_label"] == 1]
anomaly = sample[sample["anomaly_label"] == -1]

ax.scatter(normal["x"], normal["y"], normal["z"], c="blue", s=20, label="Normal")
ax.scatter(anomaly["x"], anomaly["y"], anomaly["z"], c="red", s=35, label="Anomaly")

ax.set_title(f"Protein: {sample_id} (Anomalies in Red)")
ax.legend()
plt.show()

print("\n✅ Visualization complete. You can now explore anomaly positions visually!")
