import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os

# === Config ===
FILE_PATH = "data/protein_anomalies_detected.csv"
PROTEIN_ID = "1B08"  # 🔹 Change this to visualize another protein

# === Load dataset ===
print(f"📂 Loading {FILE_PATH} ...")
df = pd.read_csv(FILE_PATH)

if "anomaly_label" not in df.columns:
    raise ValueError("❌ 'anomaly_label' column not found. Run generate_anomalies.py first!")

subset = df[df["protein_id"] == PROTEIN_ID]

if subset.empty:
    raise ValueError(f"❌ No entries found for protein {PROTEIN_ID}. Check available IDs using df['protein_id'].unique()")

# === Visualization ===
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection="3d")

# Plot normal vs anomalous residues
normals = subset[subset["anomaly_label"] == 0]
anomalies = subset[subset["anomaly_label"] == 1]

ax.scatter(normals["x"], normals["y"], normals["z"], c="blue", s=5, label="Normal")
ax.scatter(anomalies["x"], anomalies["y"], anomalies["z"], c="red", s=10, label="Anomaly")

ax.set_title(f"3D Protein Structure with Anomalies – {PROTEIN_ID}", fontsize=14)
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")
ax.legend()
plt.tight_layout()
os.makedirs("outputs", exist_ok=True)
output_path = f"outputs/{PROTEIN_ID}_anomalies_3d.png"
plt.savefig(output_path, dpi=300, bbox_inches="tight")
print(f"✅ Saved 3D anomaly plot → {output_path}")

