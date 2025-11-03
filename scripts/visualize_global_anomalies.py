import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns

# === CONFIG ===
INPUT = "data/protein_anomalies_detected.csv"  # path to your latest dataset
OUTPUT_DIR = "results"

# === LOAD DATA ===
print(f"📂 Loading {INPUT} ...")
df = pd.read_csv(INPUT)

# Ensure numeric
df = df[pd.to_numeric(df["x"], errors="coerce").notna()]
df = df[pd.to_numeric(df["y"], errors="coerce").notna()]
df = df[pd.to_numeric(df["z"], errors="coerce").notna()]

# === ANOMALY SUMMARY ===
print("\n📊 === Global Anomaly Summary ===")
summary = (
    df.groupby("protein_id")["anomaly_label"]
    .value_counts()
    .unstack(fill_value=0)
    .rename(columns={0: "Normal", 1: "Anomalous"})
)

summary["Total"] = summary["Normal"] + summary["Anomalous"]
summary["Anomaly %"] = (summary["Anomalous"] / summary["Total"]) * 100
print(summary.sort_values("Anomaly %", ascending=False).head(10))

# === BARPLOT ===
plt.figure(figsize=(12, 6))
top10 = summary.sort_values("Anomaly %", ascending=False).head(10)
sns.barplot(x=top10.index, y=top10["Anomaly %"], palette="coolwarm")
plt.title("Top 10 Proteins by Anomaly Percentage")
plt.ylabel("Anomaly Percentage (%)")
plt.xlabel("Protein ID")
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/top_anomalous_proteins.png")
print(f"📈 Saved barplot → {OUTPUT_DIR}/top_anomalous_proteins.png")

# === GLOBAL 3D SCATTER ===
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection="3d")
subset = df.sample(n=min(10000, len(df)), random_state=42)  # sample for performance
colors = subset["anomaly_label"].map({0: "blue", 1: "red"})

ax.scatter(subset["x"], subset["y"], subset["z"], c=colors, s=4, alpha=0.6)
ax.set_title("Global Protein Structure Space — Normal (blue) vs Anomalous (red)")
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/global_anomalies_3d.png")
print(f"🌍 Saved global 3D visualization → {OUTPUT_DIR}/global_anomalies_3d.png")

print("\n✅ Done! Check 'results/' for saved visualizations.")
