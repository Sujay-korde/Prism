import pandas as pd

# ======== CONFIG =========
INPUT_FILE = "data/protein_anomalies.csv"   # Path to your latest CSV
# =========================

# Load data
df = pd.read_csv(INPUT_FILE, on_bad_lines='warn', engine='python')


print("\n📊 === Anomaly Summary ===")
print(f"Total records: {len(df)}\n")

# Check label distribution
if "anomaly_label" in df.columns:
    print("🔍 Anomaly Label Counts:")
    print(df["anomaly_label"].value_counts(), "\n")
else:
    print("⚠️ No 'anomaly_label' column found.")

# Describe anomaly scores
if "anomaly_score" in df.columns:
    print("📈 Anomaly Score Statistics:")
    print(df["anomaly_score"].describe(), "\n")
else:
    print("⚠️ No 'anomaly_score' column found.")

# Quick look at extreme anomalies (lowest scores)
if "anomaly_score" in df.columns:
    print("🚨 Top 5 Most Anomalous Residues:")
    print(df.nsmallest(5, "anomaly_score")[[
        "protein_id", "chain", "residue_name", "residue_id", "anomaly_score"
    ]])
