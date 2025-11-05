#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
🌈 Aesthetic 3D Protein Visualizer with Anomaly Highlighting
────────────────────────────────────────────────────────────
Displays protein structures, highlights anomaly sites, and shows linked diseases.
"""

import streamlit as st
import streamlit.components.v1 as components
import pandas as pd
import joblib
import py3Dmol
import json
import os

# ===============================
# 🔧 File paths
# ===============================
SUMMARY_FILE = "data/processed/protein_level_summary.csv"
PRED_FILE = "data/risk_results/predicted_disease_risk.csv"
MODEL_FILE = "models/RandomForest_disease_predictor.pkl"

# ===============================
# 🎨 Streamlit Page Settings
# ===============================
st.set_page_config(page_title="Protein Anomaly Visualizer", layout="wide")
st.title("🧬 Protein Anomaly & Disease Risk Visualizer")
st.caption("Explore predicted disease risks with interactive, high-fidelity 3D protein structures.")

# ===============================
# 📂 Load Data
# ===============================
@st.cache_data
def load_data():
    summary = pd.read_csv(SUMMARY_FILE)
    preds = pd.read_csv(PRED_FILE)
    merged = summary.merge(preds, on="protein_id", how="left")
    return merged

df = load_data()
model = joblib.load(MODEL_FILE)

# ===============================
# 🧩 Sidebar Selection
# ===============================
st.sidebar.header("🔍 Protein Selection")
protein_id = st.sidebar.selectbox("Choose a Protein (PDB ID)", sorted(df["protein_id"].unique()))

protein = df[df["protein_id"] == protein_id].iloc[0]
pred_prob = protein.get("predicted_disease_prob", 0)
actual_label = "Disease-associated" if protein.get("actual_disease_label", 0) == 1 else "No disease association"

# ===============================
# 🧠 Display Info
# ===============================
col1, col2 = st.columns([2, 1])

with col1:
    st.markdown(f"### 🧪 Predicted Disease Risk: **{pred_prob * 100:.2f}%**")
    st.markdown(f"#### 💉 True Label: :red[{actual_label}]")

    summary_data = {
        "Mean Anomaly Score": round(float(protein["mean_anomaly_score"]), 4),
        "Std Dev of Anomaly": round(float(protein["std_anomaly_score"]), 4),
        "Residues": int(protein["num_residues"]),
        "Chains": int(protein["num_chains"]),
        "Fraction Anomalous": round(float(protein["fraction_anomalous"]), 4)
    }
    st.markdown("#### 📊 Protein Summary")
    st.table(pd.DataFrame(summary_data.items(), columns=["Metric", "Value"]))

with col2:
    st.markdown("### 🧫 Disease Information")
    try:
        diseases = json.loads(protein["diseases"].replace("'", "\""))
        for d in diseases:
            st.markdown(f"**🧠 {d.get('disease_name','Unknown')}**")
            if d.get("description"):
                st.caption(d["description"])
    except Exception:
        st.info("No disease annotations found for this protein.")

# ===============================
# 🔬 3D Structure Visualization
# ===============================
st.markdown("### 🎥 Interactive 3D Protein Structure")

# Create py3Dmol viewer
view = py3Dmol.view(query=f"pdb:{protein_id}")
view.setBackgroundColor("0xF7F9FC")  # light modern background

# Style protein — mix ribbon and surface for aesthetics
view.setStyle({
    "cartoon": {"color": "spectrum", "arrows": True, "opacity": 0.9},
    "surface": {"opacity": 0.25, "color": "white"}
})

# Highlight anomaly regions if available
if "anomaly_residues" in df.columns and isinstance(protein["anomaly_residues"], str):
    try:
        anomaly_residues = json.loads(protein["anomaly_residues"])
        for res in anomaly_residues[:50]:  # highlight top 50 anomaly sites
            chain = res.get("chain", "A")
            resi = int(res.get("residue", 0))
            view.addStyle(
                {"chain": chain, "resi": resi},
                {"stick": {"color": "magenta", "radius": 0.4}}
            )
    except Exception:
        st.warning("⚠️ Unable to parse anomaly site data.")

# Auto zoom and render
view.zoomTo()
html = view._make_html()

components.html(html, height=600, width=900)

# ===============================
# ⚠️ Anomaly Alert
# ===============================
frac = protein.get("fraction_anomalous", 0)
if frac > 0.15:
    st.warning(f"⚠️ High anomaly fraction ({frac:.2f}) — possible structural instability detected!")

# ===============================
# 👣 Footer
# ===============================
st.markdown("---")
st.caption("Developed by Prism Research Team • AI-driven Structural Bioinformatics Platform.")
