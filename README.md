# 🧬 CRISPR-to-Structure Pipeline

## Overview
This application is a high-throughput screening tool designed to bridge the gap between CRISPR-Cas9 genome editing and protein structural biology. It automates the process of gRNA design, protein translation simulation, and 3D structural impact assessment using the ESMFold API.

This tool is specifically designed for R&D professionals in biotechnology and drug discovery who need to evaluate how specific genetic edits might alter protein folding and function.

## Key Features
- **Intelligent gRNA Design**: Automatically searches for Sense (+) and Antisense (-) PAM sites (NGG/CCN).
- **In-silico Editing & Translation**: Simulates knock-in sequences and predicts the resulting amino acid sequence.
- **Structural Impact Prediction**:
    - Fetches 3D structures via **ESMFold API**.
    - Calculates **RMSD (Root Mean Square Deviation)** to quantify structural changes.
    - Provides interactive 3D visualization using **3Dmol.js**.
- **Off-target Risk Assessment**: Estimates potential off-target risks based on sequence motifs.
- **Professional Reporting**: Generates downloadable PDF reports including gRNA metrics, protein sequences, and structural analysis results.

## Tech Stack
- **Language**: Python 3.11+
- **Frontend**: Streamlit
- **Bioinformatics**: Biopython (Seq, PDB)
- **3D Visualization**: 3Dmol.js (integrated via custom HTML/JS)
- **Reporting**: FPDF
- **Structure Prediction**: ESMFold API (Meta AI)

## Installation & Usage
1. Clone the repository:
   ```bash
   git clone [https://github.com/YOUR_USERNAME/CRISPR-Structure-App.git](https://github.com/YOUR_USERNAME/CRISPR-Structure-App.git)
   cd CRISPR-Structure-App
Create and activate a virtual environment:

Bash
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate
Install dependencies:

Bash
pip install -r requirements.txt
Run the application:

Bash
streamlit run app.py
Author
Tsubaki

Ph.D. in Agriculture

R&D Researcher in Antibody Drug Development