# TP53-SVE: Structural Variation Engine
### AI-Powered Structural Pathogenicity Analysis of TP53 Mutations

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![AlphaFold3](https://img.shields.io/badge/Powered%20By-AlphaFold3-blue)](https://alphafoldserver.com/)

**TP53-SVE** (Structural Variation Engine) is a novel bioinformatics pipeline designed to quantify the structural impact of cancer-associated missense mutations in the TP53 protein using **AlphaFold3** predictions and interpretable biophysical descriptors.

---

## 🚀 The Core Novelty

While most pathogenicity predictors rely on sequence conservation or global structural metrics (like RMSD), this project demonstrates that:
1. **Global Metrics Fail**: Standard metrics like backbone RMSD and graph centrality often fail to discriminate between pathogenic and benign variants in AlphaFold models (p=0.96).
2. **Micro-Environmental Success**: Localized biophysical perturbations—such as **SASA exposure**, **contact network rewiring**, and **Zinc-binding pocket displacement**—contain strong signals for pathogenicity.
3. **Interpretable ML**: Uses high-performance but transparent models (Fisher LDA & Kolmogorov-Arnold Networks) to provide mechanistic explanations for every classification.

---

## 🏗️ Repository Architecture

- **`src/`**: Core logic divided into three phases:
    - `phase1_rmsd/`: Initial structural deviation analysis.
    - `phase2_features/`: Multi-dimensional feature extraction (SASA, H-bonds, SS-changes).
    - `phase3_kan_ml/`: Advanced scoring models and the novel KAN classifier.
- **`dashboard/`**: Source code for the interactive React/Streamlit web visualizer.
- **`docs/manuscripts/`**: Full-length research papers and technical reports.
- **`data/samples/`**: 3D structural models (.cif) for various TP53 mutants.
- **`output/`**: Curated visualizations, severeity rankings, and SHAP explainability plots.

---

## 🛠️ Getting Started

### Prerequisites
- Python 3.10+
- [ChimeraX](https://www.rbvi.ucsf.edu/chimerax/) (for 3D visualization scripts)

### Installation
```bash
git clone https://github.com/your-username/tp53-sve.git
cd tp53-sve
pip install -r requirements.txt
```

### Running the Dashboard
```bash
streamlit run antigravity_webapp.py
```

## 📦 External Data Hosting (1.3GB)

To keep this repository lightweight, the full AlphaFold3 structural database (128 variants, ~1.3GB) is hosted externally. 

1. **Download**: [Placeholder: Link to Google Drive / Zenodo]
2. **Setup**: Extracted structural files should be placed in `data/structures/`.
3. **Reproducibility**: Once the data is present, you can rerun the full Phase 1-3 pipeline using `src/utils/run_full_pipeline.py`.

---

## 📊 Key Results

- **Classification Performance**: Achieved **100% LOOCV Accuracy** on confirmed pathogenic/benign subsets.
- **Scaling**: Successfully analyzed **128 clinical variants** from COSMIC & ClinVar.
- **Explainability**: Identified Zinc-coordination disruption and Hydrophobic core exposure as the primary drivers of p53 structural collapse.

---

## 📄 Documentation

For detailed methodology, see the final research manuscript:
[TP53_SVE_ULTIMATE_RESEARCH_PAPER.md](docs/manuscripts/TP53_SVE_ULTIMATE_RESEARCH_PAPER.md)

## 🎓 Citation

If you use this software or the Structural Pathogenicity Index (SPI) in your research, please cite:

```bibtex
@article{VortexQuasarX2026,
  title={TP53-SVE: Interpretable Pathogenicity Prediction via Structural Variation Ensembles},
  author={VortexQuasarX and Gemini AI},
  journal={Bioinformatics Exploration},
  year={2026},
  url={https://github.com/VortexQuasarX/tp53-sve}
}
```

---

## 📜 License
Distributed under the MIT License. See `LICENSE` for more information.
