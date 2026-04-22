# Results Index

This document provides a complete map of all results and output files organized by analysis phase.

---

## Phase 1: Global Structural Geometry

| File | Description |
|------|-------------|
| `rmsd_scores.csv` | Global RMSD values for all 50 mutations ranked by severity |
| `rmsd_scores_all50.csv` | Complete RMSD dataset with classifications |
| `plots/rmsd_ranking.png` | Bar chart ranking all mutations by Global RMSD |
| `plots/mutation_distribution.png` | Distribution of mutation positions along p53 sequence |
| `plots/backbone_visualization.png` | 3D backbone trace comparison |

## Phase 2: Domain-Specific Analysis

| File | Description |
|------|-------------|
| `domain_rmsd.csv` | RMSD decomposed by structural domain (N-term, DBD, TET, C-term) |
| `rmsd_scores_classified.csv` | RMSD scores with severity classifications |
| `tool_comparison.csv` | Comparison against SIFT and PolyPhen-2 predictions |
| `heatmap_matrix.csv` | Per-residue displacement heatmap data |
| `dashboard.html` | Interactive web dashboard for Phase 2 results |
| `plots/domain_rmsd_chart.png` | Multi-bar chart of domain-isolated RMSD values |
| `plots/severity_classification.png` | Color-coded severity classification chart |
| `plots/tool_comparison.png` | RMSD vs SIFT vs PolyPhen-2 comparison chart |
| `plots/mutation_heatmap.png` | Heatmap of per-residue displacements |

## Phase 3: Advanced Biophysical Metrics

### CSV Data Files

| File | Description |
|------|-------------|
| `tm_scores.csv` | TM-Score fold conservation for all 50 mutations |
| `sasa_analysis.csv` | Solvent Accessible Surface Area changes |
| `secondary_structure.csv` | Secondary structure element changes (helix/sheet/coil) |
| `contact_analysis.csv` | Contact network analysis (8A cutoff) |
| `p53_dbca.csv` | DNA-Binding Competence Assessment (5-probe functional evaluation) |
| `ares_scores.csv` | Atomistic Residue Energy Scoring (thermodynamic destabilization) |
| `sve_scores.csv` | Structural Variant Evaluator (34-feature classifier) |
| `composite_scores.csv` | Composite Structural Impact Score |
| `spi_scores.csv` | Structural Pathogenicity Index |
| `local_global_impact.csv` | Local vs Global displacement ratios |
| `clustering_pca.csv` | PCA dimensionality reduction results |
| `compactness_torsion.csv` | Radius of Gyration and Ramachandran strain |

### Visualization Plots

| File | Description |
|------|-------------|
| `plots/tm_score_ranking.png` | TM-Score ranking bar chart |
| `plots/tm_vs_rmsd.png` | TM-Score vs RMSD scatter plot |
| `plots/tm_vs_rmsd_rank.png` | TM rank vs RMSD rank comparison |
| `plots/sasa_change.png` | SASA change bar chart |
| `plots/hydrophobic_exposure.png` | Hydrophobic surface exposure |
| `plots/ss_changes.png` | Secondary structure changes |
| `plots/ss_transitions.png` | Secondary structure transition types |
| `plots/contact_changes.png` | Contact network changes |
| `plots/contact_preservation.png` | Contact preservation rates |
| `plots/dbd_contact_loss.png` | DBD-specific contact losses |
| `plots/dbca_ranking.png` | DBCA score ranking |
| `plots/dbca_components.png` | DBCA 5-probe component breakdown |
| `plots/dbca_vs_rmsd.png` | DBCA vs RMSD comparison |
| `plots/dbca_zinc_detail.png` | Zinc coordination detail |
| `plots/ares_ranking.png` | ARES energy score ranking |
| `plots/ares_dde_landscape.png` | ARES contact energy landscape |
| `plots/ares_vs_rmsd.png` | ARES vs RMSD comparison |
| `plots/ares_wavefront.png` | ARES disruption wavefront propagation |
| `plots/sve_ranking.png` | SVE final pathogenicity ranking |
| `plots/sve_roc.png` | SVE ROC curve (AUC = 1.0000) |
| `plots/sve_distribution.png` | SVE score distribution |
| `plots/sve_feature_importance.png` | SVE Fisher discriminant feature weights |
| `plots/csis_ranking.png` | Composite score ranking |
| `plots/spi_ranking.png` | SPI score ranking |
| `plots/local_global_ratio.png` | Local vs Global impact ratios |
| `plots/local_vs_global_scatter.png` | Local vs Global scatter plot |
| `plots/pca_scatter.png` | PCA cluster visualization |
| `plots/pca_scree.png` | PCA variance explained |
| `plots/dendrogram.png` | Hierarchical clustering dendrogram |
| `plots/radius_of_gyration.png` | Protein compactness |

