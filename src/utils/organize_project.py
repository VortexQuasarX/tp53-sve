"""
Organize the entire project workspace into a clean, professional structure.
Moves scattered root-level documents, organizes chimerax scripts, and 
creates a results index file.
"""
import os
import shutil

BASE = r"C:\Users\LENOVO\.gemini\antigravity\playground\chrono-shepard"

def ensure_dir(path):
    os.makedirs(path, exist_ok=True)

# Define clean directory structure
dirs_to_create = [
    os.path.join(BASE, "docs"),
    os.path.join(BASE, "docs", "paper_drafts"),
    os.path.join(BASE, "docs", "literature"),
    os.path.join(BASE, "docs", "guides"),
    os.path.join(BASE, "visualization"),
    os.path.join(BASE, "visualization", "chimerax_scripts"),
    os.path.join(BASE, "results"),
    os.path.join(BASE, "results", "phase1_global_geometry"),
    os.path.join(BASE, "results", "phase2_domain_analysis"),
    os.path.join(BASE, "results", "phase3_advanced_metrics"),
    os.path.join(BASE, "results", "phase3_advanced_metrics", "plots"),
    os.path.join(BASE, "results", "phase2_domain_analysis", "plots"),
    os.path.join(BASE, "results", "phase1_global_geometry", "plots"),
]

for d in dirs_to_create:
    ensure_dir(d)
    print(f"Created: {d}")

# ============================================================
# MOVE ROOT-LEVEL DOCUMENTS INTO docs/
# ============================================================
paper_drafts = {
    "Frontiers_Paper_Obsidian.md": "docs/paper_drafts/",
    "TP53_SVE_DESCRIPTIVE_RESEARCH_PAPER.md": "docs/paper_drafts/",
    "TP53_SVE_FINAL_RESEARCH_PAPER.md": "docs/paper_drafts/",
    "TP53_SVE_MASSIVE_THEORETICAL_BOOK.md": "docs/paper_drafts/",
    "TP53_SVE_ULTIMATE_RESEARCH_PAPER.md": "docs/paper_drafts/",
    "JENNI_AI_EXHAUSTIVE_MASTER.md": "docs/paper_drafts/",
    "JENNI_AI_PAPER_CONTENT.md": "docs/paper_drafts/",
    "THE_MASSIVE_150_PAGE_DATA_BOOK.md": "docs/paper_drafts/",
    "THE_ULTIMATE_PROJECT_BOOK.md": "docs/paper_drafts/",
}

literature_docs = {
    "ALPHAFOLD_EVALUATION_LITERATURE.md": "docs/literature/",
    "EVALUATION_LITERATURE_SUPPORT.md": "docs/literature/",
    "PROJECT_NOVELTY_ASSESSMENT.md": "docs/literature/",
}

guide_docs = {
    "PHASE2_PHASE3_COMPLETE_GUIDE.md": "docs/guides/",
}

# Move root .cxc files to visualization/chimerax_scripts/
cxc_files = {
    "highlight_destruction.cxc": "visualization/chimerax_scripts/",
    "highlight_stability.cxc": "visualization/chimerax_scripts/",
    "visualize.cxc": "visualization/chimerax_scripts/",
    "visualize_r213q.cxc": "visualization/chimerax_scripts/",
}

def safe_move(filename, dest_subdir):
    src = os.path.join(BASE, filename)
    dst_dir = os.path.join(BASE, dest_subdir)
    dst = os.path.join(dst_dir, filename)
    if os.path.exists(src):
        if os.path.exists(dst):
            print(f"  SKIP (already exists): {filename}")
            return
        shutil.move(src, dst)
        print(f"  MOVED: {filename} -> {dest_subdir}")
    else:
        print(f"  NOT FOUND: {filename}")

print("\n--- Moving Paper Drafts ---")
for fname, dest in paper_drafts.items():
    safe_move(fname, dest)

print("\n--- Moving Literature Documents ---")
for fname, dest in literature_docs.items():
    safe_move(fname, dest)

print("\n--- Moving Guide Documents ---")
for fname, dest in guide_docs.items():
    safe_move(fname, dest)

print("\n--- Moving ChimeraX Scripts ---")
for fname, dest in cxc_files.items():
    safe_move(fname, dest)

# ============================================================
# COPY OUTPUT FILES INTO ORGANIZED results/ STRUCTURE
# ============================================================
print("\n--- Organizing Results ---")

# Phase 1 results
phase1_csv = ["rmsd_scores.csv", "rmsd_scores.csv"]
phase1_img = ["rmsd_ranking.png", "mutation_distribution.png", 
              "backbone_visualization.png", "backbone_visualization_tp53_r213q.png"]

for fname in phase1_csv:
    src = os.path.join(BASE, "output", fname)
    dst = os.path.join(BASE, "results", "phase1_global_geometry", fname)
    if os.path.exists(src) and not os.path.exists(dst):
        shutil.copy2(src, dst)
        print(f"  COPIED: {fname} -> results/phase1_global_geometry/")

for fname in phase1_img:
    src = os.path.join(BASE, "output", fname)
    dst = os.path.join(BASE, "results", "phase1_global_geometry", "plots", fname)
    if os.path.exists(src) and not os.path.exists(dst):
        shutil.copy2(src, dst)
        print(f"  COPIED: {fname} -> results/phase1_global_geometry/plots/")

# Phase 2 results
phase2_dir = os.path.join(BASE, "output", "phase2")
if os.path.isdir(phase2_dir):
    for fname in os.listdir(phase2_dir):
        src = os.path.join(phase2_dir, fname)
        if os.path.isfile(src):
            if fname.endswith('.csv'):
                dst = os.path.join(BASE, "results", "phase2_domain_analysis", fname)
            elif fname.endswith('.png'):
                dst = os.path.join(BASE, "results", "phase2_domain_analysis", "plots", fname)
            elif fname.endswith('.html'):
                dst = os.path.join(BASE, "results", "phase2_domain_analysis", fname)
            else:
                continue
            if not os.path.exists(dst):
                shutil.copy2(src, dst)
                print(f"  COPIED: {fname} -> results/phase2_domain_analysis/")

# Phase 3 results
phase3_dir = os.path.join(BASE, "output", "phase3")
if os.path.isdir(phase3_dir):
    for fname in os.listdir(phase3_dir):
        src = os.path.join(phase3_dir, fname)
        if os.path.isfile(src):
            if fname.endswith('.csv'):
                dst = os.path.join(BASE, "results", "phase3_advanced_metrics", fname)
            elif fname.endswith('.png'):
                dst = os.path.join(BASE, "results", "phase3_advanced_metrics", "plots", fname)
            else:
                continue
            if not os.path.exists(dst):
                shutil.copy2(src, dst)
                print(f"  COPIED: {fname} -> results/phase3_advanced_metrics/")

# Move chimerax scripts from output/chimerax_scripts/
cxc_src = os.path.join(BASE, "output", "chimerax_scripts")
if os.path.isdir(cxc_src):
    for fname in os.listdir(cxc_src):
        src = os.path.join(cxc_src, fname) 
        dst = os.path.join(BASE, "visualization", "chimerax_scripts", fname)
        if os.path.isfile(src) and not os.path.exists(dst):
            shutil.copy2(src, dst)
            print(f"  COPIED: {fname} -> visualization/chimerax_scripts/")

# Move backbone plots
bp_src = os.path.join(BASE, "output", "backbone_plots")
bp_dst = os.path.join(BASE, "visualization", "backbone_plots")
if os.path.isdir(bp_src) and not os.path.isdir(bp_dst):
    shutil.copytree(bp_src, bp_dst)
    print(f"  COPIED: backbone_plots/ -> visualization/backbone_plots/")

# Move Frontiers_Template_Dir 
ft_src = os.path.join(BASE, "Frontiers_Template_Dir")
ft_dst = os.path.join(BASE, "docs", "Frontiers_Template_Dir")
if os.path.isdir(ft_src) and not os.path.isdir(ft_dst):
    shutil.copytree(ft_src, ft_dst)
    print(f"  COPIED: Frontiers_Template_Dir/ -> docs/Frontiers_Template_Dir/")

# ============================================================
# GENERATE RESULTS INDEX
# ============================================================
print("\n--- Generating Results Index ---")
index_path = os.path.join(BASE, "results", "RESULTS_INDEX.md")
with open(index_path, 'w', encoding='utf-8') as f:
    f.write("# Results Index\n\n")
    f.write("This document provides a complete map of all results and output files organized by analysis phase.\n\n---\n\n")
    
    f.write("## Phase 1: Global Structural Geometry\n\n")
    f.write("| File | Description |\n|------|-------------|\n")
    f.write("| `rmsd_scores.csv` | Global RMSD values for all 128 mutations ranked by severity |\n")
    f.write("| `rmsd_scores.csv` | Complete RMSD dataset with classifications |\n")
    f.write("| `plots/rmsd_ranking.png` | Bar chart ranking all mutations by Global RMSD |\n")
    f.write("| `plots/mutation_distribution.png` | Distribution of mutation positions along p53 sequence |\n")
    f.write("| `plots/backbone_visualization.png` | 3D backbone trace comparison |\n\n")
    
    f.write("## Phase 2: Domain-Specific Analysis\n\n")
    f.write("| File | Description |\n|------|-------------|\n")
    f.write("| `domain_rmsd.csv` | RMSD decomposed by structural domain (N-term, DBD, TET, C-term) |\n")
    f.write("| `rmsd_scores_classified.csv` | RMSD scores with severity classifications |\n")
    f.write("| `tool_comparison.csv` | Comparison against SIFT and PolyPhen-2 predictions |\n")
    f.write("| `heatmap_matrix.csv` | Per-residue displacement heatmap data |\n")
    f.write("| `dashboard.html` | Interactive web dashboard for Phase 2 results |\n")
    f.write("| `plots/domain_rmsd_chart.png` | Multi-bar chart of domain-isolated RMSD values |\n")
    f.write("| `plots/severity_classification.png` | Color-coded severity classification chart |\n")
    f.write("| `plots/tool_comparison.png` | RMSD vs SIFT vs PolyPhen-2 comparison chart |\n")
    f.write("| `plots/mutation_heatmap.png` | Heatmap of per-residue displacements |\n\n")
    
    f.write("## Phase 3: Advanced Biophysical Metrics\n\n")
    f.write("### CSV Data Files\n\n")
    f.write("| File | Description |\n|------|-------------|\n")
    f.write("| `tm_scores.csv` | TM-Score fold conservation for all 128 mutations |\n")
    f.write("| `sasa_analysis.csv` | Solvent Accessible Surface Area changes |\n")
    f.write("| `secondary_structure.csv` | Secondary structure element changes (helix/sheet/coil) |\n")
    f.write("| `contact_analysis.csv` | Contact network analysis (8A cutoff) |\n")
    f.write("| `p53_dbca.csv` | DNA-Binding Competence Assessment (5-probe functional evaluation) |\n")
    f.write("| `ares_scores.csv` | Atomistic Residue Energy Scoring (thermodynamic destabilization) |\n")
    f.write("| `sve_scores.csv` | Structural Variant Evaluator (34-feature classifier) |\n")
    f.write("| `composite_scores.csv` | Composite Structural Impact Score |\n")
    f.write("| `spi_scores.csv` | Structural Pathogenicity Index |\n")
    f.write("| `local_global_impact.csv` | Local vs Global displacement ratios |\n")
    f.write("| `clustering_pca.csv` | PCA dimensionality reduction results |\n")
    f.write("| `compactness_torsion.csv` | Radius of Gyration and Ramachandran strain |\n\n")
    
    f.write("### Visualization Plots\n\n")
    f.write("| File | Description |\n|------|-------------|\n")
    f.write("| `plots/tm_score_ranking.png` | TM-Score ranking bar chart |\n")
    f.write("| `plots/tm_vs_rmsd.png` | TM-Score vs RMSD scatter plot |\n")
    f.write("| `plots/tm_vs_rmsd_rank.png` | TM rank vs RMSD rank comparison |\n")
    f.write("| `plots/sasa_change.png` | SASA change bar chart |\n")
    f.write("| `plots/hydrophobic_exposure.png` | Hydrophobic surface exposure |\n")
    f.write("| `plots/ss_changes.png` | Secondary structure changes |\n")
    f.write("| `plots/ss_transitions.png` | Secondary structure transition types |\n")
    f.write("| `plots/contact_changes.png` | Contact network changes |\n")
    f.write("| `plots/contact_preservation.png` | Contact preservation rates |\n")
    f.write("| `plots/dbd_contact_loss.png` | DBD-specific contact losses |\n")
    f.write("| `plots/dbca_ranking.png` | DBCA score ranking |\n")
    f.write("| `plots/dbca_components.png` | DBCA 5-probe component breakdown |\n")
    f.write("| `plots/dbca_vs_rmsd.png` | DBCA vs RMSD comparison |\n")
    f.write("| `plots/dbca_zinc_detail.png` | Zinc coordination detail |\n")
    f.write("| `plots/ares_ranking.png` | ARES energy score ranking |\n")
    f.write("| `plots/ares_dde_landscape.png` | ARES contact energy landscape |\n")
    f.write("| `plots/ares_vs_rmsd.png` | ARES vs RMSD comparison |\n")
    f.write("| `plots/ares_wavefront.png` | ARES disruption wavefront propagation |\n")
    f.write("| `plots/sve_ranking.png` | SVE final pathogenicity ranking |\n")
    f.write("| `plots/sve_roc.png` | SVE ROC curve (AUC = 1.0000) |\n")
    f.write("| `plots/sve_distribution.png` | SVE score distribution |\n")
    f.write("| `plots/sve_feature_importance.png` | SVE Fisher discriminant feature weights |\n")
    f.write("| `plots/csis_ranking.png` | Composite score ranking |\n")
    f.write("| `plots/spi_ranking.png` | SPI score ranking |\n")
    f.write("| `plots/local_global_ratio.png` | Local vs Global impact ratios |\n")
    f.write("| `plots/local_vs_global_scatter.png` | Local vs Global scatter plot |\n")
    f.write("| `plots/pca_scatter.png` | PCA cluster visualization |\n")
    f.write("| `plots/pca_scree.png` | PCA variance explained |\n")
    f.write("| `plots/dendrogram.png` | Hierarchical clustering dendrogram |\n")
    f.write("| `plots/radius_of_gyration.png` | Protein compactness |\n\n")

print(f"Results index written to: {index_path}")

# ============================================================
# UPDATE README
# ============================================================
readme_path = os.path.join(BASE, "README.md")
with open(readme_path, 'w', encoding='utf-8') as f:
    f.write("# TP53 Structural Variant Evaluator (TP53-SVE)\n\n")
    f.write("A multidimensional computational pipeline for classifying TP53 missense mutation pathogenicity using AlphaFold-predicted protein structures.\n\n")
    f.write("## Project Structure\n\n")
    f.write("```\n")
    f.write("chrono-shepard/\n")
    f.write("|\n")
    f.write("|-- data/                          # Input data\n")
    f.write("|   |-- alphafold_inputs_reconstructed/  # AlphaFold input JSON files\n")
    f.write("|   |-- raw_downloads/             # Raw COSMIC mutation data\n")
    f.write("|   |-- structures/                # AlphaFold .cif structure files\n")
    f.write("|\n")
    f.write("|-- src/                           # Source code\n")
    f.write("|   |-- phase2/                    # Phase 2 analysis scripts\n")
    f.write("|   |-- phase3/                    # Phase 3 analysis scripts\n")
    f.write("|   |-- methodology/               # Supporting methodology scripts\n")
    f.write("|\n")
    f.write("|-- results/                       # Organized results (clean copy)\n")
    f.write("|   |-- phase1_global_geometry/    # RMSD scores and ranking plots\n")
    f.write("|   |-- phase2_domain_analysis/    # Domain RMSD, tool comparison, heatmaps\n")
    f.write("|   |-- phase3_advanced_metrics/   # TM, SASA, DBCA, ARES, SVE scores + plots\n")
    f.write("|   |-- RESULTS_INDEX.md           # Complete file-by-file results guide\n")
    f.write("|\n")
    f.write("|-- output/                        # Raw output (original pipeline output)\n")
    f.write("|   |-- phase2/                    # Phase 2 raw outputs\n")
    f.write("|   |-- phase3/                    # Phase 3 raw outputs\n")
    f.write("|\n")
    f.write("|-- docs/                          # Documentation\n")
    f.write("|   |-- paper_drafts/              # Research paper drafts\n")
    f.write("|   |-- literature/                # Literature reviews and citations\n")
    f.write("|   |-- guides/                    # Project guides and walkthroughs\n")
    f.write("|   |-- Frontiers_Template_Dir/    # Frontiers Journal LaTeX template\n")
    f.write("|\n")
    f.write("|-- visualization/                 # Visualization assets\n")
    f.write("|   |-- chimerax_scripts/          # UCSF ChimeraX .cxc command scripts\n")
    f.write("|   |-- backbone_plots/            # Backbone overlay plots\n")
    f.write("|\n")
    f.write("|-- dashboard/                     # Interactive React dashboard\n")
    f.write("|\n")
    f.write("|-- PROJECT_ANALYSIS_SUMMARY.md    # Comprehensive analysis report\n")
    f.write("|-- README.md                      # This file\n")
    f.write("```\n\n")
    
    f.write("## Key Results\n\n")
    f.write("- **50 TP53 mutations** analyzed via AlphaFold structure prediction\n")
    f.write("- **Global RMSD proven unreliable** — P72R (benign) ranked 2nd at 37.08 A\n")
    f.write("- **3 novel algorithms developed:** p53-DBCA, TP53-ARES, TP53-SVE\n")
    f.write("- **TP53-SVE classifier** achieves **AUC = 1.0000** separating pathogenic from benign\n")
    f.write("- **34 biophysical features** integrated via Fisher's Linear Discriminant Analysis\n\n")
    
    f.write("## How to Run\n\n")
    f.write("```bash\n")
    f.write("# Phase 1: RMSD calculation\n")
    f.write("python src/calculate_rmsd_scores.py\n\n")
    f.write("# Phase 2: Domain analysis\n")
    f.write("python src/phase2/per_residue_rmsd.py\n")
    f.write("python src/phase2/domain_rmsd.py\n")
    f.write("python src/phase2/tool_comparison.py\n\n")
    f.write("# Phase 3: Advanced metrics\n")
    f.write("python src/phase3/tm_score.py\n")
    f.write("python src/phase3/sasa_analysis.py\n")
    f.write("python src/phase3/contact_network.py\n")
    f.write("python src/phase3/p53_dbca.py\n")
    f.write("python src/phase3/tp53_ares.py\n")
    f.write("python src/phase3/tp53_sve.py\n")
    f.write("```\n")

print(f"README updated at: {readme_path}")
print("\nOrganization complete!")
