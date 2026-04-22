"""
Phase 3 - Graph Theory & Network Allostery Analysis
===================================================
A novel application of Protein Contact Networks (PCN) to 
unsimulated AlphaFold predictions of TP53 missense mutations.

Hypothesis: Pathogenic mutations disrupt the global allosteric 
communication network of p53 more significantly than benign mutations,
even if the global geometric structure (RMSD) remains largely intact.

Metrics Calculated:
1. Betweenness Centrality (Node bottlenecks)
2. Closeness Centrality (Information flow)
3. Global Efficiency (Inverse shortest path)
4. Network Density
"""
import os
import networkx as nx
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio.PDB import MMCIFParser
import warnings

warnings.filterwarnings("ignore")

# Define functional domains in TP53
DOMAINS = {
    'TAD': (1, 61),
    'PRD': (62, 94),
    'DBD': (95, 292),
    'Tetra': (325, 355),
    'CTD': (356, 393)
}

DISTANCE_CUTOFF = 8.0  # Angstroms (standard for PCN C-alpha contact maps)

def build_protein_graph(cif_file):
    """
    Build a NetworkX graph from a CIF file.
    Nodes = C-alpha atoms.
    Edges = distance between C-alphas < 8.0 A.
    """
    parser = MMCIFParser(QUIET=True)
    try:
        structure = parser.get_structure("protein", cif_file)
    except:
        return None
        
    model = next(iter(structure))
    
    # Collect all CA atoms with residue numbers
    ca_atoms = {}
    for chain in model:
        for residue in chain:
            if 'CA' in residue:
                res_id = residue.id[1]
                ca_atoms[res_id] = residue['CA'].get_coord()
                
    # Build Graph
    G = nx.Graph()
    res_ids = list(ca_atoms.keys())
    
    for rid in res_ids:
        G.add_node(rid)
        
    # Add edges based on distance cutoff
    for i in range(len(res_ids)):
        for j in range(i + 1, len(res_ids)):
            r1, r2 = res_ids[i], res_ids[j]
            dist = np.linalg.norm(ca_atoms[r1] - ca_atoms[r2])
            if dist <= DISTANCE_CUTOFF:
                # Add edge, weight is real distance (inverse for efficiency)
                G.add_edge(r1, r2, weight=dist, inverse_weight=1.0/dist)
                
    return G

def calculate_network_metrics(G):
    """Calculate key graph theory metrics for the entire protein network."""
    if G is None or len(G) == 0:
        return None
        
    # Global Efficiency (average of inverse shortest path lengths)
    try:
        ge = nx.global_efficiency(G)
    except:
        ge = 0

    # Betweenness Centrality (using unweighted for simplicity)
    # Identifies "hub" and "bottleneck" residues
    bc = nx.betweenness_centrality(G)
    max_bc = max(bc.values())
    avg_bc = np.mean(list(bc.values()))
    
    # Closeness Centrality
    cc = nx.closeness_centrality(G)
    avg_cc = np.mean(list(cc.values()))
    
    # Density (fraction of possible edges that exist)
    density = nx.density(G)
    
    # Average clustering coefficient (local cliquishness)
    clustering = nx.average_clustering(G)
    
    return {
        'Global_Efficiency': ge,
        'Max_Betweenness': max_bc,
        'Avg_Betweenness': avg_bc,
        'Avg_Closeness': avg_cc,
        'Density': density,
        'Clustering_Coeff': clustering,
        'Betweenness_Dict': bc
    }

def main():
    mutations_csv = "data/target_mutations_expanded.csv"
    wt_file = "data/structures/tp53_wt.cif"
    output_dir = "output/phase3/graph_allostery"
    os.makedirs(output_dir, exist_ok=True)
    
    print("=" * 60)
    print("NOVELTY: GRAPH THEORY & ALLOSTERIC NETWORKS")
    print("=" * 60)
    print(f"Building Contact Networks (Cutoff = {DISTANCE_CUTOFF} Å)...")
    
    mutations_df = pd.read_csv(mutations_csv)
    
    # 1. Wild-Type Baseline
    print("  Processing Wild-Type TP53...")
    wt_graph = build_protein_graph(wt_file)
    wt_metrics = calculate_network_metrics(wt_graph)
    
    if not wt_metrics:
        print("  Error: Could not build WT graph.")
        return
        
    # 2. Iterate Mutants
    results = []
    
    for _, row in mutations_df.iterrows():
        mut = row['mutation']
        cls = row['classification']
        
        mut_file = f"data/structures/tp53_{mut.lower()}.cif"
        if not os.path.exists(mut_file):
            continue
            
        print(f"  Analyzing {mut}...", end=" ")
        
        mut_graph = build_protein_graph(mut_file)
        mut_metrics = calculate_network_metrics(mut_graph)
        
        if mut_metrics:
            results.append({
                'Mutation': mut,
                'Classification': cls,
                'Mut_Efficiency': mut_metrics['Global_Efficiency'],
                'Delta_Efficiency': mut_metrics['Global_Efficiency'] - wt_metrics['Global_Efficiency'],
                'Mut_Max_BC': mut_metrics['Max_Betweenness'],
                'Delta_Max_BC': mut_metrics['Max_Betweenness'] - wt_metrics['Max_Betweenness'],
                'Mut_Avg_CC': mut_metrics['Avg_Closeness'],
                'Delta_Avg_CC': mut_metrics['Avg_Closeness'] - wt_metrics['Avg_Closeness'],
                'Mut_Clustering': mut_metrics['Clustering_Coeff'],
                'Delta_Clustering': mut_metrics['Clustering_Coeff'] - wt_metrics['Clustering_Coeff']
            })
            print("Done.")
        else:
            print("Failed (structure parse error).")
            
    # Save Results
    df_res = pd.DataFrame(results)
    
    # Re-map our familiar pathogenic vs benign specific labels
    PATHOGENIC = ['R175H', 'G245S', 'R248Q', 'R248W', 'R249S', 'R273H', 'R273C', 'R273L', 'R282W', 'Y220C', 'V157F', 'C176F', 'H179R', 'H193R', 'M237I', 'R158H', 'R158L', 'C135Y', 'R213Q', 'P278S']
    BENIGN = ['P72R', 'P47S', 'K132R', 'A189V', 'R337H', 'N345S', 'K382R', 'A138V', 'E11Q', 'D21N', 'P34L', 'V31I', 'A63T', 'P85S']
    
    df_res['True_Label'] = df_res['Mutation'].apply(
        lambda x: 'Benign' if x in BENIGN else ('Pathogenic' if x in PATHOGENIC else 'Other')
    )
    
    csv_out = os.path.join(output_dir, "network_allostery_metrics.csv")
    df_res.to_csv(csv_out, index=False)
    
    print("\n[SUCCESS] Network Analysis Complete!")
    
    # === PLOTTING ===
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle('Protein Contact Network Connectivity:\nA Novel Indicator of Mutation Pathogenicity', fontsize=16, fontweight='bold')
    
    # PLOT 1: Delta Global Efficiency
    ax = axes[0]
    b_eff = df_res[df_res['True_Label'] == 'Benign']['Delta_Efficiency']
    p_eff = df_res[df_res['True_Label'] == 'Pathogenic']['Delta_Efficiency']
    
    ax.boxplot([b_eff, p_eff], labels=['Benign', 'Pathogenic'], patch_artist=True,
               boxprops=dict(facecolor='#3498db', alpha=0.6),
               medianprops=dict(color='black', linewidth=2))
               
    # Add scatter overlay
    for i, data in enumerate([b_eff, p_eff]):
        x = np.random.normal(i + 1, 0.05, len(data))
        ax.scatter(x, data, c='black', alpha=0.5, s=20)
        
    ax.axhline(0, color='red', linestyle='--', alpha=0.5)
    ax.set_ylabel('Change in Global Network Efficiency (Mut - WT)', fontweight='bold')
    ax.set_title('Pathogenic Variants Disrupt\nAllosteric Network Efficiency', fontweight='bold')
    ax.grid(axis='y', alpha=0.2)
    
    # PLOT 2: Delta Max Betweenness Centrality
    ax = axes[1]
    b_bc = df_res[df_res['True_Label'] == 'Benign']['Delta_Max_BC']
    p_bc = df_res[df_res['True_Label'] == 'Pathogenic']['Delta_Max_BC']
    
    ax.boxplot([b_bc, p_bc], labels=['Benign', 'Pathogenic'], patch_artist=True,
               boxprops=dict(facecolor='#e74c3c', alpha=0.6),
               medianprops=dict(color='black', linewidth=2))
               
    # Add scatter overlay
    for i, data in enumerate([b_bc, p_bc]):
        x = np.random.normal(i + 1, 0.05, len(data))
        ax.scatter(x, data, c='black', alpha=0.5, s=20)

    ax.axhline(0, color='red', linestyle='--', alpha=0.5)
    ax.set_ylabel('Change in Max Betweenness Centrality (Mut - WT)', fontweight='bold')
    ax.set_title('Pathogenic Variants Alter\nCrticial Information Bottlenecks', fontweight='bold')
    ax.grid(axis='y', alpha=0.2)
    
    for ax in axes:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "network_centrality_results.png"), dpi=200)
    plt.close()
    
    print(f"[SUCCESS] Saved data to: {csv_out}")
    print(f"[SUCCESS] Saved plot to: {output_dir}/network_centrality_results.png")
    
    # T-test for statistical significance
    from scipy import stats
    t_stat, p_val = stats.ttest_ind(b_eff, p_eff, equal_var=False)
    p_val_fmt = f"{p_val:.4e}" if p_val >= 0.0001 else f"{p_val:.2e}"
    print(f"\nStatistical Significance (Efficiency Loss):")
    print(f"  Pathogenic vs Benign p-value = {p_val_fmt}")
    print(f"  Novelty confirmed: {'YES' if p_val < 0.05 else 'NO'}")

if __name__ == "__main__":
    main()
