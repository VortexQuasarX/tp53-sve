import pandas as pd
import matplotlib.pyplot as plt
import os

def plot_distribution():
    # 1. Load the Data
    csv_file = "data/target_mutations_expanded.csv"
    if not os.path.exists(csv_file):
        print("Error: File not found.")
        return
    
    df = pd.read_csv(csv_file)
    
    # 2. Setup the Plot
    plt.figure(figsize=(12, 6))
    
    # Draw the "Gene" line (TP53 is ~393 residues long)
    plt.hlines(y=0, xmin=0, xmax=393, color='gray', linewidth=4, zorder=1)
    
    # 3. Plot Each Mutation as a "Lollipop"
    # We will use the 'position' column
    # If multiple mutations are at the same position, we stack them slightly or just make the marker bigger
    
    counts = df['position'].value_counts()
    
    for pos, count in counts.items():
        # Get the label (e.g., "R175H")
        muts_at_pos = df[df['position'] == pos]['mutation'].tolist()
        label = "\n".join(muts_at_pos)
        
        # Plot the stem
        plt.vlines(x=pos, ymin=0, ymax=count, color='#d62728', linewidth=2, zorder=2)
        
        # Plot the head (marker)
        plt.plot(pos, count, 'o', color='#d62728', markersize=10, zorder=3)
        
        # Add label txt
        plt.text(pos, count + 0.1, label, ha='center', va='bottom', fontsize=9, rotation=0)

    # 4. Formatting
    plt.title("TP53 Mutation Distribution (The 'Hotspots')", fontsize=14, fontweight='bold')
    plt.xlabel("Residue Position (Amino Acid Index)", fontsize=12)
    plt.ylabel("Frequency in Dataset", fontsize=12)
    plt.xlim(0, 400)
    plt.ylim(0, counts.max() + 1)
    plt.grid(axis='y', linestyle='--', alpha=0.5)
    
    # Add domains (approximate) for context
    # DNA Binding Domain is roughly 100-300
    plt.axvspan(102, 292, color='blue', alpha=0.1, label='DNA Binding Domain')
    plt.legend(loc='upper right')

    plt.tight_layout()
    
    # 5. Save
    output_path = "output/mutation_distribution.png"
    plt.savefig(output_path, dpi=300)
    print(f"Plot saved to: {output_path}")

if __name__ == "__main__":
    plot_distribution()
