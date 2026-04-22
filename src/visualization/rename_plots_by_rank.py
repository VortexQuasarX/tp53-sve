import pandas as pd
import os
import shutil

def main():
    csv_file = "output/rmsd_scores.csv"
    plots_dir = "output/backbone_plots"
    
    if not os.path.exists(csv_file):
        print("RMSD Scores file not found.")
        return

    # 1. Load and Sort Data
    df = pd.read_csv(csv_file)
    # Sort DESCENDING (Highest RMSD is #1)
    df = df.sort_values(by="RMSD (Angstroms)", ascending=False)
    
    # 2. Rename Files
    print(f"Renaming files in {plots_dir} based on damage ranking...")
    
    for rank, (index, row) in enumerate(df.iterrows(), start=1):
        mutation = row['Mutation'] # e.g. R175H
        rmsd = row['RMSD (Angstroms)']
        
        # Construct original filename pattern
        # Our plots are named: plot_tp53_r175h.png
        gene_lower = "tp53"
        mut_lower = mutation.lower()
        original_filename = f"plot_{gene_lower}_{mut_lower}.png"
        original_path = os.path.join(plots_dir, original_filename)
        
        # Construct new filename
        # Format: #01_plot_tp53_r175h.png (Zero padded for sorting)
        new_filename = f"#{rank:02d}_{original_filename}"
        new_path = os.path.join(plots_dir, new_filename)
        
        if os.path.exists(original_path):
            # Check if already renamed to avoid double prefixing if run twice
            if not os.path.exists(new_path):
                os.rename(original_path, new_path)
                print(f"Rank {rank}: {original_filename} -> {new_filename} (RMSD: {rmsd})")
            else:
                 print(f"Skipping {original_filename}, target {new_filename} already exists.")
        else:
            # Maybe it was already renamed? Try to find existing numbered file
            pass

    print("\n[SUCCESS] Files renamed. Sorted by destruction level.")

if __name__ == "__main__":
    main()
