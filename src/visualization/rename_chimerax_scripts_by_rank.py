import pandas as pd
import os
import shutil

def main():
    csv_file = "output/rmsd_scores.csv"
    scripts_dir = "output/chimerax_scripts"
    
    if not os.path.exists(csv_file):
        print("RMSD Scores file not found.")
        return

    # 1. Load and Sort Data
    df = pd.read_csv(csv_file)
    # Sort DESCENDING (Highest RMSD is #1)
    df = df.sort_values(by="RMSD (Angstroms)", ascending=False)
    
    # 2. Rename Files
    print(f"Renaming files in {scripts_dir} based on damage ranking...")
    
    for rank, (index, row) in enumerate(df.iterrows(), start=1):
        mutation = row['Mutation'] # e.g. R175H
        
        # Construct original filename pattern
        # Our scripts are named: visualize_tp53_r175h.cxc
        gene_lower = "tp53"
        mut_lower = mutation.lower()
        original_filename = f"visualize_{gene_lower}_{mut_lower}.cxc"
        original_path = os.path.join(scripts_dir, original_filename)
        
        # Construct new filename
        # Format: #01_visualize_tp53_r175h.cxc (Zero padded for sorting)
        new_filename = f"#{rank:02d}_{original_filename}"
        new_path = os.path.join(scripts_dir, new_filename)
        
        if os.path.exists(original_path):
            if not os.path.exists(new_path):
                os.rename(original_path, new_path)
                print(f"Rank {rank}: {original_filename} -> {new_filename}")
            else:
                 print(f"Skipping {original_filename}, target {new_filename} already exists.")
        else:
             # Check if it was already renamed (starts with #rank)
             # This simple check avoids errors if run multiple times
             pass

    print("\n[SUCCESS] ChimeraX scripts renamed and sorted.")

if __name__ == "__main__":
    main()
