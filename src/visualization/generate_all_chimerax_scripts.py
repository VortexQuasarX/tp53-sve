import pandas as pd
import os

def generate_scripts():
    # 1. Load the Data
    csv_file = "data/target_mutations_expanded.csv"
    output_dir = "output/chimerax_scripts"
    os.makedirs(output_dir, exist_ok=True)
    
    if not os.path.exists(csv_file):
        print("Error: data/target_mutations_expanded.csv not found")
        return

    df = pd.read_csv(csv_file)
    
    # 2. Template for the ChimeraX Script
    # Upgraded to include "Sphere Highlighting" (Green vs Yellow)
    template = """# ChimeraX Script for {mutation_name}
# 1. Open the Structures
open ../../data/structures/tp53_wt.cif
open ../../data/structures/{mutant_filename}

# 2. Align them (Match Maker)
matchmaker #2 to #1

# 3. Styling
# We use default Ribbon style for the backbone (Cleanest view)
# Blue = Wild Type, Red = Mutant
color #1 blue
color #2 red

# 4. Focus on the Mutation Site
# Residue {residue_number} on Mutant Chain (Model #2)
show #2:{residue_number} atoms
view #2:{residue_number}

# 5. Highlight the Change (Sphere Style)
# This matches the "Destruction vs Stability" demo style
style #2:{residue_number} sphere
color #2:{residue_number} {sphere_color}

# 6. Coloring the Backbones (Already done in Step 3)
# Blue = Wild Type, Red = Mutant

# Label the site
2dlabels create title text "{mutation_code}: {classification}" xpos 0.5 ypos 0.9 color white size 30

# Lighting
lighting soft
view orient

# 7. Save Image (Optional)
# save ../visualization_{mutation_name}.png width 1200 height 1200
"""

    # Load RMSD scores for color grading
    rmsd_df = pd.read_csv("output/rmsd_scores.csv")
    rmsd_map = dict(zip(rmsd_df['Mutation'], rmsd_df['RMSD (Angstroms)']))

    print(f"Generating scripts in {output_dir}...")

    # 3. Loop through all 20 mutations
    for _, row in df.iterrows():
        # Parse data
        gene = row['gene']
        mutation = row['mutation'] # e.g. R175H
        residue_num = row['position'] # e.g. 175
        classification = row['classification'] # e.g. Hotspot
        
        # Get RMSD
        rmsd = rmsd_map.get(mutation, 0.0)
        
        # Determine Color Scheme based on Damage
        # < 15 A = Green (Stable)
        # 15 - 30 A = Yellow (Moderate)
        # > 30 A = Red (Destructive)
        
        sphere_color = "yellow"
        if rmsd < 15:
            sphere_color = "green" 
        elif rmsd > 30:
            sphere_color = "red"
        
        # Construct filenames
        mutant_name = f"{gene.lower()}_{mutation.lower()}" # tp53_r175h
        mutant_filename = f"{mutant_name}.cif" # tp53_r175h.cif
        script_filename = f"visualize_{mutant_name}.cxc"
        
        # Fill template
        script_content = template.format(
            mutation_name=mutation,
            mutant_filename=mutant_filename,
            residue_number=residue_num,
            mutation_code=mutation,
            classification=classification,
            sphere_color=sphere_color
        )
        
        # Write file
        full_path = os.path.join(output_dir, script_filename)
        with open(full_path, 'w') as f:
            f.write(script_content)
        
        print(f"Created: {script_filename} (RMSD: {rmsd:.1f} -> {sphere_color})")

    print(f"\n[SUCCESS] Generated {len(df)} scripts with dynamic coloring.")

if __name__ == "__main__":
    generate_scripts()
