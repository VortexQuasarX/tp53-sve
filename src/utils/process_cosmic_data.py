import pandas as pd
import os

def process_cosmic_data():
    # 1. Define the Raw Source File (Downloaded from COSMIC Website)
    raw_file = "data/raw_downloads/CosmicMutantExport_Sample.csv"
    output_file = "data/target_mutations_from_cosmic.csv"
    
    print(f"Reading raw COSMIC database dump: {raw_file}...")
    
    # 2. Load the Data
    df = pd.read_csv(raw_file)
    print(f"Total rows in raw file: {len(df)}")
    
    # 3. Filter for TP53 Gene ONLY
    print("Filtering for Gene = 'TP53'...")
    tp53_df = df[df['Gene name'] == 'TP53'].copy()
    
    # 4. Filter for our specific target mutations (The 'Frequent' ones)
    # in a real scenario, you might filter by frequency count
    print("Extracting relevant columns (Mutation AA, Mutation CDS, ID)...")
    
    # Clean up the 'Mutation AA' column (remove 'p.')
    tp53_df['clean_mutation'] = tp53_df['Mutation AA'].str.replace('p.', '', regex=False)
    
    cleaned_list = tp53_df[['Gene name', 'clean_mutation', 'Mutation AA', 'Mutation ID']].drop_duplicates()
    
    # 5. Save the Clean "Target List"
    cleaned_list.to_csv(output_file, index=False)
    print(f"Successfully extracted {len(cleaned_list)} mutations.")
    print(f"Saved processed list to: {output_file}")
    
    # Show the result
    print("\n--- Extracted Data Preview ---")
    print(cleaned_list.to_string(index=False))

if __name__ == "__main__":
    process_cosmic_data()
