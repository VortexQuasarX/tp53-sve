import json
import os

def show_mutation_data(mutant_name="TP53_P278S", position=278):
    # 1. Load the Mutant File
    filepath = f"data/alphafold_inputs_reconstructed/{mutant_name}.json"
    if not os.path.exists(filepath):
        print(f"File not found: {filepath}")
        return

    with open(filepath, 'r') as f:
        data = json.load(f)

    # 2. Extract the Sequence
    # JSON structure: names -> sequences -> proteinChain -> sequence
    full_seq = data['sequences'][0]['proteinChain']['sequence']
    
    # 3. Define Window to Show (e.g., +/- 10 residues around mutation)
    # Note: position is 1-indexed, python string is 0-indexed
    idx = position - 1
    start = max(0, idx - 10)
    end = min(len(full_seq), idx + 11)
    
    snippet = full_seq[start:end]
    
    # 4. Print the "Visual Data"
    print(f"\n--- DATA INSIDE FILE: {mutant_name}.json ---")
    print(f"Full Sequence Length: {len(full_seq)}")
    print(f"Focusing on Residue {position}...\n")
    
    print(f"Sequence Context (Residues {start+1}-{end}):")
    print(f"  {snippet}")
    
    # Create a pointer to the mutation
    pointer = " " * (idx - start) + "^"
    print(f"  {pointer}")
    print(f"  Mutation Site: {full_seq[idx]} (Serine)")
    print("------------------------------------------------")
    print("This 'S' is the only byte different from the Wild Type.\n")

if __name__ == "__main__":
    # Demo for P278S
    show_mutation_data("TP53_P278S", 278)
    
    # Demo for R273H
    show_mutation_data("TP53_R273H", 273)
