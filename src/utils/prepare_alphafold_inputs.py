import pandas as pd
import json
import os

def main():
    # 1. Load the Target List (The "Source Data")
    csv_file = "data/target_mutations_expanded.csv"
    if not os.path.exists(csv_file):
        print("Error: Target mutation list not found!")
        return
    
    df = pd.read_csv(csv_file)
    print(f"Loaded {len(df)} mutations from {csv_file}")

    # 2. Define the Wild-Type Sequence (TP53 DNA Binding Domain)
    # (Restored from project history)
    wt_seq = "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELPPGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPGGSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD"
    
    output_dir = "data/alphafold_inputs_reconstructed"
    os.makedirs(output_dir, exist_ok=True)

    # 3. Generate Inputs for each Mutation
    for _, row in df.iterrows():
        mutant_name = f"{row['gene']}_{row['mutation']}"
        pos = row['position'] - 1 # 0-indexed
        mut_res = row['mut_residue']
        
        # Apply mutation
        mut_seq = list(wt_seq)
        if pos < len(mut_seq):
            mut_seq[pos] = mut_res
            mut_seq_str = "".join(mut_seq)
            
            # Create the JSON structure AlphaFold expects
            job_data = {
                "name": mutant_name,
                "modelSeeds": [],
                "sequences": [
                    {
                        "proteinChain": {
                            "sequence": mut_seq_str,
                            "count": 1
                        }
                    }
                ]
            }
            
            # Save format
            filename = os.path.join(output_dir, f"{mutant_name}.json")
            with open(filename, 'w') as f:
                json.dump(job_data, f, indent=2)
            
            print(f"Generated input file: {filename}")
        else:
            print(f"Error: Position {pos} out of range for {mutant_name}")

    print("\n[SUCCESS] All AlphaFold input files have been reconstructed.")

if __name__ == "__main__":
    main()
