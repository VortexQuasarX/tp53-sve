"""
Fix AlphaFold JSON format - OFFICIAL AlphaFold Server format.
Based on: https://github.com/google-deepmind/alphafold/tree/main/server

Required format:
  - Top level: array of job objects [...]
  - Each job needs: name, modelSeeds, sequences, dialect, version
  - dialect: "alphafoldserver"
  - version: 1
  - Protein uses "proteinChain" with "sequence" and "count"
"""
import json
import os

# Wild-type p53 sequence (393 aa)
WT_SEQUENCE = (
    "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP"
    "DEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAK"
    "SVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHE"
    "RCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNS"
    "SCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELP"
    "PGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPG"
    "GSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD"
)

MUTATIONS = [
    # A: Same-Position Variants (7)
    ("R175G", 175, "G"), ("R175C", 175, "C"), ("G245D", 245, "D"),
    ("R282Q", 282, "Q"), ("R248L", 248, "L"), ("Y220S", 220, "S"),
    ("R249G", 249, "G"),
    # B: Benign Controls (5)
    ("P72R", 72, "R"), ("P47S", 47, "S"), ("A189V", 189, "V"),
    ("R337H", 337, "H"), ("K132R", 132, "R"),
    # C: Non-DBD Domain (6)
    ("L344R", 344, "R"), ("R342P", 342, "P"), ("L22F", 22, "F"),
    ("W23R", 23, "R"), ("N345S", 345, "S"), ("K382R", 382, "R"),
    # D: Gain-of-Function (4)
    ("R280K", 280, "K"), ("V272M", 272, "M"), ("D281G", 281, "G"),
    ("S241F", 241, "F"),
    # E: Temperature-Sensitive (3)
    ("V143A", 143, "A"), ("A138V", 138, "V"), ("I195T", 195, "T"),
    # F: Rare Pathogenic (5)
    ("E285K", 285, "K"), ("N239D", 239, "D"), ("T125M", 125, "M"),
    ("L194R", 194, "R"), ("N247D", 247, "D"),
]

def mutate(seq, pos, new_aa):
    s = list(seq)
    s[pos - 1] = new_aa
    return "".join(s)

def make_server_job(name, sequence):
    """Create a job in the OFFICIAL AlphaFold Server format."""
    return {
        "name": name,
        "modelSeeds": [],
        "sequences": [
            {
                "proteinChain": {
                    "sequence": sequence,
                    "count": 1
                }
            }
        ],
        "dialect": "alphafoldserver",
        "version": 1
    }

def main():
    output_dir = "data/alphafold_inputs_phase2"
    os.makedirs(output_dir, exist_ok=True)
    
    # Clean old JSON files
    for f in os.listdir(output_dir):
        if f.endswith(".json"):
            os.remove(os.path.join(output_dir, f))
    
    all_jobs = []
    
    print("=" * 60)
    print("ALPHAFOLD SERVER FORMAT (OFFICIAL)")
    print("  dialect: alphafoldserver")
    print("  version: 1")
    print("=" * 60)
    
    for i, (mutation, pos, mut_aa) in enumerate(MUTATIONS, 1):
        seq = mutate(WT_SEQUENCE, pos, mut_aa)
        name = f"TP53_{mutation}"
        job = make_server_job(name, seq)
        
        # Save individual file (wrapped in array as required)
        path = os.path.join(output_dir, f"{name}.json")
        with open(path, "w") as f:
            json.dump([job], f, indent=2)
        
        all_jobs.append(job)
        print(f"  [{i:2d}/30] {name}")
    
    # Save ALL 30 as ONE batch file (array of 30 jobs)
    batch_path = os.path.join(output_dir, "_BATCH_ALL_30.json")
    with open(batch_path, "w") as f:
        json.dump(all_jobs, f, indent=2)
    
    print(f"\n{'='*60}")
    print(f"DONE! All 30 files in OFFICIAL format.")
    print(f"{'='*60}")
    print(f"  Individual: {output_dir}/TP53_*.json")
    print(f"  BATCH:      {output_dir}/_BATCH_ALL_30.json")
    print(f"\n  Upload _BATCH_ALL_30.json → Submit 30 jobs as drafts")

if __name__ == "__main__":
    main()
