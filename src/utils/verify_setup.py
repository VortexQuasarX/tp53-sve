import os
import sys

def verify_dataset():
    print("="*60)
    print("TP53-SVE: Dataset Integrity & Setup Verifier")
    print("="*60)

    base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    struct_dir = os.path.join(base_dir, "data", "structures")
    
    # Core baseline samples that should be present
    baseline_samples = ["tp53_wt.cif", "tp53_r175h.cif", "tp53_y220c.cif"]
    missing_baseline = []

    if not os.path.exists(struct_dir):
        print(f"[!] Critical: Data directory not found at {struct_dir}")
        os.makedirs(struct_dir, exist_ok=True)
        print("[+] Created data/structures/ directory.")
    
    for sample in baseline_samples:
        if not os.path.exists(os.path.join(struct_dir, sample)):
            missing_baseline.append(sample)

    if missing_baseline:
        print(f"[!] Missing Baseline Samples: {', '.join(missing_baseline)}")
        print("    Please ensure you have cloned the repository correctly.")
    else:
        print("[✔] Baseline Samples (WT, R175H, Y220C) are PRESENT.")
        print("    You can run the dashboard and basic analysis immediately.")

    # Check for the full 1.3GB database (approx 128 variants)
    all_files = [f for f in os.listdir(struct_dir) if f.endswith('.cif')] if os.path.exists(struct_dir) else []
    
    print("-" * 60)
    print(f"Current Structure Count: {len(all_files)} / 128 variants detected.")
    
    if len(all_files) < 100:
        print("\n[INFO] Full Structural Database (1.3GB) is currently MISSING.")
        print("To run analysis on all 128 clinical variants, follow these steps:")
        print("1. Download the database from: [INSERT_LINK_HERE_IN_README]")
        print("2. Extract all .cif files into: data/structures/")
        print("3. Rerun this script to verify.")
    else:
        print("\n[✔] Full Structural Database is PRESENT.")
        print("    You have full capability to rerun the scientific analysis pipeline.")
    
    print("-" * 60)

if __name__ == "__main__":
    verify_dataset()
