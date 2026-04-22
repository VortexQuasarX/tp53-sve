"""
AlphaFold Server Batch Helper
Opens a browser tab for each mutation with the sequence ready to paste.
Run this, then for each tab: paste the sequence, name the job, click Submit.
"""
import json
import os
import webbrowser
import time

def main():
    input_dir = "data/alphafold_inputs_phase2"
    
    json_files = sorted([f for f in os.listdir(input_dir) if f.endswith(".json")])
    
    print("=" * 70)
    print("ALPHAFOLD SERVER BATCH SUBMISSION HELPER")
    print("=" * 70)
    print(f"\nTotal jobs to submit: {len(json_files)}")
    print(f"\nThis will open the AlphaFold Server page.")
    print(f"For each mutation, you'll see the name and sequence below.")
    print(f"Steps for each:")
    print(f"  1. Click 'Add protein' on the AlphaFold Server page")
    print(f"  2. Paste the sequence (copied to clipboard)")
    print(f"  3. Name the job")
    print(f"  4. Click 'Submit' / 'Run'")
    print(f"  5. Press Enter here to get the next one")
    
    # Open AlphaFold Server
    print(f"\n--- Opening AlphaFold Server ---")
    webbrowser.open("https://alphafoldserver.com/jobs")
    time.sleep(2)
    
    for i, json_file in enumerate(json_files, 1):
        filepath = os.path.join(input_dir, json_file)
        with open(filepath) as f:
            data = json.load(f)
        
        name = data["name"]
        sequence = data["sequences"][0]["proteinChain"]["sequence"]
        
        print(f"\n{'='*70}")
        print(f"  JOB {i}/30: {name}")
        print(f"{'='*70}")
        print(f"  Sequence length: {len(sequence)} aa")
        print(f"  Sequence (copy this):")
        print(f"  {sequence}")
        
        # Try to copy to clipboard
        try:
            import subprocess
            process = subprocess.Popen(['clip'], stdin=subprocess.PIPE)
            process.communicate(sequence.encode('utf-8'))
            print(f"\n  ✅ Sequence copied to clipboard! Just paste it.")
        except:
            print(f"\n  ⚠️  Couldn't auto-copy. Please copy the sequence above manually.")
        
        if i < len(json_files):
            input(f"\n  Press ENTER when you've submitted {name} (→ next: {json_files[i].replace('.json','')})...")
        else:
            input(f"\n  Press ENTER when you've submitted the LAST job ({name})...")
    
    print(f"\n{'='*70}")
    print(f"ALL 30 JOBS SUBMITTED! 🎉")
    print(f"{'='*70}")
    print(f"Wait for AlphaFold to complete all jobs (usually 1-5 minutes each).")
    print(f"Then download each .cif file and place in data/structures/")

if __name__ == "__main__":
    main()
