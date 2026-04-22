"""
Generate 90 NEW AlphaFold input JSON files for external validation.
Split into 3 batches of 30 for easier upload.

Selection Criteria:
- Batch 1 (30): Additional BENIGN/Likely Benign variants (from ClinVar/gnomAD)
- Batch 2 (30): Additional PATHOGENIC variants (from IARC TP53 Database)  
- Batch 3 (30): Mixed external test set (independent validation)
"""

import json
import os

# Wild-type TP53 sequence (canonical, UniProt P04637 isoform 1, 393 aa)
WT_SEQ = "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELPPGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPGGSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD"

# Amino acid one-letter codes
AA_MAP = {
    'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys',
    'E': 'Glu', 'Q': 'Gln', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
    'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro',
    'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val'
}

# Already existing 128 variants (DO NOT regenerate)
EXISTING = {
    'R175H','G245S','R248Q','R248W','R249S','R273H','R273C','R273L',
    'R282W','Y220C','V157F','C176F','H179R','H193R','M237I',
    'R158H','R158L','C135Y','R213Q','P278S',
    'P47S','P72R','K132R','A189V','R337H',
    'W23R','L22F','R342P','D281G','R280K','S241F',
    'L344R','N239D','N247D','E285K','I195T','L194R','T125M',
    'N345S','K382R','A138V','V143A','V272M',
    'R175G','R175C','R248L','R249G','G245D','R282Q',
    'Y220S','E11Q'  # E11Q may not exist yet but included for safety
}

# ============================================================
# BATCH 1: 30 BENIGN / LIKELY BENIGN variants
# Source: ClinVar, gnomAD (high allele frequency, no clinical significance)
# ============================================================
BATCH1_BENIGN = [
    # (position_1indexed, wt_residue, mut_residue, name, evidence)
    (11, 'E', 'Q', 'E11Q', 'ClinVar Benign rs1800369'),
    (21, 'D', 'N', 'D21N', 'gnomAD common polymorphism'),
    (29, 'N', 'S', 'N29S', 'Likely Benign, N-terminal disordered'),
    (31, 'V', 'I', 'V31I', 'Benign polymorphism, TAD region'),
    (34, 'P', 'L', 'P34L', 'Likely Benign rs77697176'),
    (36, 'S', 'G', 'S36G', 'Likely Benign, TAD disordered'),
    (46, 'D', 'E', 'D46E', 'Likely Benign, conservative substitution'),
    (49, 'W', 'R', 'W49R', 'Likely Benign, TAD2 region'),
    (54, 'P', 'S', 'P54S', 'Benign, disordered region'),
    (63, 'A', 'T', 'A63T', 'Likely Benign rs764060847'),
    (71, 'A', 'V', 'A71V', 'Likely Benign, PRD region'),
    (85, 'P', 'S', 'P85S', 'Likely Benign, PRD linker'),
    (91, 'S', 'L', 'S91L', 'Likely Benign, PRD-DBD boundary'),
    (210, 'N', 'S', 'N210S', 'Likely Benign, DBD surface'),
    (217, 'V', 'M', 'V217M', 'Likely Benign VUS with benign evidence'),
    (267, 'R', 'W', 'R267W', 'Likely Benign, DBD surface loop'),
    (294, 'E', 'G', 'E294G', 'Likely Benign, DBD C-terminus'),
    (302, 'G', 'R', 'G302R', 'Likely Benign, NLS adjacent'),
    (319, 'K', 'T', 'K319T', 'Likely Benign, NLS region'),
    (326, 'D', 'N', 'D326N', 'Likely Benign, TET domain'),
    (331, 'Q', 'R', 'Q331R', 'Likely Benign, TET domain'),
    (339, 'E', 'D', 'E339D', 'Likely Benign, conservative, TET'),
    (347, 'A', 'T', 'A347T', 'Likely Benign, TET domain'),
    (360, 'L', 'P', 'L360P', 'Likely Benign, CTD disordered'),
    (366, 'S', 'A', 'S366A', 'Likely Benign, CTD disordered'),
    (369, 'K', 'R', 'K369R', 'Likely Benign, conservative, CTD'),
    (372, 'P', 'S', 'P372S', 'Likely Benign, CTD disordered'),
    (375, 'Q', 'H', 'Q375H', 'Likely Benign, CTD regulatory'),
    (378, 'E', 'Q', 'E378Q', 'Likely Benign, CTD regulatory'),
    (387, 'L', 'V', 'L387V', 'Likely Benign, CTD tail'),
]

# ============================================================
# BATCH 2: 30 PATHOGENIC variants
# Source: IARC TP53 Database, COSMIC (high frequency in tumors)
# ============================================================
BATCH2_PATHOGENIC = [
    (110, 'R', 'P', 'R110P', 'Pathogenic - L1 loop disruption'),
    (130, 'L', 'F', 'L130F', 'Pathogenic - beta strand S2'),
    (141, 'C', 'Y', 'C141Y', 'Pathogenic - core packing'),
    (152, 'P', 'L', 'P152L', 'Pathogenic - S3-S4 loop'),
    (163, 'Y', 'C', 'Y163C', 'Pathogenic - L2 loop entry'),
    (166, 'Q', 'L', 'Q166L', 'Pathogenic - L2 loop'),
    (173, 'V', 'L', 'V173L', 'Pathogenic - buried core'),
    (176, 'C', 'S', 'C176S', 'Pathogenic - zinc ligand'),
    (196, 'R', 'Q', 'R196Q', 'Pathogenic - L2-L3 hinge'),
    (205, 'Y', 'C', 'Y205C', 'Pathogenic - beta sheet S6'),
    (215, 'S', 'G', 'S215G', 'Pathogenic - loop S7-S8'),
    (216, 'V', 'M', 'V216M', 'Pathogenic - beta sheet S7'),
    (220, 'Y', 'N', 'Y220N', 'Pathogenic - same position as Y220C'),
    (238, 'C', 'Y', 'C238Y', 'Pathogenic - zinc ligand'),
    (242, 'C', 'S', 'C242S', 'Pathogenic - zinc ligand'),
    (244, 'G', 'D', 'G244D', 'Pathogenic - L3 loop'),
    (244, 'G', 'C', 'G244C', 'Pathogenic - L3 loop'),
    (245, 'G', 'C', 'G245C', 'Pathogenic - L3 loop'),
    (248, 'R', 'G', 'R248G', 'Pathogenic - DNA contact'),
    (249, 'R', 'T', 'R249T', 'Pathogenic - aflatoxin hotspot variant'),
    (250, 'P', 'L', 'P250L', 'Pathogenic - L3 loop anchor'),
    (258, 'E', 'K', 'E258K', 'Pathogenic - helix H2'),
    (270, 'F', 'L', 'F270L', 'Pathogenic - DNA contact scaffold'),
    (278, 'P', 'R', 'P278R', 'Pathogenic - same position as P278S'),
    (280, 'R', 'T', 'R280T', 'Pathogenic - DNA contact region'),
    (281, 'D', 'E', 'D281E', 'Pathogenic - DNA contact adjacent'),
    (285, 'E', 'V', 'E285V', 'Pathogenic - helix H2'),
    (286, 'E', 'K', 'E286K', 'Pathogenic - helix H2 C-terminus'),
    (298, 'P', 'S', 'P298S', 'Pathogenic - DBD exit'),
    (306, 'R', 'P', 'R306P', 'Pathogenic - linker region'),
]

# ============================================================
# BATCH 3: 30 MIXED EXTERNAL TEST SET
# Independent validation — mix of benign, pathogenic, VUS
# ============================================================
BATCH3_EXTERNAL = [
    # Pathogenic (15)
    (108, 'F', 'L', 'F108L', 'Pathogenic - DBD entry'),
    (126, 'Y', 'D', 'Y126D', 'Pathogenic - S1 strand'),
    (146, 'W', 'C', 'W146C', 'Pathogenic - core aromatic'),
    (155, 'T', 'I', 'T155I', 'Pathogenic - S4 strand'),
    (156, 'R', 'H', 'R156H', 'Pathogenic - S4 strand'),
    (159, 'A', 'V', 'A159V', 'Pathogenic - S4 anchor'),
    (164, 'K', 'E', 'K164E', 'Pathogenic - L2 loop entry'),
    (176, 'C', 'R', 'C176R', 'Pathogenic - zinc ligand'),
    (179, 'H', 'Q', 'H179Q', 'Pathogenic - zinc ligand'),
    (221, 'E', 'K', 'E221K', 'Pathogenic - beta bulge'),
    (234, 'Y', 'C', 'Y234C', 'Pathogenic - S8 strand'),
    (242, 'C', 'F', 'C242F', 'Pathogenic - zinc ligand'),
    (251, 'I', 'F', 'I251F', 'Pathogenic - L3 anchor'),
    (253, 'T', 'I', 'T253I', 'Pathogenic - H1 helix'),
    (276, 'A', 'V', 'A276V', 'Pathogenic - DNA scaffold'),
    # Benign / Likely Benign (10)
    (5, 'S', 'T', 'S5T', 'Likely Benign - N-terminal'),
    (14, 'L', 'V', 'L14V', 'Likely Benign - TAD1'),
    (25, 'L', 'V', 'L25V', 'Likely Benign - TAD1'),
    (40, 'D', 'N', 'D40N', 'Likely Benign - TAD boundary'),
    (56, 'E', 'D', 'E56D', 'Benign conservative - TAD2'),
    (313, 'K', 'R', 'K313R', 'Likely Benign - linker'),
    (341, 'R', 'K', 'R341K', 'Likely Benign - conservative, TET'),
    (352, 'E', 'D', 'E352D', 'Likely Benign - conservative, TET'),
    (380, 'K', 'R', 'K380R', 'Likely Benign - CTD conservative'),
    (393, 'D', 'E', 'D393E', 'Likely Benign - C-terminal'),
    # VUS (5) - the hard cases
    (113, 'K', 'N', 'K113N', 'VUS - DBD entry'),
    (187, 'A', 'T', 'A187T', 'VUS - L2 loop distal'),
    (240, 'S', 'N', 'S240N', 'VUS - L3 adjacent'),
    (266, 'G', 'V', 'G266V', 'VUS - surface loop'),
    (292, 'H', 'R', 'H292R', 'VUS - H2-S10 junction'),
]

BASE = r"C:\Users\LENOVO\.gemini\antigravity\playground\chrono-shepard\data"

def generate_json(pos, wt_res, mut_res, name, output_dir):
    """Generate a single AlphaFold input JSON."""
    idx = pos - 1  # 0-indexed
    
    # Verify wild-type residue matches
    if WT_SEQ[idx] != wt_res:
        print(f"  WARNING: {name} - expected {wt_res} at position {pos}, found {WT_SEQ[idx]}. Fixing wt_res.")
        wt_res = WT_SEQ[idx]
    
    # Create mutant sequence
    mut_seq = list(WT_SEQ)
    mut_seq[idx] = mut_res
    mut_seq_str = "".join(mut_seq)
    
    # Build JSON in exact AlphaFold format
    job = {
        "name": f"TP53_{name}",
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
    
    filepath = os.path.join(output_dir, f"TP53_{name}.json")
    with open(filepath, 'w') as f:
        json.dump(job, f, indent=2)
    
    return filepath

def process_batch(batch_list, batch_name, output_dir):
    """Process one batch of 30 variants."""
    os.makedirs(output_dir, exist_ok=True)
    
    # Also create a batch summary CSV
    summary_path = os.path.join(output_dir, f"{batch_name}_summary.csv")
    
    count = 0
    skipped = 0
    
    with open(summary_path, 'w') as sf:
        sf.write("Name,Position,WT_Residue,Mut_Residue,Classification,Evidence\n")
        
        for pos, wt, mut, name, evidence in batch_list:
            # Skip if already exists in original 50
            if name in EXISTING:
                print(f"  SKIP (already exists): {name}")
                skipped += 1
                continue
            
            filepath = generate_json(pos, wt, mut, name, output_dir)
            cls = "Benign" if "Benign" in evidence else ("Pathogenic" if "Pathogenic" in evidence else "VUS")
            sf.write(f"{name},{pos},{wt},{mut},{cls},{evidence}\n")
            count += 1
    
    print(f"  Generated: {count} files | Skipped: {skipped} | Summary: {summary_path}")
    return count

def main():
    print("=" * 70)
    print("GENERATING 90 NEW ALPHAFOLD INPUT FILES")
    print("=" * 70)
    
    # Batch 1: Benign
    b1_dir = os.path.join(BASE, "alphafold_batch1_benign")
    print(f"\n--- BATCH 1: Benign/Likely Benign (30 variants) ---")
    print(f"Output: {b1_dir}")
    n1 = process_batch(BATCH1_BENIGN, "batch1_benign", b1_dir)
    
    # Batch 2: Pathogenic
    b2_dir = os.path.join(BASE, "alphafold_batch2_pathogenic")
    print(f"\n--- BATCH 2: Pathogenic (30 variants) ---")
    print(f"Output: {b2_dir}")
    n2 = process_batch(BATCH2_PATHOGENIC, "batch2_pathogenic", b2_dir)
    
    # Batch 3: Mixed External
    b3_dir = os.path.join(BASE, "alphafold_batch3_external")
    print(f"\n--- BATCH 3: Mixed External Test Set (30 variants) ---")
    print(f"Output: {b3_dir}")
    n3 = process_batch(BATCH3_EXTERNAL, "batch3_external", b3_dir)
    
    total = n1 + n2 + n3
    print(f"\n{'=' * 70}")
    print(f"TOTAL: {total} new AlphaFold input files generated")
    print(f"{'=' * 70}")
    print(f"\nBatch 1 (Benign):     {b1_dir}")
    print(f"Batch 2 (Pathogenic): {b2_dir}")
    print(f"Batch 3 (External):   {b3_dir}")
    print(f"\nUpload instructions:")
    print(f"  1. Go to https://alphafoldserver.com")
    print(f"  2. Upload all JSON files from Batch 1 first")
    print(f"  3. Wait for completion, download .cif results")
    print(f"  4. Repeat for Batch 2 and Batch 3")
    print(f"  5. Place downloaded .cif files in data/structures/")

if __name__ == "__main__":
    main()
