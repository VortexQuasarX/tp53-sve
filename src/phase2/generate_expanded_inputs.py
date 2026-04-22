"""
Phase 2 - Dataset Expansion: Generate 30 AlphaFold 3 Input JSONs
Creates input files for 30 NEW scientifically-selected mutations.

Selection Criteria:
  A. Same-Position Variants (7): Test substitution-specific effects
  B. Benign/Neutral Controls (5): Validate method discrimination
  C. Non-DBD Domain Mutations (6): Cross-domain comparison
  D. Gain-of-Function (4): GOF vs LOF structural signatures
  E. Temperature-Sensitive (3): Partial-activity mutants
  F. Rare Pathogenic (5): Non-hotspot but confirmed pathogenic
"""
import json
import os
import csv

# Wild-type p53 sequence (NP_000537.3, isoform a, 393 aa)
WT_SEQUENCE = (
    "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP"
    "DEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAK"
    "SVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHE"
    "RCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNS"
    "SCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELP"
    "PGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPG"
    "GSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD"
)

AA_MAP = {
    'A': 'Alanine', 'R': 'Arginine', 'N': 'Asparagine', 'D': 'Aspartic Acid',
    'C': 'Cysteine', 'E': 'Glutamic Acid', 'Q': 'Glutamine', 'G': 'Glycine',
    'H': 'Histidine', 'I': 'Isoleucine', 'L': 'Leucine', 'K': 'Lysine',
    'M': 'Methionine', 'F': 'Phenylalanine', 'P': 'Proline', 'S': 'Serine',
    'T': 'Threonine', 'W': 'Tryptophan', 'Y': 'Tyrosine', 'V': 'Valine'
}

CRITERION_LABELS = {
    "A": "Same-Position Variant",
    "B": "Benign/Neutral Control",
    "C": "Non-DBD Domain",
    "D": "Gain-of-Function",
    "E": "Temperature-Sensitive",
    "F": "Rare Pathogenic"
}

# ========================================================
# 30 NEW MUTATIONS - 6 Scientific Criteria Categories
# ========================================================
NEW_MUTATIONS = [
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # CRITERION A: Same-Position Variants (7 mutations)
    # Purpose: Does the SPECIFIC amino acid substitution matter, or just the position?
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    {"mutation": "R175G", "wt": "R", "pos": 175, "mut": "G",
     "criterion": "A", "classification": "Likely Oncogenic",
     "rationale": "Same position as Phase1 R175H. Glycine (smallest AA, no side chain) vs Histidine - tests size effect."},
    {"mutation": "R175C", "wt": "R", "pos": 175, "mut": "C",
     "criterion": "A", "classification": "Likely Oncogenic",
     "rationale": "Same position as Phase1 R175H. Cysteine introduces thiol group - different chemistry than Histidine."},
    {"mutation": "G245D", "wt": "G", "pos": 245, "mut": "D",
     "criterion": "A", "classification": "Likely Oncogenic",
     "rationale": "Same position as Phase1 G245S. Aspartate (charged) vs Serine (neutral) - tests charge effect."},
    {"mutation": "R282Q", "wt": "R", "pos": 282, "mut": "Q",
     "criterion": "A", "classification": "Likely Oncogenic",
     "rationale": "Same position as Phase1 R282W. Glutamine (small) vs Tryptophan (largest AA) - tests size effect."},
    {"mutation": "R248L", "wt": "R", "pos": 248, "mut": "L",
     "criterion": "A", "classification": "Likely Oncogenic",
     "rationale": "Same position as Phase1 R248Q/R248W. Leucine (hydrophobic) - tests 3rd substitution at critical DNA contact."},
    {"mutation": "Y220S", "wt": "Y", "pos": 220, "mut": "S",
     "criterion": "A", "classification": "Likely Oncogenic",
     "rationale": "Same position as Phase1 Y220C. Serine vs Cysteine - both small, but different chemistry."},
    {"mutation": "R249G", "wt": "R", "pos": 249, "mut": "G",
     "criterion": "A", "classification": "Likely Oncogenic",
     "rationale": "Same position as Phase1 R249S. Glycine (no side chain) vs Serine - tests how minimal substitution compares."},

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # CRITERION B: Benign/Neutral Controls (5 mutations)
    # Purpose: Can our RMSD method DISTINGUISH pathogenic from benign?
    # If benign mutations show LOW RMSD, our method is VALIDATED.
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    {"mutation": "P72R", "wt": "P", "pos": 72, "mut": "R",
     "criterion": "B", "classification": "Benign",
     "rationale": "THE most common TP53 polymorphism. Billions of people carry it. NOT cancer-causing. KEY validation control."},
    {"mutation": "P47S", "wt": "P", "pos": 47, "mut": "S",
     "criterion": "B", "classification": "Benign",
     "rationale": "Functional polymorphism in African populations. Mild effect on apoptosis but not pathogenic."},
    {"mutation": "A189V", "wt": "A", "pos": 189, "mut": "V",
     "criterion": "B", "classification": "VUS",
     "rationale": "Variant of Uncertain Significance. Conservative substitution (Ala→Val). Tests sensitivity threshold."},
    {"mutation": "R337H", "wt": "R", "pos": 337, "mut": "H",
     "criterion": "B", "classification": "Low-Penetrance",
     "rationale": "Brazilian founder mutation. Low cancer risk, pH-dependent tetramerization defect. Borderline case."},
    {"mutation": "K132R", "wt": "K", "pos": 132, "mut": "R",
     "criterion": "B", "classification": "VUS",
     "rationale": "Ultra-conservative substitution (Lys→Arg, both positive charge). Should show minimal RMSD if method works."},

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # CRITERION C: Non-DBD Domain Mutations (6 mutations)
    # Purpose: Do mutations OUTSIDE the DNA-binding domain cause
    #          different structural impacts?
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    {"mutation": "L344R", "wt": "L", "pos": 344, "mut": "R",
     "criterion": "C", "classification": "Likely Oncogenic",
     "rationale": "Tetramerization domain (325-355). Hydrophobic→Charged disrupts oligomerization interface."},
    {"mutation": "R342P", "wt": "R", "pos": 342, "mut": "P",
     "criterion": "C", "classification": "Likely Oncogenic",
     "rationale": "Tetramerization domain. Proline is a helix-breaker, should disrupt alpha-helix at tetramer interface."},
    {"mutation": "L22F", "wt": "L", "pos": 22, "mut": "F",
     "criterion": "C", "classification": "Likely Oncogenic",
     "rationale": "Transactivation Domain (TAD, 1-61). Affects transcriptional activation, NOT DNA binding."},
    {"mutation": "W23R", "wt": "W", "pos": 23, "mut": "R",
     "criterion": "C", "classification": "Likely Oncogenic",
     "rationale": "TAD domain. W23 is ESSENTIAL for MDM2 binding. Disrupts p53 degradation regulation."},
    {"mutation": "N345S", "wt": "N", "pos": 345, "mut": "S",
     "criterion": "C", "classification": "Likely Oncogenic",
     "rationale": "Tetramerization domain. Tests conservative substitution in non-DBD domain."},
    {"mutation": "K382R", "wt": "K", "pos": 382, "mut": "R",
     "criterion": "C", "classification": "Likely Oncogenic",
     "rationale": "C-terminal regulatory domain. K382 is acetylation site - post-translational regulation region."},

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # CRITERION D: Gain-of-Function Mutations (4 mutations)
    # Purpose: GOF mutations gain NEW oncogenic activities.
    #          Do they have different structural signatures than LOF?
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    {"mutation": "R280K", "wt": "R", "pos": 280, "mut": "K",
     "criterion": "D", "classification": "Gain-of-Function",
     "rationale": "Known GOF. Promotes invasion and metastasis. Conservative change (Arg→Lys) but catastrophic effect."},
    {"mutation": "V272M", "wt": "V", "pos": 272, "mut": "M",
     "criterion": "D", "classification": "Gain-of-Function",
     "rationale": "Known GOF. Enhances cancer aggressiveness. Tests if subtle substitutions can still be detected."},
    {"mutation": "D281G", "wt": "D", "pos": 281, "mut": "G",
     "criterion": "D", "classification": "Gain-of-Function",
     "rationale": "Known GOF. Glycine adds backbone flexibility. Tests if increased flexibility = detectable RMSD."},
    {"mutation": "S241F", "wt": "S", "pos": 241, "mut": "F",
     "criterion": "D", "classification": "Gain-of-Function",
     "rationale": "GOF mutation in L3 loop. Large hydrophobic substitution near DNA contact surface."},

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # CRITERION E: Temperature-Sensitive Mutations (3 mutations)
    # Purpose: These mutants are INACTIVE at 37°C but REGAIN
    #          function at 32°C. Structurally intermediate?
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    {"mutation": "V143A", "wt": "V", "pos": 143, "mut": "A",
     "criterion": "E", "classification": "Temperature-Sensitive",
     "rationale": "THE classic temp-sensitive p53 mutant. Inactive at 37°C, active at 32°C. Structurally rescued by cooling."},
    {"mutation": "A138V", "wt": "A", "pos": 138, "mut": "V",
     "criterion": "E", "classification": "Temperature-Sensitive",
     "rationale": "Temperature-sensitive mutant. Can be reactivated. Tests if structural damage is indeed moderate."},
    {"mutation": "I195T", "wt": "I", "pos": 195, "mut": "T",
     "criterion": "E", "classification": "Temperature-Sensitive",
     "rationale": "Partial activity mutant near zinc-binding region. Ile→Thr changes hydrophobicity dramatically."},

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # CRITERION F: Rare but Confirmed Pathogenic (5 mutations)
    # Purpose: Are rare mutations less structurally disruptive
    #          than common hotspots? Frequency vs severity.
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    {"mutation": "E285K", "wt": "E", "pos": 285, "mut": "K",
     "criterion": "F", "classification": "Likely Oncogenic",
     "rationale": "Confirmed pathogenic but rare. Charge reversal (Glu-→Lys+). Outside classic hotspot."},
    {"mutation": "N239D", "wt": "N", "pos": 239, "mut": "D",
     "criterion": "F", "classification": "Likely Oncogenic",
     "rationale": "Rare pathogenic in beta-sandwich. Asn→Asp introduces charge. Non-hotspot structural mutation."},
    {"mutation": "T125M", "wt": "T", "pos": 125, "mut": "M",
     "criterion": "F", "classification": "Likely Oncogenic",
     "rationale": "DBD non-hotspot. Thr→Met changes polarity. Tests if non-hotspot DBD mutations are less damaging."},
    {"mutation": "L194R", "wt": "L", "pos": 194, "mut": "R",
     "criterion": "F", "classification": "Likely Oncogenic",
     "rationale": "Rare pathogenic near zinc-binding. Hydrophobic→charged in structural core. Should be very disruptive."},
    {"mutation": "N247D", "wt": "N", "pos": 247, "mut": "D",
     "criterion": "F", "classification": "Likely Oncogenic",
     "rationale": "Adjacent to hotspot R248. Tests if neighboring residue mutations cause comparable damage."},
]

def mutate_sequence(wt_seq, position, new_aa):
    """Apply a point mutation to the wild-type sequence."""
    seq_list = list(wt_seq)
    idx = position - 1
    original = seq_list[idx]
    seq_list[idx] = new_aa
    return "".join(seq_list), original

def create_alphafold_json(name, sequence):
    """Create AlphaFold 3 Server input JSON."""
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
        ]
    }

def main():
    output_dir = "data/alphafold_inputs_phase2"
    os.makedirs(output_dir, exist_ok=True)
    
    # Clean previous files
    for f in os.listdir(output_dir):
        if f.endswith(".json"):
            os.remove(os.path.join(output_dir, f))
    
    print("=" * 70)
    print("PHASE 2: GENERATING 30 NEW ALPHAFOLD INPUT FILES")
    print("=" * 70)
    print(f"Wild-type sequence length: {len(WT_SEQUENCE)} amino acids")
    print(f"Total new mutations: {len(NEW_MUTATIONS)}\n")
    
    csv_rows = []
    criteria_counts = {}
    errors = []
    
    for i, mut_info in enumerate(NEW_MUTATIONS, 1):
        mutation = mut_info["mutation"]
        wt_aa = mut_info["wt"]
        pos = mut_info["pos"]
        mut_aa = mut_info["mut"]
        criterion = mut_info["criterion"]
        classification = mut_info["classification"]
        rationale = mut_info["rationale"]
        
        # Verify wild-type amino acid matches
        actual_wt = WT_SEQUENCE[pos - 1]
        if actual_wt != wt_aa:
            print(f"  ⚠️  {mutation}: Expected WT={wt_aa} at pos {pos}, found {actual_wt}. ADJUSTING.")
            wt_aa = actual_wt
            mutation = f"{actual_wt}{pos}{mut_aa}"
            mut_info["mutation"] = mutation
            mut_info["wt"] = wt_aa
        
        # Generate mutant sequence
        mut_sequence, _ = mutate_sequence(WT_SEQUENCE, pos, mut_aa)
        assert mut_sequence[pos-1] == mut_aa, f"Mutation failed for {mutation}"
        assert len(mut_sequence) == len(WT_SEQUENCE), "Sequence length changed!"
        
        # Create and save JSON
        json_name = f"TP53_{mutation}"
        json_data = create_alphafold_json(json_name, mut_sequence)
        json_path = os.path.join(output_dir, f"{json_name}.json")
        with open(json_path, "w") as f:
            json.dump(json_data, f, indent=2)
        
        csv_rows.append({
            "gene": "TP53",
            "mutation": mutation,
            "wt_residue": wt_aa,
            "position": pos,
            "mut_residue": mut_aa,
            "classification": classification,
            "criterion": criterion,
            "rationale": rationale
        })
        
        criteria_counts[criterion] = criteria_counts.get(criterion, 0) + 1
        
        label = CRITERION_LABELS.get(criterion, "Unknown")
        print(f"  [{i:2d}/30] ✅ {mutation:<8} | {criterion}: {label:<22} | "
              f"{AA_MAP.get(wt_aa, wt_aa)} → {AA_MAP.get(mut_aa, mut_aa)} @ pos {pos}")
    
    # Save Phase 2 CSV
    csv_path = "data/target_mutations_phase2.csv"
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["gene", "mutation", "wt_residue", "position",
                                                "mut_residue", "classification", "criterion", "rationale"])
        writer.writeheader()
        writer.writerows(csv_rows)
    
    # Create combined CSV (Phase 1 + Phase 2)
    import pandas as pd
    phase1_df = pd.read_csv("data/target_mutations_expanded.csv")
    phase1_df["criterion"] = "Phase1"
    phase1_df["rationale"] = "COSMIC/TCGA hotspot"
    
    phase2_df = pd.DataFrame(csv_rows)
    combined_df = pd.concat([phase1_df, phase2_df], ignore_index=True)
    combined_df.to_csv("data/target_mutations_expanded.csv", index=False)
    
    print(f"\n{'=' * 70}")
    print(f"SUMMARY")
    print(f"{'=' * 70}")
    print(f"\n  Criterion Breakdown:")
    for crit in sorted(criteria_counts.keys()):
        label = CRITERION_LABELS.get(crit, "Unknown")
        print(f"    {crit}. {label:<25}: {criteria_counts[crit]} mutations")
    
    print(f"\n  Phase 1 (existing):  {len(phase1_df)} mutations (all Likely Oncogenic)")
    print(f"  Phase 2 (NEW):       {len(csv_rows)} mutations (scientifically diverse)")
    print(f"  Combined TOTAL:      {len(combined_df)} mutations")
    
    # Classification breakdown
    print(f"\n  Classification Breakdown (Phase 2):")
    for cls in sorted(set(r["classification"] for r in csv_rows)):
        count = sum(1 for r in csv_rows if r["classification"] == cls)
        print(f"    {cls:<22}: {count}")
    
    print(f"\n  Files:")
    print(f"    JSON inputs: {output_dir}/ ({len(csv_rows)} files)")
    print(f"    Phase 2 CSV: {csv_path}")
    print(f"    Combined CSV: data/target_mutations_expanded.csv")
    
    print(f"\n{'=' * 70}")
    print(f"NEXT STEPS FOR YOU:")
    print(f"{'=' * 70}")
    print(f"  1. Open https://alphafoldserver.com")
    print(f"  2. For EACH of the 30 JSON files in {output_dir}/:")
    print(f"     → Click 'New Job' → Upload the JSON → Submit")
    print(f"  3. Download each resulting .cif file")
    print(f"  4. Rename and place in data/structures/:")
    print(f"     Example: tp53_r175g.cif, tp53_p72r.cif, etc.")
    print(f"  5. Come back - I'll re-run the ENTIRE analysis pipeline")
    print(f"{'=' * 70}")

if __name__ == "__main__":
    main()
