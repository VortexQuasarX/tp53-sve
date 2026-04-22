"""
Generate JSON data files for the React dashboard.
Reads all CSVs and outputs JSON to dashboard/public/data/
"""
import json, csv, os

BASE = "C:/Users/LENOVO/.gemini/antigravity/playground/chrono-shepard"
OUT = f"{BASE}/dashboard/public/data"

def read_csv(path):
    with open(path, encoding='utf-8') as f:
        return list(csv.DictReader(f))

def main():
    os.makedirs(OUT, exist_ok=True)

    # 1. RMSD scores
    rmsd = read_csv(f"{BASE}/output/rmsd_scores_all50.csv")
    rmsd_data = [{
        "m": r["Mutation"], "rmsd": round(float(r["RMSD (Angstroms)"]), 4),
        "cls": r["Classification"], "crit": r["Criterion"],
        "pos": int(r["Position"]), "wt": r["WT_Residue"], "mut": r["Mut_Residue"]
    } for r in rmsd]
    with open(f"{OUT}/mutations.json", "w") as f:
        json.dump(rmsd_data, f)
    print(f"  mutations.json: {len(rmsd_data)} entries")

    # 2. Per-residue summary
    summary = read_csv(f"{BASE}/output/phase2/per_residue_rmsd/per_residue_summary.csv")
    sum_data = {}
    for s in summary:
        sum_data[s["Mutation"]] = {
            "mean": round(float(s["Mean_Displacement"]), 4),
            "max": round(float(s["Max_Displacement"]), 4),
            "maxRes": int(s["Max_Residue"]),
            "above5": int(s["Residues_Above_5A"]),
            "above10": int(s["Residues_Above_10A"])
        }
    with open(f"{OUT}/summary.json", "w") as f:
        json.dump(sum_data, f)
    print(f"  summary.json: {len(sum_data)} entries")

    # 3. Per-residue profiles (all 50)
    per_res = {}
    for r in rmsd:
        m = r["Mutation"]
        p = f"{BASE}/output/phase2/per_residue_rmsd/{m}_per_residue.csv"
        if os.path.exists(p):
            rows = read_csv(p)
            per_res[m] = [round(float(row["Displacement_Angstrom"]), 2) for row in rows]
    with open(f"{OUT}/perResidue.json", "w") as f:
        json.dump(per_res, f)
    print(f"  perResidue.json: {len(per_res)} profiles ({len(next(iter(per_res.values())))} residues each)")

    # 4. Tool comparison
    tool_path = f"{BASE}/output/phase2/tool_comparison.csv"
    tool_data = [{
        "m": t["Mutation"], "rmsd": round(float(t["RMSD"]), 4),
        "sev": t["RMSD_Severity"],
        "sift": float(t["SIFT_Score"]), "pp2": float(t["PolyPhen2_Score"]),
        "rRank": int(t["RMSD_Rank"]), "sRank": int(t["SIFT_Rank"]), "pRank": int(t["PP2_Rank"])
    } for t in read_csv(tool_path)]
    with open(f"{OUT}/toolComparison.json", "w") as f:
        json.dump(tool_data, f)
    print(f"  toolComparison.json: {len(tool_data)} entries")

    # 5. Domain RMSD
    domain = read_csv(f"{BASE}/output/phase2/domain_rmsd.csv")
    dom_data = [{
        "m": d["Mutation"],
        "full": round(float(d["Full Protein"]), 4),
        "nterm": round(float(d["N-Terminal (TAD+PRD)"]), 4),
        "dbd": round(float(d["DNA-Binding Domain"]), 4),
        "tetra": round(float(d["Tetramerization"]), 4),
        "cterm": round(float(d["C-Terminal"]), 4)
    } for d in domain]
    with open(f"{OUT}/domains.json", "w") as f:
        json.dump(dom_data, f)
    print(f"  domains.json: {len(dom_data)} entries")

    # 6. Enriched data (correlations)
    enr_path = f"{BASE}/output/phase2/correlation_plots/enriched_mutation_data.csv"
    enr = read_csv(enr_path)
    enr_data = [{
        "m": e["Mutation"], "rmsd": round(float(e["RMSD (Angstroms)"]), 4),
        "plddt": round(float(e["pLDDT_Change"]), 4),
        "charge": round(float(e["Charge_Change"]), 2),
        "mass": round(float(e["Mass_Change"]), 2),
        "hydro": round(float(e["Hydrophobicity_Change"]), 2),
        "dbdDist": int(e["Distance_From_DBD_Center"])
    } for e in enr]
    with open(f"{OUT}/correlations.json", "w") as f:
        json.dump(enr_data, f)
    print(f"  correlations.json: {len(enr_data)} entries")

    print(f"\n[DONE] All data files written to {OUT}")

if __name__ == "__main__":
    main()
