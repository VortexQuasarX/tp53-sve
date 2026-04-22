"""
Generate Ultimate Interactive Dashboard
Reads ALL analysis data and embeds it into a premium HTML dashboard.
Uses string template replacement (not f-strings) to avoid JS brace conflicts.
"""
import json, csv, os

BASE = "C:/Users/LENOVO/.gemini/antigravity/playground/chrono-shepard"

def read_csv(path):
    with open(path, encoding='utf-8') as f:
        return list(csv.DictReader(f))

def load_per_residue(mutation):
    p = f"{BASE}/output/phase2/per_residue_rmsd/{mutation}_per_residue.csv"
    if not os.path.exists(p): return []
    rows = read_csv(p)
    return [{"r":int(r["Residue_Number"]),"d":round(float(r["Displacement_Angstrom"]),2)} for r in rows]

def main():
    rmsd = read_csv(f"{BASE}/output/rmsd_scores.csv")
    summary = read_csv(f"{BASE}/output/phase2/per_residue_rmsd/per_residue_summary.csv")
    tool_cmp = read_csv(f"{BASE}/output/phase2/tool_comparison.csv")
    domain = read_csv(f"{BASE}/output/phase2/domain_rmsd.csv")
    
    rmsd_data = []
    for r in rmsd:
        rmsd_data.append({
            "m": r["Mutation"], "rmsd": float(r["RMSD (Angstroms)"]),
            "cls": r["Classification"], "crit": r["Criterion"],
            "pos": int(r["Position"]), "wt": r["WT_Residue"], "mut": r["Mut_Residue"]
        })
    
    sum_data = {}
    for s in summary:
        sum_data[s["Mutation"]] = {
            "mean": float(s["Mean_Displacement"]), "max": float(s["Max_Displacement"]),
            "maxRes": int(s["Max_Residue"]),
            "above5": int(s["Residues_Above_5A"]), "above10": int(s["Residues_Above_10A"])
        }
    
    per_res = {}
    for r in rmsd:
        data = load_per_residue(r["Mutation"])
        if data:
            per_res[r["Mutation"]] = data
    
    tool_data = []
    for t in tool_cmp:
        tool_data.append({
            "m": t["Mutation"], "rmsd": float(t["RMSD"]),
            "sift": float(t["SIFT_Score"]), "pp2": float(t["PolyPhen2_Score"]),
            "rRank": int(t["RMSD_Rank"])
        })
    
    dom_data = []
    for d in domain:
        dom_data.append({
            "m": d["Mutation"], "full": float(d["Full Protein"]),
            "nterm": float(d["N-Terminal (TAD+PRD)"]),
            "dbd": float(d["DNA-Binding Domain"]),
            "tetra": float(d["Tetramerization"]),
            "cterm": float(d["C-Terminal"])
        })

    print(f"  Loaded {len(rmsd_data)} mutations")
    print(f"  Loaded {len(per_res)} per-residue profiles")
    print(f"  Loaded {len(tool_data)} tool comparisons")
    print(f"  Loaded {len(dom_data)} domain records")

    # Read template
    tpl_path = f"{BASE}/src/phase2/dashboard_template.html"
    with open(tpl_path, 'r', encoding='utf-8') as f:
        html = f.read()
    
    # Replace data placeholders
    html = html.replace('/*__DATA__*/', json.dumps(rmsd_data))
    html = html.replace('/*__SUM__*/', json.dumps(sum_data))
    html = html.replace('/*__PR__*/', json.dumps(per_res))
    html = html.replace('/*__TOOL__*/', json.dumps(tool_data))
    html = html.replace('/*__DOM__*/', json.dumps(dom_data))
    
    out_path = f"{BASE}/output/phase2/dashboard.html"
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(html)
    
    print(f"\n[SUCCESS] Dashboard: {out_path}")
    print(f"  Size: {os.path.getsize(out_path)/1024:.0f} KB")

if __name__ == "__main__":
    main()
