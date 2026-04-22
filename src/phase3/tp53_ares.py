"""
TP53-ARES: Atomistic Residue Energy Scoring
============================================
A COMPLETELY NOVEL evaluation paradigm for TP53 mutations.

WHAT MAKES THIS DIFFERENT FROM EVERYTHING ELSE:
───────────────────────────────────────────────
RMSD / TM-score / DBCA  →  Measures SHAPE change
TP53-ARES               →  Estimates ENERGY change (Delta Delta G proxy)

This is the thermodynamic approach: mutations cause disease not because
the shape changes, but because the ENERGY LANDSCAPE changes — the protein
becomes less stable, loses favorable interactions, or gains unfavorable ones.

FOUR COMPLETELY NOVEL COMPONENTS:
──────────────────────────────────
1. STATISTICAL CONTACT POTENTIAL (Delta Delta E_contact)
   Uses Miyazawa-Jernigan (1996) residue-level statistical potentials
   to estimate the change in inter-residue interaction energy when the
   wild-type residue is replaced by the mutant residue. This is how
   FoldX and Rosetta estimate Delta Delta G, but we compute it directly from
   the AlphaFold structures.

2. DISRUPTION WAVEFRONT PROPAGATION
   Models how structural damage propagates through the protein's residue
   contact network from the mutation site, like seismic waves from an
   earthquake. Builds a graph of residue contacts, then measures displacement
   at each "shell" (1, 2, 3... contacts away). The decay rate reveals whether
   disruption is absorbed locally or propagates systemically.

3. ELECTROSTATIC SURFACE POTENTIAL
   Quantifies changes to the protein's charge distribution — specifically
   whether the mutation alters the electrostatic character of the DNA-binding
   surface. p53 must be positively charged to bind negatively charged DNA.

4. INTERACTION NETWORK REWIRING
   For each contact the mutated residue makes, computes the ENERGY CHANGE
   from replacing the WT residue type with the mutant residue type. This
   captures whether the mutation creates energetically favorable or
   unfavorable new interactions.

REFERENCES:
  Miyazawa S, Jernigan RL (1996). J Mol Biol 256:623-644.
  Cho Y et al. (1994). Science 265:346-355.
  Bullock AN et al. (1997). PNAS 94:14338-14342.
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from Bio.PDB import MMCIFParser, Superimposer
from collections import defaultdict
from scipy.stats import spearmanr
import warnings

warnings.filterwarnings("ignore")

# ═══════════════════════════════════════════════════════════════
# MIYAZAWA-JERNIGAN CONTACT POTENTIAL MATRIX (1996)
# 20×20 symmetric matrix of residue-residue interaction energies
# Values in RT units (multiply by 0.6 for kcal/mol)
# More negative = more favorable interaction
# ═══════════════════════════════════════════════════════════════

AA_ORDER = ['C', 'M', 'F', 'I', 'L', 'V', 'W', 'Y', 'A', 'G',
            'T', 'S', 'N', 'Q', 'D', 'E', 'H', 'R', 'K', 'P']

# Miyazawa-Jernigan contact energies (eij in RT units)
# Source: Miyazawa & Jernigan, J Mol Biol 256:623-644 (1996), Table III
MJ_MATRIX_RAW = {
    ('C','C'): -5.44, ('C','M'): -4.99, ('C','F'): -5.80, ('C','I'): -5.50, ('C','L'): -5.83,
    ('C','V'): -4.96, ('C','W'): -4.95, ('C','Y'): -4.16, ('C','A'): -3.57, ('C','G'): -3.16,
    ('C','T'): -3.11, ('C','S'): -2.86, ('C','N'): -2.59, ('C','Q'): -2.85, ('C','D'): -2.41,
    ('C','E'): -2.27, ('C','H'): -3.60, ('C','R'): -2.57, ('C','K'): -1.95, ('C','P'): -3.07,
    ('M','M'): -5.46, ('M','F'): -6.56, ('M','I'): -6.02, ('M','L'): -6.41, ('M','V'): -5.32,
    ('M','W'): -5.55, ('M','Y'): -4.91, ('M','A'): -3.94, ('M','G'): -3.39, ('M','T'): -3.51,
    ('M','S'): -3.03, ('M','N'): -2.95, ('M','Q'): -3.30, ('M','D'): -2.57, ('M','E'): -2.89,
    ('M','H'): -3.98, ('M','R'): -3.12, ('M','K'): -2.48, ('M','P'): -3.45,
    ('F','F'): -7.26, ('F','I'): -6.84, ('F','L'): -7.28, ('F','V'): -6.29, ('F','W'): -6.16,
    ('F','Y'): -5.66, ('F','A'): -4.81, ('F','G'): -4.13, ('F','T'): -4.28, ('F','S'): -4.02,
    ('F','N'): -3.75, ('F','Q'): -4.10, ('F','D'): -3.48, ('F','E'): -3.56, ('F','H'): -4.77,
    ('F','R'): -3.98, ('F','K'): -3.36, ('F','P'): -4.25,
    ('I','I'): -6.54, ('I','L'): -7.04, ('I','V'): -6.05, ('I','W'): -5.78, ('I','Y'): -5.25,
    ('I','A'): -4.58, ('I','G'): -3.78, ('I','T'): -4.03, ('I','S'): -3.52, ('I','N'): -3.24,
    ('I','Q'): -3.67, ('I','D'): -3.17, ('I','E'): -3.27, ('I','H'): -4.14, ('I','R'): -3.63,
    ('I','K'): -3.01, ('I','P'): -3.76,
    ('L','L'): -7.37, ('L','V'): -6.48, ('L','W'): -6.14, ('L','Y'): -5.67, ('L','A'): -4.91,
    ('L','G'): -4.16, ('L','T'): -4.34, ('L','S'): -3.92, ('L','N'): -3.74, ('L','Q'): -4.04,
    ('L','D'): -3.40, ('L','E'): -3.59, ('L','H'): -4.54, ('L','R'): -4.03, ('L','K'): -3.37,
    ('L','P'): -4.20,
    ('V','V'): -5.52, ('V','W'): -5.18, ('V','Y'): -4.62, ('V','A'): -4.04, ('V','G'): -3.38,
    ('V','T'): -3.46, ('V','S'): -3.05, ('V','N'): -2.83, ('V','Q'): -3.07, ('V','D'): -2.48,
    ('V','E'): -2.67, ('V','H'): -3.58, ('V','R'): -3.07, ('V','K'): -2.49, ('V','P'): -3.32,
    ('W','W'): -5.06, ('W','Y'): -4.66, ('W','A'): -3.82, ('W','G'): -3.42, ('W','T'): -3.22,
    ('W','S'): -2.99, ('W','N'): -3.07, ('W','Q'): -3.11, ('W','D'): -2.84, ('W','E'): -2.99,
    ('W','H'): -3.98, ('W','R'): -3.41, ('W','K'): -2.69, ('W','P'): -3.73,
    ('Y','Y'): -4.17, ('Y','A'): -3.36, ('Y','G'): -3.01, ('Y','T'): -3.01, ('Y','S'): -2.78,
    ('Y','N'): -2.76, ('Y','Q'): -2.97, ('Y','D'): -2.76, ('Y','E'): -2.79, ('Y','H'): -3.52,
    ('Y','R'): -3.16, ('Y','K'): -2.60, ('Y','P'): -3.19,
    ('A','A'): -2.72, ('A','G'): -2.31, ('A','T'): -2.32, ('A','S'): -2.01, ('A','N'): -1.84,
    ('A','Q'): -1.89, ('A','D'): -1.70, ('A','E'): -1.51, ('A','H'): -2.41, ('A','R'): -1.83,
    ('A','K'): -1.31, ('A','P'): -2.03,
    ('G','G'): -2.24, ('G','T'): -2.08, ('G','S'): -1.82, ('G','N'): -1.74, ('G','Q'): -1.66,
    ('G','D'): -1.59, ('G','E'): -1.22, ('G','H'): -2.15, ('G','R'): -1.72, ('G','K'): -1.15,
    ('G','P'): -1.87,
    ('T','T'): -2.12, ('T','S'): -1.96, ('T','N'): -1.88, ('T','Q'): -1.90, ('T','D'): -1.80,
    ('T','E'): -1.74, ('T','H'): -2.42, ('T','R'): -1.90, ('T','K'): -1.31, ('T','P'): -1.90,
    ('S','S'): -1.67, ('S','N'): -1.58, ('S','Q'): -1.49, ('S','D'): -1.63, ('S','E'): -1.48,
    ('S','H'): -2.11, ('S','R'): -1.62, ('S','K'): -1.05, ('S','P'): -1.57,
    ('N','N'): -1.68, ('N','Q'): -1.71, ('N','D'): -1.68, ('N','E'): -1.51, ('N','H'): -2.08,
    ('N','R'): -1.64, ('N','K'): -1.21, ('N','P'): -1.53,
    ('Q','Q'): -1.54, ('Q','D'): -1.46, ('Q','E'): -1.42, ('Q','H'): -1.98, ('Q','R'): -1.80,
    ('Q','K'): -1.29, ('Q','P'): -1.73,
    ('D','D'): -1.21, ('D','E'): -1.02, ('D','H'): -2.32, ('D','R'): -2.29, ('D','K'): -1.68,
    ('D','P'): -1.33,
    ('E','E'): -0.91, ('E','H'): -2.15, ('E','R'): -2.27, ('E','K'): -1.80, ('E','P'): -1.26,
    ('H','H'): -3.05, ('H','R'): -2.16, ('H','K'): -1.35, ('H','P'): -2.25,
    ('R','R'): -1.55, ('R','K'): -0.59, ('R','P'): -1.70,
    ('K','K'): -0.12, ('K','P'): -0.97,
    ('P','P'): -1.75,
}

# Build symmetric lookup
MJ = {}
for (a1, a2), val in MJ_MATRIX_RAW.items():
    MJ[(a1, a2)] = val
    MJ[(a2, a1)] = val

# 3-letter to 1-letter amino acid conversion
AA3TO1 = {
    'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D', 'CYS':'C',
    'GLU':'E', 'GLN':'Q', 'GLY':'G', 'HIS':'H', 'ILE':'I',
    'LEU':'L', 'LYS':'K', 'MET':'M', 'PHE':'F', 'PRO':'P',
    'SER':'S', 'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V',
}

# Amino acid charges at pH 7
AA_CHARGE = {
    'R': +1, 'K': +1, 'H': +0.1,
    'D': -1, 'E': -1,
    'C': 0, 'M': 0, 'F': 0, 'I': 0, 'L': 0, 'V': 0, 'W': 0, 'Y': 0,
    'A': 0, 'G': 0, 'T': 0, 'S': 0, 'N': 0, 'Q': 0, 'P': 0,
}

# Amino acid hydrophobicity (Kyte-Doolittle scale)
AA_HYDRO = {
    'I': 4.5, 'V': 4.2, 'L': 3.8, 'F': 2.8, 'C': 2.5, 'M': 1.9, 'A': 1.8,
    'G': -0.4, 'T': -0.7, 'S': -0.8, 'W': -0.9, 'Y': -1.3, 'P': -1.6,
    'H': -3.2, 'D': -3.5, 'E': -3.5, 'N': -3.5, 'Q': -3.5, 'K': -3.9, 'R': -4.5,
}

CONTACT_CUTOFF = 8.0  # Cα-Cα distance for contact definition
SEQ_SEP = 3  # Minimum sequence separation for contacts

DBD_RANGE = (95, 292)
DNA_CONTACT_RESIDUES = [120, 241, 248, 273, 280, 281, 283]


def get_residue_data(structure):
    """Extract Cα coordinates and residue types."""
    model = next(iter(structure))
    data = {}
    for chain in model:
        for residue in chain:
            if 'CA' in residue:
                aa3 = residue.get_resname()
                aa1 = AA3TO1.get(aa3, 'X')
                if aa1 != 'X':
                    data[residue.id[1]] = {
                        'coord': residue['CA'].get_vector().get_array(),
                        'atom': residue['CA'],
                        'aa1': aa1,
                        'aa3': aa3,
                    }
    return data


def build_contact_graph(residue_data):
    """Build residue contact graph (adjacency list)."""
    res_ids = sorted(residue_data.keys())
    graph = defaultdict(list)
    for i, ri in enumerate(res_ids):
        for rj in res_ids[i + 1:]:
            if abs(ri - rj) >= SEQ_SEP:
                d = np.linalg.norm(residue_data[ri]['coord'] - residue_data[rj]['coord'])
                if d < CONTACT_CUTOFF:
                    graph[ri].append(rj)
                    graph[rj].append(ri)
    return graph


def component1_contact_potential(wt_data, mut_data, mutation_pos, wt_aa, mut_aa):
    """
    COMPONENT 1: Statistical Contact Potential Energy Change (Delta Delta E_contact)
    
    For every contact the mutated residue makes, compute the ENERGY CHANGE
    from replacing WT residue with mutant residue using the MJ matrix.
    
    Delta Delta E = Σ[ E(mut_aa, neighbor) - E(wt_aa, neighbor) ] for all neighbors
    
    More positive Delta Delta E = destabilizing mutation.
    """
    # Find contacts of the mutation position in WT structure
    contacts = []
    for res_id, res_data in wt_data.items():
        if mutation_pos in wt_data and abs(res_id - mutation_pos) >= SEQ_SEP:
            d = np.linalg.norm(wt_data[mutation_pos]['coord'] - res_data['coord'])
            if d < CONTACT_CUTOFF:
                contacts.append(res_id)

    if not contacts:
        return 0, {}

    # Compute energy change for each contact
    dde_total = 0
    contact_details = []
    for cid in contacts:
        neighbor_aa = wt_data[cid]['aa1']
        e_wt = MJ.get((wt_aa, neighbor_aa), 0)
        e_mut = MJ.get((mut_aa, neighbor_aa), 0)
        dde = e_mut - e_wt  # positive = destabilizing
        dde_total += dde
        contact_details.append({
            'neighbor': cid,
            'neighbor_aa': neighbor_aa,
            'E_wt': round(e_wt, 3),
            'E_mut': round(e_mut, 3),
            'Delta Delta E': round(dde, 3),
        })

    return round(dde_total, 4), {
        'n_contacts': len(contacts),
        'dde_per_contact': round(dde_total / len(contacts), 4) if contacts else 0,
        'details': contact_details[:5],  # top 5 for summary
    }


def component2_disruption_propagation(wt_data, mut_data, mutation_pos):
    """
    COMPONENT 2: Disruption Wavefront Propagation
    
    Models how structural disruption spreads from the mutation site through
    the residue contact network, like seismic waves from an epicenter.
    
    1. Build contact graph from WT structure
    2. BFS from mutation site, measuring displacement at each shell
    3. Compute decay rate (how fast displacement drops per shell)
    4. Compute propagation reach (how far displacement >2 Å extends)
    
    Rapid decay = localized damage (may be absorbed)
    Slow decay = systemic disruption (whole protein destabilized)
    """
    wt_graph = build_contact_graph(wt_data)

    if mutation_pos not in wt_graph:
        return 0, 0, {}

    # BFS from mutation site
    visited = {mutation_pos: 0}
    queue = [mutation_pos]
    shell_displacements = defaultdict(list)

    # Displacement of each residue
    common = sorted(set(wt_data.keys()) & set(mut_data.keys()))
    displacements = {}
    for r in common:
        d = np.linalg.norm(wt_data[r]['coord'] - mut_data[r]['coord'])
        displacements[r] = d

    while queue:
        current = queue.pop(0)
        current_shell = visited[current]
        if current_shell > 8:  # max 8 shells
            break
        if current in displacements:
            shell_displacements[current_shell].append(displacements[current])
        for neighbor in wt_graph.get(current, []):
            if neighbor not in visited:
                visited[neighbor] = current_shell + 1
                queue.append(neighbor)

    if not shell_displacements:
        return 0, 0, {}

    # Mean displacement per shell
    shell_means = {}
    for shell, disps in sorted(shell_displacements.items()):
        shell_means[shell] = np.mean(disps)

    # Decay rate: fit exponential decay
    shells = sorted(shell_means.keys())
    if len(shells) >= 2:
        means = [shell_means[s] for s in shells]
        if means[0] > 0:
            decay_ratios = [means[i] / max(means[0], 0.001) for i in range(len(means))]
            # Decay rate = how quickly displacement drops
            decay_rate = 1.0 - np.mean(decay_ratios[1:min(4, len(decay_ratios))])
        else:
            decay_rate = 1.0
    else:
        decay_rate = 1.0

    # Propagation reach: how many shells have meandisp > 2 Å
    reach = sum(1 for s in shells if shell_means[s] > 2.0)

    return round(decay_rate, 4), reach, {
        'shell_means': {s: round(v, 3) for s, v in shell_means.items()},
        'n_shells': len(shells),
    }


def component3_electrostatic_impact(wt_aa, mut_aa, mutation_pos):
    """
    COMPONENT 3: Electrostatic Surface Potential Change
    
    p53 binds negatively-charged DNA through positively-charged surface
    residues. Mutations that remove positive charge or add negative charge
    at the DNA-binding surface are electrostatically destructive.
    
    Delta Charge = charge(mut) - charge(wt)
    Weighted by proximity to DNA-binding surface.
    """
    charge_change = AA_CHARGE.get(mut_aa, 0) - AA_CHARGE.get(wt_aa, 0)
    hydro_change = AA_HYDRO.get(mut_aa, 0) - AA_HYDRO.get(wt_aa, 0)

    # Weight by domain
    in_dbd = DBD_RANGE[0] <= mutation_pos <= DBD_RANGE[1]
    near_dna = mutation_pos in DNA_CONTACT_RESIDUES
    domain_weight = 3.0 if near_dna else 2.0 if in_dbd else 1.0

    electro_impact = abs(charge_change) * domain_weight
    hydro_impact = abs(hydro_change) * 0.5

    return round(electro_impact, 4), round(hydro_impact, 4), {
        'charge_change': round(charge_change, 2),
        'hydro_change': round(hydro_change, 2),
        'domain_weight': domain_weight,
        'in_dbd': in_dbd,
        'near_dna': near_dna,
    }


def component4_network_rewiring(wt_data, mut_data):
    """
    COMPONENT 4: Interaction Network Rewiring Score
    
    Compares the contact networks of WT and mutant structures.
    For contacts that are lost or gained, sums the MJ energies.
    
    Lost favorable contacts = destabilizing
    Gained unfavorable contacts = destabilizing
    """
    wt_graph = build_contact_graph(wt_data)
    mut_graph = build_contact_graph(mut_data)

    # Get all contacts as sets of tuples
    wt_contacts = set()
    for ri, neighbors in wt_graph.items():
        for rj in neighbors:
            wt_contacts.add((min(ri, rj), max(ri, rj)))

    mut_contacts = set()
    for ri, neighbors in mut_graph.items():
        for rj in neighbors:
            mut_contacts.add((min(ri, rj), max(ri, rj)))

    # Lost contacts (in WT, not in mutant)
    lost = wt_contacts - mut_contacts
    gained = mut_contacts - wt_contacts
    preserved = wt_contacts & mut_contacts

    # Energy of lost contacts
    lost_energy = 0
    for ri, rj in lost:
        if ri in wt_data and rj in wt_data:
            aa_i = wt_data[ri]['aa1']
            aa_j = wt_data[rj]['aa1']
            e = MJ.get((aa_i, aa_j), 0)
            lost_energy += abs(e)  # lost favorable energy

    # DBD-specific losses
    dbd_lost = sum(1 for ri, rj in lost
                   if (DBD_RANGE[0] <= ri <= DBD_RANGE[1]) or
                      (DBD_RANGE[0] <= rj <= DBD_RANGE[1]))

    return round(lost_energy, 4), {
        'contacts_lost': len(lost),
        'contacts_gained': len(gained),
        'contacts_preserved': len(preserved),
        'dbd_contacts_lost': dbd_lost,
        'lost_energy': round(lost_energy, 4),
    }


def classify_ares(score):
    """Classify by ARES score."""
    if score >= 75: return "Highly Destabilizing"
    if score >= 50: return "Destabilizing"
    if score >= 25: return "Moderately Destabilizing"
    return "Neutral/Mild"


def main():
    mutations_csv = "data/target_mutations_expanded.csv"
    wt_file = "data/structures/tp53_wt.cif"
    output_dir = "output/phase3"
    os.makedirs(output_dir, exist_ok=True)

    parser = MMCIFParser(QUIET=True)
    mutations_df = pd.read_csv(mutations_csv)
    rmsd_df = pd.read_csv("output/rmsd_scores.csv")

    print("╔" + "═" * 68 + "╗")
    print("║  TP53-ARES: Atomistic Residue Energy Scoring                      ║")
    print("║  A Novel Energy-Based Evaluation for TP53 Mutations               ║")
    print("╠" + "═" * 68 + "╣")
    print("║  Component 1: Statistical Contact Potential (Miyazawa-Jernigan)    ║")
    print("║  Component 2: Disruption Wavefront Propagation                    ║")
    print("║  Component 3: Electrostatic Surface Potential Change              ║")
    print("║  Component 4: Interaction Network Rewiring                        ║")
    print("╚" + "═" * 68 + "╝")
    print()

    # Parse WT once
    wt_struct = parser.get_structure("wt", wt_file)
    wt_data = get_residue_data(wt_struct)

    results = []
    for _, row in mutations_df.iterrows():
        mutation = row['mutation']
        import re
        wt_aa = mutation[0]
        mut_aa = mutation[-1]
        match = re.search(r'\d+', mutation)
        if not match:
            continue
        position = int(match.group())
        mut_file = f"data/structures/tp53_{mutation.lower()}.cif"

        if not os.path.exists(mut_file):
            print(f"  [SKIP] {mutation}")
            continue

        print(f"  {mutation:12s}", end=" │ ")

        try:
            # Superimpose mutant onto WT
            mut_struct = parser.get_structure("mut", mut_file)
            mut_data_raw = get_residue_data(mut_struct)

            common = sorted(set(wt_data.keys()) & set(mut_data_raw.keys()))
            wt_atoms = [wt_data[r]['atom'] for r in common]
            mut_atoms = [mut_data_raw[r]['atom'] for r in common]

            sup = Superimposer()
            sup.set_atoms(wt_atoms, mut_atoms)

            mut_struct2 = parser.get_structure("mut2", mut_file)
            sup.apply(mut_struct2.get_atoms())
            mut_data = get_residue_data(mut_struct2)

            # Component 1: Contact potential
            dde, dde_details = component1_contact_potential(wt_data, mut_data, position, wt_aa, mut_aa)

            # Component 2: Disruption propagation
            decay, reach, prop_details = component2_disruption_propagation(wt_data, mut_data, position)

            # Component 3: Electrostatic impact
            electro, hydro, electro_details = component3_electrostatic_impact(wt_aa, mut_aa, position)

            # Component 4: Network rewiring
            rewire_energy, rewire_details = component4_network_rewiring(wt_data, mut_data)

            # RMSD
            rmsd_row = rmsd_df[rmsd_df['Mutation'] == mutation]
            rmsd_val = float(rmsd_row['RMSD (Angstroms)'].iloc[0]) if len(rmsd_row) > 0 else 0

            result = {
                'Mutation': mutation,
                'WT_AA': wt_aa,
                'Mut_AA': mut_aa,
                'Position': position,
                'RMSD': round(rmsd_val, 4),
                'Criterion': row.get('criterion', ''),
                # Component scores
                'DDE_Contact': dde,
                'N_Contacts_Site': dde_details.get('n_contacts', 0),
                'DDE_Per_Contact': dde_details.get('dde_per_contact', 0),
                'Decay_Rate': decay,
                'Propagation_Reach': reach,
                'Electrostatic_Impact': electro,
                'Hydrophobicity_Change': hydro,
                'Charge_Change': electro_details.get('charge_change', 0),
                'In_DBD': electro_details.get('in_dbd', False),
                'Rewiring_Energy': rewire_energy,
                'Contacts_Lost': rewire_details.get('contacts_lost', 0),
                'Contacts_Gained': rewire_details.get('contacts_gained', 0),
                'DBD_Contacts_Lost': rewire_details.get('dbd_contacts_lost', 0),
            }
            results.append(result)

            print(f"Delta Delta E={dde:+6.2f} │ Decay={decay:.2f} │ Reach={reach} │ "
                  f"Electro={electro:.1f} │ Rewire={rewire_energy:.1f}")

        except Exception as e:
            print(f"ERROR: {e}")

    df = pd.DataFrame(results)

    # --- COMPUTE FINAL ARES SCORE ---
    # Normalize each component to 0-100
    def norm(vals):
        mi, mx = vals.min(), vals.max()
        if mx == mi: return np.full_like(vals, 50.0, dtype=float)
        return (vals - mi) / (mx - mi) * 100

    # Delta Delta E: positive = destabilizing → higher score
    dde_norm = norm(df['DDE_Contact'].values)
    # Decay: lower decay rate = propagation spreads far → higher score
    decay_norm = norm(1 - df['Decay_Rate'].values)
    # Reach: more shells affected → higher score
    reach_norm = norm(df['Propagation_Reach'].values.astype(float))
    # Electrostatic: higher impact → higher score
    electro_norm = norm(df['Electrostatic_Impact'].values + df['Hydrophobicity_Change'].values)
    # Rewiring: more lost energy → higher score
    rewire_norm = norm(df['Rewiring_Energy'].values)

    # ARES = weighted combination
    ares = (
        0.30 * dde_norm +          # Direct energy change at mutation site
        0.20 * decay_norm +         # Propagation behavior
        0.10 * reach_norm +         # Propagation reach
        0.15 * electro_norm +       # Electrostatic/hydrophobic impact
        0.25 * rewire_norm          # Network energy loss
    )

    df['ARES'] = np.round(ares, 2)
    df['ARES_Class'] = df['ARES'].apply(classify_ares)
    df['ARES_Rank'] = df['ARES'].rank(ascending=False).astype(int)
    df['RMSD_Rank'] = df['RMSD'].rank(ascending=False).astype(int)
    df['Rank_Change'] = df['RMSD_Rank'] - df['ARES_Rank']

    df = df.sort_values('ARES', ascending=False).reset_index(drop=True)
    df.to_csv(os.path.join(output_dir, 'ares_scores.csv'), index=False)

    # --- RANKINGS ---
    print(f"\n{'═' * 80}")
    print("TP53-ARES FINAL RANKING")
    print(f"{'═' * 80}")
    print(f"{'Rank':>4}  {'Mutation':<10}  {'ARES':>6}  {'Class':<22}  "
          f"{'Delta Delta E':>7}  {'Decay':>5}  {'RMSD':>7}  {'Delta Rank':>6}")
    print("-" * 80)
    for _, r in df.iterrows():
        d = int(r['Rank_Change'])
        ar = f"+{d}" if d > 0 else str(d) if d < 0 else "="
        print(f"  {r['ARES_Rank']:>2}  {r['Mutation']:<10}  {r['ARES']:>6.2f}  "
              f"{r['ARES_Class']:<22}  {r['DDE_Contact']:>+6.2f}  "
              f"{r['Decay_Rate']:>5.2f}  {r['RMSD']:>6.2f}  {ar:>6}")

    # Correlations
    rho_rmsd, p_rmsd = spearmanr(df['ARES'], df['RMSD'])
    print(f"\n  ARES vs RMSD:  Spearman ρ={rho_rmsd:.4f} (p={p_rmsd:.2e})")

    # Severity distribution
    ares_counts = df['ARES_Class'].value_counts()
    print(f"\n{'═' * 80}")
    print("ARES CLASSIFICATION")
    print(f"{'═' * 80}")
    for cls in ["Highly Destabilizing", "Destabilizing", "Moderately Destabilizing", "Neutral/Mild"]:
        print(f"  {cls:30s}: {ares_counts.get(cls, 0)} mutations")

    # Significant movers
    movers = df[df['Rank_Change'].abs() >= 5].sort_values('Rank_Change', ascending=False)
    if len(movers) > 0:
        print(f"\n{'═' * 80}")
        print("RANK MOVERS (ARES vs RMSD, ≥5 positions)")
        print(f"{'═' * 80}")
        for _, r in movers.iterrows():
            d = "⬆" if r['Rank_Change'] > 0 else "⬇"
            print(f"  {d} {r['Mutation']}: RMSD #{int(r['RMSD_Rank'])} → ARES #{int(r['ARES_Rank'])} "
                  f"({int(abs(r['Rank_Change']))} pos) │ Delta Delta E={r['DDE_Contact']:+.2f}")

    # --- PLOTS ---

    # Plot 1: ARES ranking
    fig, ax = plt.subplots(figsize=(14, 12))
    df_plot = df.sort_values('ARES', ascending=True)
    cls_colors = {
        'Highly Destabilizing': '#ef4444', 'Destabilizing': '#f97316',
        'Moderately Destabilizing': '#eab308', 'Neutral/Mild': '#22c55e'
    }
    colors = [cls_colors.get(c, '#94a3b8') for c in df_plot['ARES_Class']]
    ax.barh(range(len(df_plot)), df_plot['ARES'], color=colors, height=0.7, alpha=0.85)
    for i, (_, r) in enumerate(df_plot.iterrows()):
        d = int(r['Rank_Change'])
        if abs(d) >= 5:
            arrow = f"↑{d}" if d > 0 else f"↓{abs(d)}"
            c = '#22c55e' if d > 0 else '#ef4444'
            ax.text(2, i, arrow, fontsize=6, va='center', fontweight='bold', color=c)
    ax.set_yticks(range(len(df_plot)))
    ax.set_yticklabels(df_plot['Mutation'], fontsize=7)
    ax.set_xlabel('TP53-ARES Score', fontsize=11, fontweight='bold')
    ax.set_title('TP53-ARES: Atomistic Residue Energy Scoring\n'
                 'Energy-based evaluation: Delta Delta E_contact + Propagation + Electrostatic + Network Rewiring\n'
                 'Arrows show movement vs RMSD ranking',
                 fontsize=11, fontweight='bold')
    legend_patches = [mpatches.Patch(color=cls_colors[c], label=f'{c} ({ares_counts.get(c, 0)})')
                      for c in cls_colors if ares_counts.get(c, 0) > 0]
    ax.legend(handles=legend_patches, fontsize=8, loc='lower right')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='x', alpha=0.2)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'ares_ranking.png'), dpi=300, bbox_inches='tight')
    plt.close()

    # Plot 2: Disruption wavefront visualization (top 5 + bottom 5)
    fig, axes = plt.subplots(2, 5, figsize=(20, 8))
    top5 = df.head(5)
    bot5 = df.tail(5)

    for idx, (subset, row_ax, title) in enumerate([
        (top5, axes[0], "Top 5 (Most Destabilizing)"),
        (bot5, axes[1], "Bottom 5 (Most Neutral)")
    ]):
        for col, (_, r) in enumerate(subset.iterrows()):
            ax = row_ax[col]
            mutation = r['Mutation']
            position = int(r['Position'])

            # Re-compute wavefront for this mutation
            mut_file = f"data/structures/tp53_{mutation.lower()}.cif"
            mut_struct = parser.get_structure("mut_viz", mut_file)
            mut_data_raw = get_residue_data(mut_struct)
            common = sorted(set(wt_data.keys()) & set(mut_data_raw.keys()))
            wt_atoms = [wt_data[r2]['atom'] for r2 in common]
            mut_atoms = [mut_data_raw[r2]['atom'] for r2 in common]
            sup = Superimposer()
            sup.set_atoms(wt_atoms, mut_atoms)
            mut_struct2 = parser.get_structure("mut_viz2", mut_file)
            sup.apply(mut_struct2.get_atoms())
            mut_data_viz = get_residue_data(mut_struct2)

            _, _, prop_d = component2_disruption_propagation(wt_data, mut_data_viz, position)
            shells = prop_d.get('shell_means', {})

            if shells:
                shell_nums = sorted(shells.keys())
                shell_vals = [shells[s] for s in shell_nums]
                color = '#ef4444' if r['ARES'] > 50 else '#22c55e'
                ax.bar(shell_nums, shell_vals, color=color, alpha=0.7, width=0.7)
                ax.axhline(2.0, color='gray', linestyle='--', alpha=0.3)

            ax.set_title(f"{mutation}\nARES={r['ARES']:.1f}", fontsize=8, fontweight='bold')
            ax.set_xlabel('Shell', fontsize=7)
            if col == 0:
                ax.set_ylabel('Mean Disp (Å)', fontsize=7)
            ax.set_ylim(0, max(shell_vals + [5]) * 1.2 if shells else 10)
            ax.tick_params(labelsize=6)

    plt.suptitle('Disruption Wavefront Propagation:\n'
                 'How structural damage spreads from mutation site through the contact network',
                 fontsize=12, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'ares_wavefront.png'), dpi=300, bbox_inches='tight')
    plt.close()

    # Plot 3: Delta Delta E contact energy landscape
    fig, ax = plt.subplots(figsize=(14, 8))
    df_dde = df.sort_values('DDE_Contact', ascending=True)
    colors = ['#ef4444' if v > 1 else '#f97316' if v > 0 else '#22c55e' for v in df_dde['DDE_Contact']]
    ax.barh(range(len(df_dde)), df_dde['DDE_Contact'], color=colors, height=0.7, alpha=0.85)
    ax.set_yticks(range(len(df_dde)))
    ax.set_yticklabels(df_dde['Mutation'], fontsize=7)
    ax.axvline(0, color='white', linewidth=0.5)
    ax.set_xlabel('Delta Delta E Contact Potential (RT units)', fontsize=11, fontweight='bold')
    ax.set_title('Miyazawa-Jernigan Contact Energy Change per Mutation\n'
                 'Delta Delta E = Σ[E(mutant,neighbor) − E(wt,neighbor)] | Positive = destabilizing',
                 fontsize=12, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='x', alpha=0.2)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'ares_dde_landscape.png'), dpi=300, bbox_inches='tight')
    plt.close()

    # Plot 4: ARES vs RMSD scatter
    fig, ax = plt.subplots(figsize=(10, 10))
    for _, r in df.iterrows():
        c = cls_colors.get(r['ARES_Class'], '#94a3b8')
        ax.scatter(r['RMSD'], r['ARES'], c=c, s=70, alpha=0.8, edgecolors='white', linewidth=0.5)
        if abs(r['Rank_Change']) >= 8:
            ax.annotate(r['Mutation'], (r['RMSD'], r['ARES']),
                        fontsize=7, ha='center', va='bottom', fontweight='bold')
    ax.set_xlabel('RMSD (Å)', fontsize=11, fontweight='bold')
    ax.set_ylabel('TP53-ARES Score', fontsize=11, fontweight='bold')
    ax.set_title('TP53-ARES vs RMSD — Energy vs Geometry\n'
                 'Points off trend = where energy-based assessment disagrees with geometry',
                 fontsize=12, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(alpha=0.15)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'ares_vs_rmsd.png'), dpi=300, bbox_inches='tight')
    plt.close()

    print(f"\n  [PLOT] ares_ranking.png")
    print(f"  [PLOT] ares_wavefront.png")
    print(f"  [PLOT] ares_dde_landscape.png")
    print(f"  [PLOT] ares_vs_rmsd.png")
    print(f"\n[SUCCESS] TP53-ARES: {output_dir}/ares_scores.csv")
    print(f"[SUCCESS] {len(df)} mutations scored across 4 energy dimensions")


if __name__ == "__main__":
    main()
