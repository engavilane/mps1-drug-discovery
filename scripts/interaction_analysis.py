import numpy as np
import pandas as pd
from pathlib import Path

# ── Paths ─────────────────────────────────────────────────
RECEPTOR_PATH = "data/receptor/receptor.pdbqt"
RESULTS_DIR   = Path("docking/results")
SCORES_CSV    = "docking/results/docking_scores.csv"
OUTPUT_DIR    = Path("analysis/interactions")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# H-bond distance cutoff (donor-acceptor)
DISTANCE_CUTOFF = 3.5  # Angstroms

# PDBQT parser 
def parse_pdbqt(filepath):
    """Parse a PDBQT file and return atoms as list of dicts."""
    atoms = []
    with open(filepath, "r") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    atom = {
                        "name":    line[12:16].strip(),
                        "resname": line[17:20].strip(),
                        "resid":   int(line[22:26].strip()),
                        "x":       float(line[30:38]),
                        "y":       float(line[38:46]),
                        "z":       float(line[46:54]),
                    }
                    atoms.append(atom)
                except ValueError:
                    continue
    return atoms

# Load receptor once 
print("Loading receptor...")
receptor_atoms = parse_pdbqt(RECEPTOR_PATH)

# Extract Gly605 and Glu603 atoms
gly605_atoms = [a for a in receptor_atoms if a["resid"] == 605]
glu603_atoms = [a for a in receptor_atoms if a["resid"] == 603]

print(f"  Gly605 atoms found: {len(gly605_atoms)}")
print(f"  Glu603 atoms found: {len(glu603_atoms)}")

if not gly605_atoms and not glu603_atoms:
    print("\n   WARNING: Neither Gly605 nor Glu603 found!")
    print("  Checking residue IDs in receptor...")
    resids = sorted(set(a["resid"] for a in receptor_atoms))
    print(f"  Available resids: {resids[:20]} ...")
    # Find residues around expected range
    nearby = [a for a in receptor_atoms if 595 <= a["resid"] <= 615]
    nearby_resids = sorted(set(a["resid"] for a in nearby))
    print(f"  Residues 595-615: {nearby_resids}")

# Load docking scores 
scores_df = pd.read_csv(SCORES_CSV)
results = []

print(f"\nAnalysing {len(scores_df)} docked ligands\n")
print("=" * 60)

for _, row in scores_df.iterrows():
    ligand_name  = row["ligand"]
    affinity     = row["affinity_kcal"]
    output_pdbqt = row["output_file"]

    if not Path(output_pdbqt).exists():
        print(f" File not found: {output_pdbqt}")
        continue

    print(f"Analysing: {ligand_name} ({affinity} kcal/mol)")

    # Parse ligand (first pose only) 
    ligand_atoms = []
    with open(output_pdbqt, "r") as f:
        for line in f:
            if line.startswith("ENDMDL"):
                break  # only first pose
            if line.startswith(("ATOM", "HETATM")):
                try:
                    ligand_atoms.append(np.array([
                        float(line[30:38]),
                        float(line[38:46]),
                        float(line[46:54])
                    ]))
                except ValueError:
                    continue

    if not ligand_atoms:
        print(f"   No ligand atoms parsed")
        continue

    ligand_coords = np.array(ligand_atoms)

    # Calculate contacts 
    def get_contacts(res_atoms, lig_coords, cutoff):
        contacts = []
        for res_atom in res_atoms:
            res_coord = np.array([
                res_atom["x"], res_atom["y"], res_atom["z"]
            ])
            dists = np.linalg.norm(lig_coords - res_coord, axis=1)
            min_dist = np.min(dists)
            if min_dist <= cutoff:
                contacts.append({
                    "res_atom":    res_atom["name"],
                    "resid":       res_atom["resid"],
                    "distance_A":  round(float(min_dist), 2)
                })
        return contacts

    gly605_contacts = get_contacts(gly605_atoms, ligand_coords,
                                   DISTANCE_CUTOFF)
    glu603_contacts = get_contacts(glu603_atoms, ligand_coords,
                                   DISTANCE_CUTOFF)

    total = len(gly605_contacts) + len(glu603_contacts)
    binds_hinge = total > 0

    status = " HINGE CONTACT" if binds_hinge else " NO HINGE CONTACT"
    print(f"  {status} | Gly605: {len(gly605_contacts)} | "
          f"Glu603: {len(glu603_contacts)}")

    if gly605_contacts:
        best = min(gly605_contacts, key=lambda x: x["distance_A"])
        print(f"    Gly605 closest: {best['res_atom']} "
              f"@ {best['distance_A']} Å")
    if glu603_contacts:
        best = min(glu603_contacts, key=lambda x: x["distance_A"])
        print(f"    Glu603 closest: {best['res_atom']} "
              f"@ {best['distance_A']} Å")

    results.append({
        "ligand":          ligand_name,
        "affinity_kcal":   affinity,
        "gly605_contacts": len(gly605_contacts),
        "glu603_contacts": len(glu603_contacts),
        "total_contacts":  total,
        "binds_hinge":     binds_hinge,
    })

# Save results 
results_df = pd.DataFrame(results)
results_df = results_df.sort_values("affinity_kcal")
output_csv = OUTPUT_DIR / "interaction_analysis.csv"
results_df.to_csv(output_csv, index=False)

# Summary 
correct   = results_df[results_df["binds_hinge"] == True]
incorrect = results_df[results_df["binds_hinge"] == False]

print("\n" + "=" * 60)
print("INTERACTION ANALYSIS COMPLETE")
print("=" * 60)
print(f" Hinge binders:    {len(correct)}/{len(results_df)}")
print(f" No hinge contact: {len(incorrect)}/{len(results_df)}")

if len(correct) > 0:
    print(f"\nTop 5 hinge binders by affinity:")
    for _, r in correct.head(5).iterrows():
        print(f"  {r['ligand']}: {r['affinity_kcal']} kcal/mol | "
              f"Gly605: {r['gly605_contacts']} | "
              f"Glu603: {r['glu603_contacts']}")

if len(incorrect) > 0:
    print(f"\nLigands with no hinge contact:")
    for _, r in incorrect.iterrows():
        print(f"  {r['ligand']}: {r['affinity_kcal']} kcal/mol")

print(f"\nResults saved to {output_csv}")





