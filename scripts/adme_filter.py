# Usage :
#   Phase 1 : python scripts/adme_filter.py
#   Phase 2 : python scripts/adme_filter.py \
#               --ligands      data/phase2/raw \
#               --interactions analysis/phase2/interactions/interaction_analysis.csv
#               --scores       docking/phase2_results/docking_scores.csv \
#               --output       analysis/phase2/adme

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, FilterCatalog
from rdkit.Chem.FilterCatalog import FilterCatalogParams
import pandas as pd
from pathlib import Path
import os
import argparse


# Argument parser
parser = argparse.ArgumentParser(
    description="ADME filtering with RDKit"
)
parser.add_argument("--ligands",      
                    default="data/ligands/raw")
parser.add_argument("--interactions", 
                    default="analysis/interactions/interaction_analysis.csv")
parser.add_argument("--scores",       
                    default="docking/results/docking_scores.csv")
parser.add_argument("--output",       
                    default="analysis/adme")
args = parser.parse_args()


# Paths
LIGANDS_DIR      = Path(args.ligands)
INTERACTIONS_CSV = args.interactions
SCORES_CSV       = args.scores
OUTPUT_DIR       = Path(args.output)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
 

# Lipinski Ro5 + drug-likeness thresholds 
# Based on Pugh et al. 2022 and standard drug discovery practice
FILTERS = {
    "MW":           (0,    500),    # Molecular weight ≤ 500 Da
    "LogP":         (-0.4,   5),    # Lipophilicity
    "HBD":          (0,      5),    # H-bond donors ≤ 5
    "HBA":          (0,     10),    # H-bond acceptors ≤ 10
    "TPSA":         (0,    140),    # Topological polar surface area
    "RotBonds":     (0,     10),    # Rotatable bonds ≤ 10
    "HeavyAtoms":   (19,    37),    # Heavy atom count (Pugh et al.)
}

# PAINS filter setup 
params = FilterCatalogParams()
params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
pains_catalog = FilterCatalog.FilterCatalog(params)

# Load interaction results 
interactions_df = pd.read_csv(INTERACTIONS_CSV)
scores_df       = pd.read_csv(SCORES_CSV)

# Merge scores + interactions
merged_df = scores_df.merge(interactions_df[
    ["ligand", "gly605_contacts", "glu603_contacts",
     "total_contacts", "binds_hinge"]
], on="ligand", how="left")

print(f"Computing ADME descriptors for {len(merged_df)} ligands\n")
print("=" * 60)

results = []

for _, row in merged_df.iterrows():
    ligand_name = row["ligand"]
    affinity    = row["affinity_kcal"]

    # Find SDF file 
    sdf_path = LIGANDS_DIR / f"{ligand_name}.sdf"
    if not sdf_path.exists():
        print(f" SDF not found: {ligand_name}")
        continue

    # Load molecule 
    mol = Chem.MolFromMolFile(str(sdf_path))
    if mol is None:
        print(f" Could not parse: {ligand_name}")
        continue

    # Compute descriptors 
    mw        = Descriptors.MolWt(mol)
    logp      = Descriptors.MolLogP(mol)
    hbd       = rdMolDescriptors.CalcNumHBD(mol)
    hba       = rdMolDescriptors.CalcNumHBA(mol)
    tpsa      = Descriptors.TPSA(mol)
    rotbonds  = rdMolDescriptors.CalcNumRotatableBonds(mol)
    heavyatoms = mol.GetNumHeavyAtoms()
    rings     = rdMolDescriptors.CalcNumRings(mol)
    arom_rings = rdMolDescriptors.CalcNumAromaticRings(mol)

    # Lipinski Ro5 check 
    ro5_violations = sum([
        mw    > 500,
        logp  > 5,
        hbd   > 5,
        hba   > 10,
    ])
    passes_ro5 = ro5_violations <= 1  # allow max 1 violation

    # Individual filter checks
    passes_mw      = FILTERS["MW"][0]        <= mw        <= FILTERS["MW"][1]
    passes_logp    = FILTERS["LogP"][0]      <= logp      <= FILTERS["LogP"][1]
    passes_hbd     = FILTERS["HBD"][0]       <= hbd       <= FILTERS["HBD"][1]
    passes_hba     = FILTERS["HBA"][0]       <= hba       <= FILTERS["HBA"][1]
    passes_tpsa    = FILTERS["TPSA"][0]      <= tpsa      <= FILTERS["TPSA"][1]
    passes_rotb    = FILTERS["RotBonds"][0]  <= rotbonds  <= FILTERS["RotBonds"][1]
    passes_heavy   = FILTERS["HeavyAtoms"][0] <= heavyatoms <= FILTERS["HeavyAtoms"][1]

    # PAINS check 
    pains_matches = pains_catalog.GetMatches(mol)
    has_pains     = len(pains_matches) > 0
    pains_alerts  = len(pains_matches)

    # Overall drug-likeness 
    passes_adme = (passes_ro5 and
                   passes_tpsa and
                   passes_rotb and
                   not has_pains)

    # Print summary 
    status = " PASS" if passes_adme else " FAIL"
    print(f"{status} | {ligand_name} ({affinity} kcal/mol)")
    print(f"       MW={mw:.1f}  LogP={logp:.2f}  "
          f"HBD={hbd}  HBA={hba}  "
          f"TPSA={tpsa:.1f}  RotB={rotbonds}  "
          f"PAINS={pains_alerts}")

    results.append({
        "ligand":           ligand_name,
        "affinity_kcal":    affinity,
        "binds_hinge":      row.get("binds_hinge", False),
        "gly605_contacts":  row.get("gly605_contacts", 0),
        "glu603_contacts":  row.get("glu603_contacts", 0),
        # Descriptors
        "MW":               round(mw, 2),
        "LogP":             round(logp, 2),
        "HBD":              hbd,
        "HBA":              hba,
        "TPSA":             round(tpsa, 2),
        "RotBonds":         rotbonds,
        "HeavyAtoms":       heavyatoms,
        "Rings":            rings,
        "AromaticRings":    arom_rings,
        # Filter results
        "passes_ro5":       passes_ro5,
        "ro5_violations":   ro5_violations,
        "passes_tpsa":      passes_tpsa,
        "passes_rotbonds":  passes_rotb,
        "has_pains":        has_pains,
        "pains_alerts":     pains_alerts,
        "passes_adme":      passes_adme,
    })

# Save full results 
results_df = pd.DataFrame(results)
results_df = results_df.sort_values("affinity_kcal")
full_csv = OUTPUT_DIR / "adme_full.csv"
results_df.to_csv(full_csv, index=False)

# Save filtered candidates only 
# Best candidates: pass ADME + bind hinge + good affinity
candidates_df = results_df[
    (results_df["passes_adme"]   == True) &
    (results_df["binds_hinge"]   == True)
].copy()
candidates_csv = OUTPUT_DIR / "adme_candidates.csv"
candidates_df.to_csv(candidates_csv, index=False)

# Summary 
passed = results_df[results_df["passes_adme"] == True]
failed = results_df[results_df["passes_adme"] == False]

print("\n" + "=" * 60)
print("ADME FILTER COMPLETE")
print("=" * 60)
print(f" Passed ADME:              {len(passed)}/{len(results_df)}")
print(f" Failed ADME:              {len(failed)}/{len(results_df)}")
print(f" Passed ADME + hinge:      {len(candidates_df)}/{len(results_df)}")

print(f"\nFailed ligands:")
for _, r in failed.iterrows():
    reasons = []
    if not r["passes_ro5"]:
        reasons.append(f"Ro5 violations={r['ro5_violations']}")
    if not r["passes_tpsa"]:
        reasons.append(f"TPSA={r['TPSA']}")
    if not r["passes_rotbonds"]:
        reasons.append(f"RotBonds={r['RotBonds']}")
    if r["has_pains"]:
        reasons.append(f"PAINS={r['pains_alerts']}")
    print(f"  {r['ligand']}: {', '.join(reasons)}")

print(f"\nTop candidates (ADME pass + hinge binding):")
for _, r in candidates_df.head(10).iterrows():
    print(f"  {r['ligand']}: {r['affinity_kcal']} kcal/mol | "
          f"MW={r['MW']} | LogP={r['LogP']} | "
          f"Gly605={r['gly605_contacts']} | "
          f"Glu603={r['glu603_contacts']}")

print(f"\nFull results → {full_csv}")
print(f"Candidates   → {candidates_csv}")

