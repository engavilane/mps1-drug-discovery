from chembl_webresource_client.new_client import new_client
import pandas as pd
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
import numpy as np

# Paths
ADME_CSV   = "analysis/adme/adme_full.csv"
LIGANDS_DIR = Path("data/ligands/raw")
OUTPUT_DIR  = Path("analysis/ic50")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# ChEMBL target for Mps1/TTK 
MPS1_CHEMBL_ID = "CHEMBL4523"


# Load our 45 compounds 
adme_df = pd.read_csv(ADME_CSV)
print(f"Searching IC50 for {len(adme_df)} compounds\n")


# Step 1: Get all IC50 data for Mps1 from ChEMBL 
print("Fetching all Mps1 bioactivity data from ChEMBL...")
activity  = new_client.activity
mps1_data = activity.filter(
    target_chembl_id=MPS1_CHEMBL_ID,
    standard_type="IC50",
    relation="=",
).only([
    "molecule_chembl_id",
    "standard_value",
    "standard_units",
    "canonical_smiles",
    "compound_name",
    "assay_chembl_id",
    "document_chembl_id",
])

mps1_df = pd.DataFrame(list(mps1_data))
print(f"  Found {len(mps1_df)} IC50 measurements for Mps1\n")

if mps1_df.empty:
    print("No data found — check ChEMBL target ID")
    exit()

# Keep only nM measurements and valid values
mps1_df = mps1_df[mps1_df["standard_units"] == "nM"]
mps1_df = mps1_df.dropna(subset=["standard_value", "canonical_smiles"])
mps1_df["IC50_nM"] = pd.to_numeric(
    mps1_df["standard_value"], errors="coerce"
)
mps1_df = mps1_df.dropna(subset=["IC50_nM"])
mps1_df = mps1_df[mps1_df["IC50_nM"] > 0]
print(f"  Valid IC50 values (nM): {len(mps1_df)}")

# Convert to pIC50 (negative log — higher = better binder)
mps1_df["pIC50"] = -np.log10(mps1_df["IC50_nM"] * 1e-9)
print(f"  pIC50 range: {mps1_df['pIC50'].min():.2f} "
      f"to {mps1_df['pIC50'].max():.2f}\n")


# Step 2: Compute fingerprints for ChEMBL compounds 
print("Computing fingerprints for ChEMBL compounds...")
chembl_fps = []
valid_chembl = []

for idx, row in mps1_df.iterrows():
    mol = Chem.MolFromSmiles(row["canonical_smiles"])
    if mol:
        fp = AllChem.GetMorganFingerprintAsBitVect(
            mol, radius=2, nBits=1024
        )
        chembl_fps.append(fp)
        valid_chembl.append(idx)

mps1_df = mps1_df.loc[valid_chembl].reset_index(drop=True)
print(f"  {len(chembl_fps)} valid fingerprints computed\n")

# Step 3: Match our ligands to ChEMBL by similarity 
print("Matching our 45 ligands to ChEMBL compounds...")
print("(Tanimoto similarity ≥ 0.4)\n")

matches = []

for _, lig_row in adme_df.iterrows():
    ligand_name = lig_row["ligand"]
    sdf_path    = LIGANDS_DIR / f"{ligand_name}.sdf"

    if not sdf_path.exists():
        continue

    mol = Chem.MolFromMolFile(str(sdf_path))
    if mol is None:
        continue

    lig_fp = AllChem.GetMorganFingerprintAsBitVect(
        mol, radius=2, nBits=1024
    )

    # Calculate Tanimoto similarity to all ChEMBL compounds
    similarities = DataStructs.BulkTanimotoSimilarity(
        lig_fp, chembl_fps
    )

    best_idx  = int(np.argmax(similarities))
    best_sim  = similarities[best_idx]

    if best_sim >= 0.4:
        best_match = mps1_df.iloc[best_idx]
        print(f"   {ligand_name}")
        print(f"    ChEMBL ID: {best_match['molecule_chembl_id']}")
        print(f"    IC50: {best_match['IC50_nM']:.2f} nM "
              f"(pIC50={best_match['pIC50']:.2f})")
        print(f"    Tanimoto: {best_sim:.3f}")

        matches.append({
            "ligand":           ligand_name,
            "affinity_kcal":    lig_row["affinity_kcal"],
            "chembl_id":        best_match["molecule_chembl_id"],
            "tanimoto":         round(best_sim, 3),
            "IC50_nM":          best_match["IC50_nM"],
            "pIC50":            round(best_match["pIC50"], 3),
            "smiles":           best_match["canonical_smiles"],
        })
    else:
        print(f"   {ligand_name} — "
              f"best match: {best_sim:.3f} (below threshold)")


# Save results 
matches_df = pd.DataFrame(matches)

if not matches_df.empty:
    matches_df = matches_df.sort_values("IC50_nM")
    matches_df.to_csv(OUTPUT_DIR / "ic50_matches.csv", index=False)

    print(f"\n{'=' * 60}")
    print(f"RESULTS")
    print(f"{'=' * 60}")
    print(f"Matched: {len(matches_df)}/{len(adme_df)} ligands")
    print(f"\nMatched compounds with IC50:")
    for _, r in matches_df.iterrows():
        print(f"  {r['ligand']:20s} | "
              f"Vina: {r['affinity_kcal']:6.3f} kcal/mol | "
              f"IC50: {r['IC50_nM']:8.2f} nM | "
              f"pIC50: {r['pIC50']:.2f}")


    # Correlation: Vina score vs pIC50 
    if len(matches_df) >= 5:
        from scipy import stats
        r, p = stats.pearsonr(
            matches_df["affinity_kcal"],
            matches_df["pIC50"]
        )
        print(f"\nCorrelation Vina score vs pIC50:")
        print(f"  Pearson r = {r:.3f}  (p={p:.4f})")
        if abs(r) > 0.5:
            print(f"  → Moderate-strong correlation ✓")
        else:
            print(f"  → Weak correlation — "
                  f"Vina scores don't fully reflect experimental affinity")

        matches_df.to_csv(OUTPUT_DIR / "ic50_matches.csv", index=False)
        print(f"\nResults saved → {OUTPUT_DIR}/ic50_matches.csv")
else:
    print("\nNo matches found at threshold 0.4")
    print("Try lowering threshold to 0.7 in the script")


# Save full ChEMBL dataset for Phase 2 ML 
mps1_df.to_csv(OUTPUT_DIR / "chembl_mps1_ic50.csv", index=False)
print(f"Full ChEMBL dataset → {OUTPUT_DIR}/chembl_mps1_ic50.csv")
print(f"({len(mps1_df)} compounds with experimental IC50)")
