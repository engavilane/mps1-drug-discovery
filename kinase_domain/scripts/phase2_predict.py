import pandas as pd
import numpy as np
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
import joblib
import warnings
warnings.filterwarnings("ignore")


# Paths 
ADME_CSV    = "kinase_domain/analysis/adme_plip/phase2/adme_candidates.csv"
LIGANDS_DIR = Path("kinase_domain/data/phase2/raw")
MODEL_PATH  = "kinase_domain/analysis/ml_model/chembl/chembl_best_model.pkl"
SCALER_PATH = "kinase_domain/analysis/ml_model/chembl/chembl_scaler.pkl"
SELECTOR_PATH = "kinase_domain/analysis/ml_model/chembl/chembl_variance_selector.pkl"
OUTPUT_DIR  = Path("kinase_domain/analysis/phase2")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# Load model 
print("Loading ChEMBL SVR model...")
model    = joblib.load(MODEL_PATH)
scaler   = joblib.load(SCALER_PATH)
selector = joblib.load(SELECTOR_PATH)
print("   Model loaded\n")


# Load ADME candidates 
adme_df = pd.read_csv(ADME_CSV)
print(f"Predicting pIC50 for {len(adme_df)} candidates\n")


# Predict 
results = []

for _, row in adme_df.iterrows():
    name = row["ligand"]

    # Find SDF — name in adme_df may have _out suffix
    clean_name = name.replace("_out", "")
    sdf_path   = LIGANDS_DIR / f"{clean_name}.sdf"

    # Try alternate naming if not found
    if not sdf_path.exists():
        matches = list(LIGANDS_DIR.glob(f"*{clean_name}*.sdf"))
        if not matches:
            print(f"   SDF not found: {name}")
            continue
        sdf_path = matches[0]

    mol = Chem.MolFromMolFile(str(sdf_path))
    if mol is None:
        print(f"   Could not parse: {name}")
        continue


    # Compute features
    fp = list(AllChem.GetMorganFingerprintAsBitVect(
        mol, radius=2, nBits=1024
    ))
    physchem = [
        Descriptors.MolWt(mol),
        Descriptors.MolLogP(mol),
        rdMolDescriptors.CalcNumHBD(mol),
        rdMolDescriptors.CalcNumHBA(mol),
        Descriptors.TPSA(mol),
        rdMolDescriptors.CalcNumRotatableBonds(mol),
        mol.GetNumHeavyAtoms(),
        rdMolDescriptors.CalcNumRings(mol),
        rdMolDescriptors.CalcNumAromaticRings(mol),
    ]

    fp_arr      = np.array(fp).reshape(1, -1)
    fp_filtered = selector.transform(fp_arr)
    X           = np.hstack([
        np.array(physchem).reshape(1, -1),
        fp_filtered
    ])
    X_scaled    = scaler.transform(X)
    pIC50_pred  = float(model.predict(X_scaled)[0])
    ic50_nM     = 10 ** (9 - pIC50_pred)

    results.append({
        "ligand":          name,
        "seed":            name.split("_")[1] + "_" +
                           name.split("_")[2],
        "affinity_kcal":   row["affinity_kcal"],
        "gly605_contacts": row["gly605_contacts"],
        "glu603_contacts": row["glu603_contacts"],
        "MW":              row["MW"],
        "LogP":            row["LogP"],
        "HBD":             row["HBD"],
        "HBA":             row["HBA"],
        "TPSA":            row["TPSA"],
        "pIC50_pred":      round(pIC50_pred, 3),
        "IC50_pred_nM":    round(ic50_nM, 2),
    })


# Sort and save 
results_df = pd.DataFrame(results)
results_df = results_df.sort_values("pIC50_pred", ascending=False)
results_df.to_csv(
    OUTPUT_DIR / "phase2_final_candidates.csv", index=False
)


# Print final ranking 
print(f"{'=' * 70}")
print(f"PHASE 2 FINAL RESULTS — TOP 20 NOVEL MPS1 INHIBITOR CANDIDATES")
print(f"{'=' * 70}")
print(f"{'Rank':<5} {'Ligand':<35} {'Vina':>7} "
      f"{'pIC50':>7} {'IC50(nM)':>10} "
      f"{'Gly605':>7} {'Glu603':>7}")
print("-" * 75)

for rank, (_, row) in enumerate(results_df.head(20).iterrows(), 1):
    print(f"{rank:<5} {row['ligand']:<35} "
          f"{row['affinity_kcal']:>7.3f} "
          f"{row['pIC50_pred']:>7.3f} "
          f"{row['IC50_pred_nM']:>10.2f} "
          f"{row['gly605_contacts']:>7.0f} "
          f"{row['glu603_contacts']:>7.0f}")

print(f"\n{'=' * 70}")
print(f"Total candidates: {len(results_df)}")
print(f"pIC50 range: {results_df['pIC50_pred'].min():.3f} "
      f"to {results_df['pIC50_pred'].max():.3f}")
print(f"Best predicted IC50: "
      f"{results_df['IC50_pred_nM'].min():.2f} nM")
print(f"\nResults → {OUTPUT_DIR}/phase2_final_candidates.csv")



# Combined score 
# Equal weighting (0.5/0.5)
results_df["vina_norm"] = (
    results_df["affinity_kcal"].min() -
    results_df["affinity_kcal"]
) / (
    results_df["affinity_kcal"].min() -
    results_df["affinity_kcal"].max()
)
results_df["pIC50_norm"] = (
    results_df["pIC50_pred"] -
    results_df["pIC50_pred"].min()
) / (
    results_df["pIC50_pred"].max() -
    results_df["pIC50_pred"].min()
)

# Equal weighting — primary ranking
results_df["combined_score"] = (
    results_df["vina_norm"] +
    results_df["pIC50_norm"]
) / 2


# Sensitivity analysis 
# Test 3 weighting schemes to confirm ranking stability
print("\n" + "=" * 60)
print("COMBINED SCORE SENSITIVITY ANALYSIS")
print("=" * 60)

weights = {
    "Vina-heavy  (0.7/0.3)": (0.7, 0.3),
    "Equal       (0.5/0.5)": (0.5, 0.5),
    "pIC50-heavy (0.3/0.7)": (0.3, 0.7),
}

# Get top 10 under equal weighting (our primary ranking)
results_df = results_df.sort_values(
    "combined_score", ascending=False
)
top10_equal = results_df.head(10)["ligand"].tolist()

print(f"\n{'Ligand':<35} {'0.7/0.3':>8} "
      f"{'0.5/0.5':>8} {'0.3/0.7':>8}")
print("-" * 60)

rankings = {}
for w_label, (w_vina, w_pic50) in weights.items():
    col = f"score_{w_label.split()[0]}"
    results_df[col] = (
        w_vina  * results_df["vina_norm"] +
        w_pic50 * results_df["pIC50_norm"]
    )
    ranked = results_df.sort_values(col, ascending=False)
    rankings[w_label] = ranked["ligand"].tolist()

for lig in top10_equal:
    ranks = []
    for w_label in weights:
        col  = f"score_{w_label.split()[0]}"
        rank = results_df.sort_values(
            col, ascending=False
        ).reset_index(drop=True)
        rank = rank[rank["ligand"] == lig].index[0] + 1
        ranks.append(rank)
    name_short = lig.replace("phase2_", "")[:33]
    print(f"  {name_short:<33} "
          f"{ranks[0]:>8} "
          f"{ranks[1]:>8} "
          f"{ranks[2]:>8}")

# Spearman correlation between weighting schemes
from scipy.stats import spearmanr

r1 = results_df.sort_values(
    "score_Vina-heavy", ascending=False
)["ligand"].tolist()
r2 = results_df.sort_values(
    "combined_score", ascending=False
)["ligand"].tolist()
r3 = results_df.sort_values(
    "score_pIC50-heavy", ascending=False
)["ligand"].tolist()

n   = len(r1)
rho_12, _ = spearmanr(
    [r1.index(l) for l in r2],
    [r2.index(l) for l in r2]
)
rho_23, _ = spearmanr(
    [r2.index(l) for l in r3],
    [r3.index(l) for l in r3]
)
rho_13, _ = spearmanr(
    [r1.index(l) for l in r3],
    [r3.index(l) for l in r3]
)

print(f"\nSpearman ρ between weighting schemes:")
print(f"  Vina-heavy vs Equal:       {rho_12:.3f}")
print(f"  Equal vs pIC50-heavy:      {rho_23:.3f}")
print(f"  Vina-heavy vs pIC50-heavy: {rho_13:.3f}")

if min(rho_12, rho_23, rho_13) > 0.8:
    print("✓ Rankings highly consistent across weightings")
    print("  Equal weighting validated by sensitivity analysis")
elif min(rho_12, rho_23, rho_13) > 0.6:
    print("⚠ Moderate consistency — top candidates robust")
    print("  but mid-tier ranking sensitive to weighting")
else:
    print("✗ Rankings sensitive to weighting choice")
    print("  Consider reporting results for all three schemes")

# Save sensitivity results
sensitivity_df = results_df[[
    "ligand", "affinity_kcal", "pIC50_pred",
    "vina_norm", "pIC50_norm", "combined_score"
]].copy()
sensitivity_df.to_csv(
    OUTPUT_DIR / "sensitivity_analysis.csv", index=False
)
print(f"\nSensitivity results → "
      f"{OUTPUT_DIR}/sensitivity_analysis.csv")

# Drop temporary weighting columns
drop_cols = [c for c in results_df.columns
             if c.startswith("score_")]
results_df = results_df.drop(columns=drop_cols)


# Synthetic Accessibility Scores 
print(f"\nComputing Synthetic Accessibility scores...")
import sys
from rdkit.Chem import RDConfig
sys.path.append(RDConfig.RDContribDir)
from SA_Score import sascorer

sa_scores = []
for _, row in results_df.iterrows():
    name  = row["ligand"]
    clean = name.replace("_out", "")
    sdf   = LIGANDS_DIR / f"{clean}.sdf"
    if not sdf.exists():
        matches = list(LIGANDS_DIR.glob(f"*{clean}*.sdf"))
        sdf = matches[0] if matches else None

    if sdf is None:
        sa_scores.append(None)
        continue

    mol = Chem.MolFromMolFile(str(sdf))
    if mol is None:
        sa_scores.append(None)
        continue

    sa_scores.append(round(sascorer.calculateScore(mol), 3))

results_df["sa_score"] = sa_scores
results_df["sa_class"] = results_df["sa_score"].apply(
    lambda s: "easy" if s is not None and s <= 3
    else ("moderate" if s is not None and s <= 6
    else ("hard" if s is not None else "unknown"))
)

print(f"  Easy (1-3):     {(results_df['sa_class']=='easy').sum()}")
print(f"  Moderate (4-6): {(results_df['sa_class']=='moderate').sum()}")
print(f"  Hard (7-10):    {(results_df['sa_class']=='hard').sum()}")
print(f"  Mean SA score:  {results_df['sa_score'].mean():.2f}")

# Sort by combined score instead
results_df = results_df.sort_values(
    "combined_score", ascending=False
)
