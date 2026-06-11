import pandas as pd
import numpy as np
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
import joblib
import warnings
warnings.filterwarnings("ignore")


# Paths 
ADME_CSV    = "analysis/phase2/adme/adme_candidates.csv"
LIGANDS_DIR = Path("data/phase2/raw")
MODEL_PATH  = "analysis/ml_model/chembl/chembl_best_model.pkl"
SCALER_PATH = "analysis/ml_model/chembl/chembl_scaler.pkl"
SELECTOR_PATH = "analysis/ml_model/chembl/chembl_variance_selector.pkl"
OUTPUT_DIR  = Path("analysis/phase2")
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

# Add combined score
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
results_df["combined_score"] = (
    results_df["vina_norm"] +
    results_df["pIC50_norm"]
) / 2

# Sort by combined score instead
results_df = results_df.sort_values(
    "combined_score", ascending=False
)
