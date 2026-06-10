import pandas as pd
import numpy as np
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.linear_model import Ridge
from sklearn.svm import SVR
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error
from sklearn.feature_selection import VarianceThreshold
import joblib
import warnings
warnings.filterwarnings("ignore")


# Paths 
CHEMBL_CSV  = "analysis/ic50/chembl_mps1_ic50.csv"
OUTPUT_DIR  = Path("analysis/ml_model/chembl")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# Load ChEMBL data 
print("Loading ChEMBL Mps1 IC50 data...")
df = pd.read_csv(CHEMBL_CSV)
df = df.dropna(subset=["canonical_smiles", "IC50_nM", "pIC50"])
df = df[df["IC50_nM"] > 0]
df = df.drop_duplicates(subset="canonical_smiles")
print(f"  {len(df)} unique compounds after cleaning\n")


# Compute features 
print("Computing molecular features...")
fingerprints  = []
physchem_data = []
valid_indices = []

for idx, row in df.iterrows():
    mol = Chem.MolFromSmiles(row["canonical_smiles"])
    if mol is None:
        continue

    # Morgan fingerprint
    fp = AllChem.GetMorganFingerprintAsBitVect(
        mol, radius=2, nBits=1024
    )
    fingerprints.append(list(fp))

    # Physicochemical descriptors
    physchem_data.append({
        "MW":          Descriptors.MolWt(mol),
        "LogP":        Descriptors.MolLogP(mol),
        "HBD":         rdMolDescriptors.CalcNumHBD(mol),
        "HBA":         rdMolDescriptors.CalcNumHBA(mol),
        "TPSA":        Descriptors.TPSA(mol),
        "RotBonds":    rdMolDescriptors.CalcNumRotatableBonds(mol),
        "HeavyAtoms":  mol.GetNumHeavyAtoms(),
        "Rings":       rdMolDescriptors.CalcNumRings(mol),
        "AromaticRings": rdMolDescriptors.CalcNumAromaticRings(mol),
    })
    valid_indices.append(idx)

df          = df.loc[valid_indices].reset_index(drop=True)
fp_df       = pd.DataFrame(fingerprints,
                            columns=[f"fp_{i}" for i in range(1024)])
physchem_df = pd.DataFrame(physchem_data)

print(f"   {len(df)} compounds processed\n")


# Target variable 
y = df["pIC50"].values
print(f"pIC50 range: {y.min():.2f} to {y.max():.2f}")
print(f"pIC50 mean:  {y.mean():.2f} ± {y.std():.2f}\n")


# Remove low variance fingerprint features 
selector   = VarianceThreshold(threshold=0.05)
fp_filtered = selector.fit_transform(fp_df)
print(f"Fingerprint features after variance filter: "
      f"{fp_filtered.shape[1]}/1024\n")


# Feature sets 
X_physchem  = physchem_df.values
X_combined  = np.hstack([physchem_df.values, fp_filtered])

feature_sets = {
    "PhysChem only":           X_physchem,
    "PhysChem + Fingerprints": X_combined,
}


# Models 
models = {
    "Ridge Regression": Ridge(alpha=1.0),
    "Random Forest":    RandomForestRegressor(
                            n_estimators=200,
                            max_depth=10,
                            min_samples_leaf=2,
                            random_state=42,
                            n_jobs=-1
                        ),
    "SVR":              SVR(kernel="rbf", C=10.0, epsilon=0.1),
}


# Train/test split (80/20) 
# With 2794 compounds we can afford a proper split
print("=" * 60)
print("TRAIN/TEST SPLIT EVALUATION (80/20)")
print("=" * 60)

results      = {}
best_r2      = -999
best_model   = None
best_name    = None
best_X_train = None
best_X_test  = None
best_scaler  = None

for feat_name, X in feature_sets.items():
    print(f"\nFeature set: {feat_name} ({X.shape[1]} features)")
    print("-" * 40)

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42
    )

    scaler  = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled  = scaler.transform(X_test)

    for model_name, model in models.items():
        model.fit(X_train_scaled, y_train)
        y_pred = model.predict(X_test_scaled)

        r2   = r2_score(y_test, y_pred)
        rmse = np.sqrt(mean_squared_error(y_test, y_pred))
        mae  = mean_absolute_error(y_test, y_pred)

        print(f"  {model_name}:")
        print(f"    R²={r2:.3f}  RMSE={rmse:.3f}  MAE={mae:.3f}")

        results[f"{feat_name} | {model_name}"] = {
            "R2": r2, "RMSE": rmse, "MAE": mae,
            "y_true": y_test, "y_pred": y_pred
        }

        if r2 > best_r2:
            best_r2      = r2
            best_model   = model
            best_name    = f"{feat_name} | {model_name}"
            best_X_train = X_train_scaled.copy()
            best_X_test  = X_test_scaled.copy()
            best_scaler  = scaler


# Best model summary 
print(f"\n{'=' * 60}")
print(f"BEST MODEL: {best_name}")
print(f"  R²={best_r2:.3f}")
print(f"{'=' * 60}")


# Feature importance (Random Forest)
if "Random Forest" in best_name:
    print("\nTop 15 most important features:")
    physchem_names = list(physchem_df.columns)
    fp_names       = [f"fp_{i}" for i in
                      range(fp_filtered.shape[1])]

    if "PhysChem only" in best_name:
        feat_names = physchem_names
    else:
        feat_names = physchem_names + fp_names

    importances = best_model.feature_importances_
    top_idx     = np.argsort(importances)[::-1][:15]

    for rank, idx in enumerate(top_idx, 1):
        name = (feat_names[idx]
                if idx < len(feat_names) else f"feature_{idx}")
        print(f"  {rank:2}. {name}: {importances[idx]:.4f}")


# Save everything 
# Model comparison
summary = []
for name, metrics in results.items():
    summary.append({
        "model": name,
        "R2":    round(metrics["R2"],   4),
        "RMSE":  round(metrics["RMSE"], 4),
        "MAE":   round(metrics["MAE"],  4),
    })

summary_df = pd.DataFrame(summary).sort_values("R2", ascending=False)
summary_df.to_csv(OUTPUT_DIR / "chembl_model_comparison.csv",
                  index=False)

# Predictions vs actual
pred_df = pd.DataFrame({
    "pIC50_true": results[best_name]["y_true"],
    "pIC50_pred": results[best_name]["y_pred"],
    "error":      np.abs(results[best_name]["y_true"] -
                         results[best_name]["y_pred"]),
})
pred_df.to_csv(OUTPUT_DIR / "chembl_predictions.csv", index=False)

# Save best model
joblib.dump(best_model,  OUTPUT_DIR / "chembl_best_model.pkl")
joblib.dump(best_scaler, OUTPUT_DIR / "chembl_scaler.pkl")
joblib.dump(selector,    OUTPUT_DIR / "chembl_variance_selector.pkl")

print(f"\nModel saved   → {OUTPUT_DIR}/chembl_best_model.pkl")
print(f"Comparison    → {OUTPUT_DIR}/chembl_model_comparison.csv")
print(f"Predictions   → {OUTPUT_DIR}/chembl_predictions.csv")
print("\nFull model comparison:")
print(summary_df.to_string(index=False))
