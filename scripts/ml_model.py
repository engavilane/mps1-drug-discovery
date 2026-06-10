import pandas as pd
import numpy as np
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.linear_model import Ridge
from sklearn.svm import SVR
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import LeaveOneOut, cross_val_score
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error
from sklearn.inspection import permutation_importance
import joblib
import warnings
warnings.filterwarnings("ignore")


# Paths 
LIGANDS_DIR  = Path("data/ligands/raw")
ADME_CSV     = "analysis/adme/adme_full.csv"
OUTPUT_DIR   = Path("analysis/ml_model")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# Load ADME data
adme_df = pd.read_csv(ADME_CSV)
print(f"Loaded {len(adme_df)} compounds\n")


# Compute Morgan fingerprints 
print("Computing Morgan fingerprints...")
fingerprints = []
valid_indices = []

for idx, row in adme_df.iterrows():
    sdf_path = LIGANDS_DIR / f"{row['ligand']}.sdf"
    mol = Chem.MolFromMolFile(str(sdf_path))
    if mol is not None:
        # Morgan fingerprint radius=2, 1024 bits
        fp = AllChem.GetMorganFingerprintAsBitVect(
            mol, radius=2, nBits=1024
        )
        fingerprints.append(list(fp))
        valid_indices.append(idx)
    else:
        print(f"   Could not parse: {row['ligand']}")

adme_df = adme_df.loc[valid_indices].reset_index(drop=True)
fp_df   = pd.DataFrame(fingerprints,
                        columns=[f"fp_{i}" for i in range(1024)])

print(f"   {len(fingerprints)} fingerprints computed\n")


# Build feature matrix 
# Physicochemical descriptors
physchem_cols = [
    "MW", "LogP", "HBD", "HBA", "TPSA",
    "RotBonds", "HeavyAtoms", "Rings", "AromaticRings",
    "gly605_contacts", "glu603_contacts"
]

physchem_df = adme_df[physchem_cols].copy()

# Combine fingerprints + physicochemical
X_full    = pd.concat([physchem_df, fp_df], axis=1)
X_physchem = physchem_df.copy()

# Target variable
y = adme_df["affinity_kcal"].values

print(f"Feature matrix shape: {X_full.shape}")
print(f"Target range: {y.min():.3f} to {y.max():.3f} kcal/mol\n")


# Remove low variance fingerprint features 
from sklearn.feature_selection import VarianceThreshold
selector = VarianceThreshold(threshold=0.05)
X_fp_filtered = selector.fit_transform(fp_df)
print(f"Fingerprint features after variance filter: "
      f"{X_fp_filtered.shape[1]}/1024\n")

X_combined = np.hstack([
    physchem_df.values,
    X_fp_filtered
])


# Scale features 
scaler = StandardScaler()


# Define models
models = {
    "Ridge Regression": Ridge(alpha=1.0),
    "Random Forest":    RandomForestRegressor(
                            n_estimators=100,
                            max_depth=5,
                            min_samples_leaf=3,
                            random_state=42
                        ),
    "SVR":              SVR(kernel="rbf", C=1.0, epsilon=0.1),
}


# LOOCV evaluation 
print("=" * 60)
print("LEAVE-ONE-OUT CROSS VALIDATION")
print("=" * 60)

loo     = LeaveOneOut()
results = {}

# Test on physicochemical only AND combined features
feature_sets = {
    "PhysChem only":    physchem_df.values,
    "PhysChem + Fingerprints": X_combined,
}

best_r2    = -999
best_model = None
best_name  = None
best_X     = None

for feat_name, X in feature_sets.items():
    print(f"\nFeature set: {feat_name} ({X.shape[1]} features)")
    print("-" * 40)

    X_scaled = scaler.fit_transform(X)

    for model_name, model in models.items():
        # LOOCV
        y_pred_loo = np.zeros(len(y))

        for train_idx, test_idx in loo.split(X_scaled):
            X_train = X_scaled[train_idx]
            X_test  = X_scaled[test_idx]
            y_train = y[train_idx]

            model.fit(X_train, y_train)
            y_pred_loo[test_idx] = model.predict(X_test)

        # Metrics
        r2   = r2_score(y, y_pred_loo)
        rmse = np.sqrt(mean_squared_error(y, y_pred_loo))
        mae  = mean_absolute_error(y, y_pred_loo)

        print(f"  {model_name}:")
        print(f"    R²={r2:.3f}  RMSE={rmse:.3f}  MAE={mae:.3f}")

        results[f"{feat_name} | {model_name}"] = {
            "R2": r2, "RMSE": rmse, "MAE": mae,
            "y_true": y, "y_pred": y_pred_loo
        }

        if r2 > best_r2:
            best_r2    = r2
            best_model = model
            best_name  = f"{feat_name} | {model_name}"
            best_X     = X_scaled.copy()

# Train best model on full data
print(f"\n{'=' * 60}")
print(f"BEST MODEL: {best_name}")
print(f"  R²={best_r2:.3f}")
print(f"{'=' * 60}")

best_model.fit(best_X, y)

# Feature importance (Random Forest only) 
if "Random Forest" in best_name:
    print("\nTop 15 most important features:")
    feat_names = (physchem_cols +
                  [f"fp_{i}" for i in
                   range(X_fp_filtered.shape[1])])
    importances = best_model.feature_importances_

    # Adjust length if needed
    top_n = min(15, len(importances))
    top_idx = np.argsort(importances)[::-1][:top_n]

    for rank, idx in enumerate(top_idx, 1):
        name = (feat_names[idx]
                if idx < len(feat_names) else f"feature_{idx}")
        print(f"  {rank:2}. {name}: {importances[idx]:.4f}")


# Save results 
# Results summary
results_summary = []
for name, metrics in results.items():
    results_summary.append({
        "model":    name,
        "R2":       round(metrics["R2"], 4),
        "RMSE":     round(metrics["RMSE"], 4),
        "MAE":      round(metrics["MAE"], 4),
    })

summary_df = pd.DataFrame(results_summary)
summary_df = summary_df.sort_values("R2", ascending=False)
summary_df.to_csv(OUTPUT_DIR / "model_comparison.csv", index=False)

# Predictions vs actual for best model
pred_df = pd.DataFrame({
    "ligand":       adme_df["ligand"],
    "affinity_true": y,
    "affinity_pred": results[best_name]["y_pred"],
    "error":        np.abs(y - results[best_name]["y_pred"]),
    "binds_hinge":  adme_df["binds_hinge"],
    "passes_adme":  adme_df["passes_adme"],
})
pred_df = pred_df.sort_values("affinity_true")
pred_df.to_csv(OUTPUT_DIR / "predictions.csv", index=False)

# Save best model + scaler
joblib.dump(best_model, OUTPUT_DIR / "best_model.pkl")
joblib.dump(scaler,     OUTPUT_DIR / "scaler.pkl")
joblib.dump(selector,   OUTPUT_DIR / "variance_selector.pkl")

print(f"\nModel saved  → {OUTPUT_DIR}/best_model.pkl")
print(f"Results      → {OUTPUT_DIR}/model_comparison.csv")
print(f"Predictions  → {OUTPUT_DIR}/predictions.csv")
print("\nModel comparison:")
print(summary_df.to_string(index=False))


# Ridge Regression coefficients
if "Ridge" in best_name and "PhysChem only" in best_name :
    print("\n" + "=" * 60)
    print("RIDGE REGRESSION COEFFICIENTS")
    print("(what drives binding affinity)")
    print("=" * 60)

    features = physchem_cols
    coefs = best_model.coef_

    coef_df = pd.DataFrame({
        "feature": features,
        "coefficient": coefs
    })
    coef_df["abs_coef"] = coef_df["coefficient"].abs()
    coef_df = coef_df.sort_values("abs_coef", ascending=False)

    for _, row in coef_df.iterrows() :
        direction = "increases" if row["coefficient"] < 0 \
                    else "decreases"
        print(f"  {row['feature']:20s}: {row['coefficient']:+.4f}"
              f"  → higher value {direction} affinity")

    coef_df.to_csv(OUTPUT_DIR / "ridge_coefficients.csv", index=False)
    print(f"\nCoefficients saved → {OUTPUT_DIR}/ridge_coefficients.csv")







    

