"""
Y-randomisation test for Model 2 (ChEMBL SVR pIC50 prediction).

Validates that the model has learned genuine structure-activity
relationships rather than chance correlations.

Method:
  1. Take the same feature matrix X and target vector y
  2. Shuffle y randomly 100 times
  3. Rebuild the model each time with shuffled labels
  4. Compare real model R² to distribution of randomised R²

Expected result:
  Real model R² >> mean randomised R²
  If randomised models achieve similar R² → model is invalid
"""

import pandas as pd
import numpy as np
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from sklearn.svm import SVR
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score
from sklearn.feature_selection import VarianceThreshold
from sklearn.pipeline import Pipeline
import joblib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')


# Paths 
CHEMBL_CSV = "kinase_domain/analysis/ic50/chembl_mps1_ic50.csv"
MODEL_PATH = "kinase_domain/analysis/ml_model/chembl/chembl_best_model.pkl"
SCALER_PATH = "kinase_domain/analysis/ml_model/chembl/chembl_scaler.pkl"
SELECTOR_PATH = "kinase_domain/analysis/ml_model/chembl/chembl_variance_selector.pkl"
OUTPUT_DIR = Path("kinase_domain/analysis/ml_model/chembl")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

N_ITERATIONS = 100
RANDOM_SEED  = 42


# Load and prepare data 
print("Loading ChEMBL data...")
df = pd.read_csv(CHEMBL_CSV)
df = df.dropna(subset=["canonical_smiles", "pIC50"])
df = df[df["IC50_nM"] > 0]
df = df.drop_duplicates(subset="canonical_smiles")
print(f"  {len(df)} compounds\n")

print("Computing features...")
fingerprints  = []
physchem_data = []
valid_indices = []

for idx, row in df.iterrows():
    mol = Chem.MolFromSmiles(row["canonical_smiles"])
    if mol is None:
        continue
    fp = list(AllChem.GetMorganFingerprintAsBitVect(
        mol, radius=2, nBits=1024
    ))
    fingerprints.append(fp)
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

# Apply variance filter
selector    = VarianceThreshold(threshold=0.05)
fp_filtered = selector.fit_transform(fp_df)

# Combined features (best configuration)
X = np.hstack([physchem_df.values, fp_filtered])
y = df["pIC50"].values

print(f"  Features: {X.shape[1]}")
print(f"  Compounds: {len(y)}\n")


# Real model performance 
print("Computing real model performance...")
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=RANDOM_SEED
)
scaler         = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled  = scaler.transform(X_test)

real_model = SVR(kernel="rbf", C=10.0, epsilon=0.1)
real_model.fit(X_train_scaled, y_train)
y_pred     = real_model.predict(X_test_scaled)
real_r2    = r2_score(y_test, y_pred)
print(f"  Real model R²: {real_r2:.4f}\n")


# Y-randomisation
print(f"Running Y-randomisation ({N_ITERATIONS} iterations)...")
print("=" * 55)

rand_r2s = []
rng      = np.random.RandomState(RANDOM_SEED)

for i in range(N_ITERATIONS):
    # Shuffle labels
    y_shuffled = rng.permutation(y)

    # Same split indices
    X_train_r, X_test_r, y_train_r, y_test_r = train_test_split(
        X, y_shuffled,
        test_size=0.2,
        random_state=RANDOM_SEED
    )

    scaler_r        = StandardScaler()
    X_train_r_scaled = scaler_r.fit_transform(X_train_r)
    X_test_r_scaled  = scaler_r.transform(X_test_r)

    rand_model = SVR(kernel="rbf", C=10.0, epsilon=0.1)
    rand_model.fit(X_train_r_scaled, y_train_r)
    y_pred_r   = rand_model.predict(X_test_r_scaled)
    r2_r       = r2_score(y_test_r, y_pred_r)
    rand_r2s.append(r2_r)

    if (i + 1) % 10 == 0:
        print(f"  Iteration {i+1:3d}/100 — "
              f"running mean R²: {np.mean(rand_r2s):.4f}")

rand_r2s = np.array(rand_r2s)


# Statistical analysis 
mean_rand   = rand_r2s.mean()
std_rand    = rand_r2s.std()
max_rand    = rand_r2s.max()
z_score     = (real_r2 - mean_rand) / std_rand
p_value     = (rand_r2s >= real_r2).sum() / N_ITERATIONS

print(f"\n{'=' * 55}")
print(f"Y-RANDOMISATION RESULTS")
print(f"{'=' * 55}")
print(f"Real model R²:          {real_r2:.4f}")
print(f"Randomised R² mean:     {mean_rand:.4f} ± {std_rand:.4f}")
print(f"Randomised R² max:      {max_rand:.4f}")
print(f"Randomised R² range:    "
      f"[{rand_r2s.min():.4f}, {max_rand:.4f}]")
print(f"Z-score:                {z_score:.2f}")
print(f"Empirical p-value:      {p_value:.4f}")
print(f"  (fraction of random models ≥ real R²)")

if p_value < 0.05:
    print(f"\n✓ Y-RANDOMISATION PASSED")
    print(f"  Real model significantly outperforms chance")
    print(f"  (p < 0.05) — genuine SAR learned")
elif p_value < 0.10:
    print(f"\n⚠ BORDERLINE (p={p_value:.3f})")
    print(f"  Model marginally better than chance")
else:
    print(f"\n✗ Y-RANDOMISATION FAILED")
    print(f"  Model does not outperform random — invalid QSAR")


# Save results 
results_df = pd.DataFrame({
    "iteration":   range(1, N_ITERATIONS + 1),
    "r2_random":   rand_r2s
})
results_df.to_csv(
    OUTPUT_DIR / "y_randomisation_results.csv", index=False
)

summary_df = pd.DataFrame({
    "metric": [
        "real_r2", "mean_random_r2", "std_random_r2",
        "max_random_r2", "z_score", "p_value"
    ],
    "value": [
        real_r2, mean_rand, std_rand,
        max_rand, z_score, p_value
    ]
})
summary_df.to_csv(
    OUTPUT_DIR / "y_randomisation_summary.csv", index=False
)


# Plot
fig, ax = plt.subplots(figsize=(8, 5))

ax.hist(rand_r2s, bins=20, color='steelblue',
        alpha=0.7, edgecolor='white',
        label=f'Randomised models (n={N_ITERATIONS})')
ax.axvline(real_r2, color='red', linewidth=2.5,
           linestyle='--',
           label=f'Real model R²={real_r2:.3f}')
ax.axvline(mean_rand, color='navy', linewidth=1.5,
           linestyle=':',
           label=f'Mean random R²={mean_rand:.3f}')

ax.set_xlabel('R² (test set)', fontsize=13)
ax.set_ylabel('Count', fontsize=13)
ax.set_title('Y-Randomisation Test — Model 2 (SVR, ChEMBL pIC50)',
             fontsize=13)
ax.legend(fontsize=11)
ax.text(0.98, 0.95,
        f'Z={z_score:.1f}\np={p_value:.3f}',
        transform=ax.transAxes,
        ha='right', va='top', fontsize=11,
        bbox=dict(boxstyle='round', facecolor='white',
                  alpha=0.8))

plt.tight_layout()
plt.savefig(
    OUTPUT_DIR / "y_randomisation_plot.png",
    dpi=300, bbox_inches='tight'
)
print(f"\nPlot → {OUTPUT_DIR}/y_randomisation_plot.png")
print(f"Results → {OUTPUT_DIR}/y_randomisation_results.csv")
print(f"Summary → {OUTPUT_DIR}/y_randomisation_summary.csv")

