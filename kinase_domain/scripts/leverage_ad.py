"""
Leverage-based Applicability Domain analysis.

Williams plot: leverage (h) vs standardised residuals
Compounds with h > h* are outside the AD.

h* = 3(k+1)/n
where k = number of features
      n = number of training compounds

Reference: Tropsha et al. 2003, QSAR Comb Sci 22:69-77
"""

import pandas as pd
import numpy as np
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import VarianceThreshold
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

CHEMBL_CSV     = "kinase_domain/analysis/ic50/expanded_mps1_activity.csv"
CANDIDATES_CSV = "kinase_domain/analysis/phase2/phase2_final_candidates.csv"
PHASE2B_CSV    = "kinase_domain/analysis/phase2b/admet_candidates_clean.csv"
OUTPUT_DIR     = Path("kinase_domain/analysis/ml_model/chembl")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

print("Leverage-based Applicability Domain Analysis")
print("=" * 55)

def get_features(smiles_list):
    fps       = []
    physchem  = []
    valid_idx = []

    for i, smi in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(str(smi))
        if mol is None:
            continue
        fp = list(AllChem.GetMorganFingerprintAsBitVect(
            mol, radius=2, nBits=1024
        ))
        fps.append(fp)
        physchem.append([
            Descriptors.MolWt(mol),
            Descriptors.MolLogP(mol),
            rdMolDescriptors.CalcNumHBD(mol),
            rdMolDescriptors.CalcNumHBA(mol),
            Descriptors.TPSA(mol),
            rdMolDescriptors.CalcNumRotatableBonds(mol),
            mol.GetNumHeavyAtoms(),
            rdMolDescriptors.CalcNumRings(mol),
            rdMolDescriptors.CalcNumAromaticRings(mol),
        ])
        valid_idx.append(i)

    X = np.hstack([
        np.array(physchem),
        np.array(fps)
    ])
    return X, valid_idx

# Load training set
print("Loading ChEMBL training set...")
train_df = pd.read_csv(CHEMBL_CSV)
train_df = train_df.dropna(
    subset=['smiles', 'pIC50']
).drop_duplicates('smiles')
print(f"  Training compounds: {len(train_df)}")

X_train, _ = get_features(train_df['smiles'].tolist())
print(f"  Features: {X_train.shape[1]}")

# Variance filter
selector    = VarianceThreshold(threshold=0.05)
X_train_filt = selector.fit_transform(X_train)
print(f"  After variance filter: {X_train_filt.shape[1]}")

# Scale
scaler       = StandardScaler()
X_train_sc   = scaler.fit_transform(X_train_filt)

n, k         = X_train_sc.shape
h_star       = 3 * (k + 1) / n
print(f"\n  n={n}, k={k}")
print(f"  AD threshold h* = {h_star:.4f}")

# Hat matrix diagonal (leverage)
# h_i = x_i^T (X^T X)^{-1} x_i
print("\nComputing leverage for training set...")
XtX     = X_train_sc.T @ X_train_sc
XtX_inv = np.linalg.pinv(XtX)

train_h = np.array([
    float(X_train_sc[i] @ XtX_inv @ X_train_sc[i])
    for i in range(n)
])

print(f"  Training leverage range: "
      f"{train_h.min():.4f} - {train_h.max():.4f}")
print(f"  Within AD: "
      f"{(train_h <= h_star).sum()}/{n} "
      f"({100*(train_h<=h_star).mean():.1f}%)")

def compute_leverage(smiles_list, name="candidates"):
    print(f"\nComputing leverage for {name}...")
    X, valid_idx = get_features(smiles_list)
    if len(X) == 0:
        return None, None, []

    X_filt = selector.transform(X)
    X_sc   = scaler.transform(X_filt)

    leverages = np.array([
        float(X_sc[i] @ XtX_inv @ X_sc[i])
        for i in range(len(X_sc))
    ])

    in_ad = leverages <= h_star
    print(f"  Within AD (h ≤ {h_star:.4f}): "
          f"{in_ad.sum()}/{len(leverages)} "
          f"({100*in_ad.mean():.1f}%)")
    print(f"  Leverage range: "
          f"{leverages.min():.4f} - {leverages.max():.4f}")

    return leverages, in_ad, valid_idx

# Phase 2A candidates
results_2a = {}
if Path(CANDIDATES_CSV).exists():
    cand_df    = pd.read_csv(CANDIDATES_CSV)
    lev, in_ad, vidx = compute_leverage(
        cand_df['smiles'].tolist()
        if 'smiles' in cand_df.columns
        else [], "Phase 2A"
    )
    if lev is not None:
        cand_df['leverage']   = np.nan
        cand_df['in_ad_leverage'] = False
        for i, orig_i in enumerate(vidx):
            cand_df.loc[orig_i, 'leverage']       = lev[i]
            cand_df.loc[orig_i, 'in_ad_leverage'] = in_ad[i]
        cand_df.to_csv(
            OUTPUT_DIR / "phase2a_leverage_ad.csv",
            index=False
        )
        results_2a = {
            'leverage': lev,
            'in_ad':    in_ad
        }

# Phase 2B candidates
results_2b = {}
if Path(PHASE2B_CSV).exists():
    cand_df    = pd.read_csv(PHASE2B_CSV)
    lev, in_ad, vidx = compute_leverage(
        cand_df['smiles'].tolist()
        if 'smiles' in cand_df.columns
        else [], "Phase 2B"
    )
    if lev is not None:
        cand_df['leverage']       = np.nan
        cand_df['in_ad_leverage'] = False
        for i, orig_i in enumerate(vidx):
            cand_df.loc[orig_i, 'leverage']       = lev[i]
            cand_df.loc[orig_i, 'in_ad_leverage'] = in_ad[i]
        cand_df.to_csv(
            OUTPUT_DIR / "phase2b_leverage_ad.csv",
            index=False
        )
        results_2b = {
            'leverage': lev,
            'in_ad':    in_ad
        }

# Williams plot
print("\nGenerating Williams plot...")
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

for ax, results, title in [
    (axes[0], results_2a, 'Phase 2A'),
    (axes[1], results_2b, 'Phase 2B'),
]:
    if not results:
        ax.text(0.5, 0.5, 'No data',
                ha='center', va='center')
        ax.set_title(title)
        continue

    lev   = results['leverage']
    in_ad = results['in_ad']

    ax.scatter(
        lev[in_ad], np.zeros(in_ad.sum()),
        c='steelblue', alpha=0.5, s=20,
        label=f'Within AD ({in_ad.sum()})'
    )
    ax.scatter(
        lev[~in_ad], np.zeros((~in_ad).sum()),
        c='red', alpha=0.5, s=20,
        label=f'Outside AD ({(~in_ad).sum()})'
    )
    ax.axvline(h_star, color='red', linestyle='--',
               linewidth=2, label=f'h* = {h_star:.3f}')
    ax.set_xlabel('Leverage (h)', fontsize=12)
    ax.set_title(f'Leverage AD — {title}', fontsize=12)
    ax.legend(fontsize=10)

plt.tight_layout()
plt.savefig(
    OUTPUT_DIR / "leverage_ad_plot.png",
    dpi=300, bbox_inches='tight'
)

print(f"\n{'=' * 55}")
print(f"LEVERAGE AD COMPLETE")
print(f"{'=' * 55}")
print(f"  h* threshold: {h_star:.4f}")
print(f"  Plot → {OUTPUT_DIR}/leverage_ad_plot.png")
if results_2a:
    print(f"  Phase 2A → {OUTPUT_DIR}/phase2a_leverage_ad.csv")
if results_2b:
    print(f"  Phase 2B → {OUTPUT_DIR}/phase2b_leverage_ad.csv")
