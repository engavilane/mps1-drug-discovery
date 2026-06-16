"""
Applicability Domain Analysis for Model 2 (ChEMBL SVR).

For each Phase 2 candidate, computes the Tanimoto similarity
to its nearest neighbour in the ChEMBL training set.

Compounds outside the applicability domain (Tanimoto < threshold)
have unreliable pIC50 predictions.

Standard threshold: Tanimoto ≥ 0.3 (Tropsha et al. 2003)
Strict threshold:   Tanimoto ≥ 0.4
"""

import pandas as pd
import numpy as np
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')


# Paths 
CHEMBL_CSV    = "kinase_domain/analysis/ic50/chembl_mps1_ic50.csv"
CANDIDATES_CSV = "kinase_domain/analysis/phase2/phase2_final_candidates.csv"
PHASE2_SDF    = Path("kinase_domain/data/phase2/raw")
OUTPUT_DIR    = Path("kinase_domain/analysis/ml_model/chembl")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# Thresholds
STANDARD_THRESHOLD = 0.3
STRICT_THRESHOLD   = 0.4


# Load training set
print("Loading ChEMBL training set...")
train_df = pd.read_csv(CHEMBL_CSV)
train_df = train_df.dropna(subset=["canonical_smiles", "pIC50"])
train_df = train_df[train_df["IC50_nM"] > 0]
train_df = train_df.drop_duplicates(subset="canonical_smiles")
print(f"  {len(train_df)} training compounds")

# Compute training set fingerprints
print("  Computing training fingerprints...")
train_fps = []
valid_train = []
for idx, row in train_df.iterrows():
    mol = Chem.MolFromSmiles(row["canonical_smiles"])
    if mol is None:
        continue
    fp = AllChem.GetMorganFingerprintAsBitVect(
        mol, radius=2, nBits=1024
    )
    train_fps.append(fp)
    valid_train.append(idx)

train_df = train_df.loc[valid_train].reset_index(drop=True)
print(f"  {len(train_fps)} valid training fingerprints\n")


# Load Phase 2 candidates
print("Loading Phase 2 candidates...")
cand_df = pd.read_csv(CANDIDATES_CSV)
print(f"  {len(cand_df)} candidates\n")


# Compute Tanimoto to nearest neighbour
print("Computing applicability domain...")
print("(Tanimoto similarity to nearest training compound)")
print("=" * 55)

results = []
outside_standard = []
outside_strict   = []

for idx, row in cand_df.iterrows():
    ligand_name = row["ligand"]

    # Find SDF file
    clean_name = ligand_name.replace("_out", "")
    sdf_path   = PHASE2_SDF / f"{clean_name}.sdf"

    if not sdf_path.exists():
        matches = list(PHASE2_SDF.glob(f"*{clean_name}*.sdf"))
        if not matches:
            continue
        sdf_path = matches[0]

    mol = Chem.MolFromMolFile(str(sdf_path))
    if mol is None:
        continue

    # Compute fingerprint
    fp = AllChem.GetMorganFingerprintAsBitVect(
        mol, radius=2, nBits=1024
    )

    # Tanimoto to all training compounds
    similarities = DataStructs.BulkTanimotoSimilarity(
        fp, train_fps
    )
    similarities = np.array(similarities)

    max_sim  = similarities.max()
    mean_sim = similarities.mean()
    nn_idx   = similarities.argmax()
    nn_smiles = train_df.iloc[nn_idx]["canonical_smiles"]
    nn_pic50  = train_df.iloc[nn_idx]["pIC50"]

    in_standard = max_sim >= STANDARD_THRESHOLD
    in_strict   = max_sim >= STRICT_THRESHOLD

    results.append({
        "ligand":          ligand_name,
        "affinity_kcal":   row["affinity_kcal"],
        "pIC50_pred":      row["pIC50_pred"],
        "IC50_pred_nM":    row["IC50_pred_nM"],
        "combined_score":  row.get("combined_score", None),
        "max_tanimoto":    round(max_sim, 4),
        "mean_tanimoto":   round(mean_sim, 4),
        "nn_pIC50":        round(nn_pic50, 3),
        "in_ad_standard":  in_standard,
        "in_ad_strict":    in_strict,
    })

    if not in_standard:
        outside_standard.append(ligand_name)
    if not in_strict:
        outside_strict.append(ligand_name)


# Results 
results_df = pd.DataFrame(results)
results_df = results_df.sort_values(
    "max_tanimoto", ascending=False
)

n_total    = len(results_df)
n_standard = results_df["in_ad_standard"].sum()
n_strict   = results_df["in_ad_strict"].sum()

print(f"\n{'=' * 55}")
print(f"APPLICABILITY DOMAIN RESULTS")
print(f"{'=' * 55}")
print(f"Total candidates:              {n_total}")
print(f"Within AD (Tanimoto ≥ 0.3):   "
      f"{n_standard}/{n_total} "
      f"({100*n_standard/n_total:.1f}%)")
print(f"Within AD (Tanimoto ≥ 0.4):   "
      f"{n_strict}/{n_total} "
      f"({100*n_strict/n_total:.1f}%)")
print(f"\nTanimoto statistics:")
print(f"  Mean:   {results_df['max_tanimoto'].mean():.3f}")
print(f"  Std:    {results_df['max_tanimoto'].std():.3f}")
print(f"  Min:    {results_df['max_tanimoto'].min():.3f}")
print(f"  Max:    {results_df['max_tanimoto'].max():.3f}")
print(f"  Median: {results_df['max_tanimoto'].median():.3f}")

if outside_standard:
    print(f"\nCompounds outside AD (Tanimoto < 0.3):")
    for l in outside_standard:
        row = results_df[results_df["ligand"]==l].iloc[0]
        print(f"  {l}: Tanimoto={row['max_tanimoto']:.3f} "
              f"pIC50_pred={row['pIC50_pred']:.3f}")
else:
    print(f"\n✓ All candidates within standard AD (≥ 0.3)")

if outside_strict:
    print(f"\nCompounds outside strict AD (Tanimoto < 0.4):")
    for l in outside_strict[:10]:
        row = results_df[results_df["ligand"]==l].iloc[0]
        print(f"  {l}: Tanimoto={row['max_tanimoto']:.3f} "
              f"pIC50_pred={row['pIC50_pred']:.3f}")
    if len(outside_strict) > 10:
        print(f"  ... and {len(outside_strict)-10} more")


# Top 10 candidates with AD info
print(f"\nTop 10 candidates by combined score + AD status:")
print(f"{'Ligand':<32} {'Tanimoto':>9} {'pIC50':>6} "
      f"{'AD(0.3)':>8} {'AD(0.4)':>8}")
print("-" * 65)

if "combined_score" in results_df.columns:
    results_df["combined_score"] = pd.to_numeric(
        results_df["combined_score"], errors="coerce"
    )
    top10 = results_df.nlargest(10, "combined_score")
else:
    top10 = results_df.head(10)

for _, r in top10.iterrows():
    name = r["ligand"].replace("phase2_", "")[:30]
    ad3  = "✓" if r["in_ad_standard"] else "✗"
    ad4  = "✓" if r["in_ad_strict"]   else "✗"
    print(f"  {name:<30} "
          f"{r['max_tanimoto']:>9.3f} "
          f"{r['pIC50_pred']:>6.3f} "
          f"{ad3:>8} "
          f"{ad4:>8}")


# Save 
results_df.to_csv(
    OUTPUT_DIR / "applicability_domain.csv", index=False
)


# Plot 
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Plot 1 — Tanimoto distribution
axes[0].hist(
    results_df["max_tanimoto"], bins=30,
    color='steelblue', alpha=0.7, edgecolor='white'
)
axes[0].axvline(
    STANDARD_THRESHOLD, color='orange',
    linewidth=2, linestyle='--',
    label=f'Standard threshold ({STANDARD_THRESHOLD})'
)
axes[0].axvline(
    STRICT_THRESHOLD, color='red',
    linewidth=2, linestyle='--',
    label=f'Strict threshold ({STRICT_THRESHOLD})'
)
axes[0].set_xlabel('Max Tanimoto to training set', fontsize=12)
axes[0].set_ylabel('Count', fontsize=12)
axes[0].set_title('Applicability Domain\n(Tanimoto to nearest training compound)',
                   fontsize=12)
axes[0].legend(fontsize=10)

# Plot 2 — Tanimoto vs predicted pIC50
colours = results_df["in_ad_strict"].map(
    {True: 'steelblue', False: 'red'}
)
axes[1].scatter(
    results_df["max_tanimoto"],
    results_df["pIC50_pred"],
    c=colours, alpha=0.6, s=20
)
axes[1].axvline(
    STRICT_THRESHOLD, color='red',
    linewidth=1.5, linestyle='--',
    label='Strict AD threshold (0.4)'
)
axes[1].set_xlabel('Max Tanimoto to training set', fontsize=12)
axes[1].set_ylabel('Predicted pIC50', fontsize=12)
axes[1].set_title('Predicted pIC50 vs Applicability Domain',
                   fontsize=12)
axes[1].legend(fontsize=10)

# Highlight top candidate
top_cand = results_df[
    results_df["ligand"].str.contains("142416385")
]
if len(top_cand) > 0:
    axes[1].scatter(
        top_cand["max_tanimoto"],
        top_cand["pIC50_pred"],
        c='gold', s=100, zorder=5,
        label='CID 142416385', marker='*'
    )
    axes[1].legend(fontsize=10)

plt.tight_layout()
plt.savefig(
    OUTPUT_DIR / "applicability_domain_plot.png",
    dpi=300, bbox_inches='tight'
)

print(f"\nResults → {OUTPUT_DIR}/applicability_domain.csv")
print(f"Plot    → {OUTPUT_DIR}/applicability_domain_plot.png")
