"""
Prepare multi-task training data for GNN retraining.

Tasks:
  - Task 0: Mps1 pIC50 (primary target)
  - Task 1: Aurora B pIC50 (selectivity)
  - Task 2: hERG pIC50 (cardiac safety)

Compounds without activity for a task get NaN (masked loss).
"""

import pandas as pd
import numpy as np
from pathlib import Path
from rdkit import Chem
from sklearn.model_selection import train_test_split
import warnings
warnings.filterwarnings('ignore')

OUTPUT_DIR = Path("kinase_domain/analysis/gnn_model")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

IC50_DIR = Path("kinase_domain/analysis/ic50")

print("Preparing multi-task GNN training data")
print("=" * 55)

# Load all datasets
mps1 = pd.read_csv(
    IC50_DIR / "expanded_mps1_activity.csv"
)[['smiles', 'pIC50']].rename(
    columns={'pIC50': 'mps1_pIC50'}
)

aurora = pd.read_csv(
    IC50_DIR / "aurora_b_activity.csv"
)[['smiles', 'pIC50']].rename(
    columns={'pIC50': 'aurora_pIC50'}
)

herg = pd.read_csv(
    IC50_DIR / "herg_activity.csv"
)[['smiles', 'pIC50']].rename(
    columns={'pIC50': 'herg_pIC50'}
)

print(f"Mps1:    {len(mps1)} compounds")
print(f"Aurora B: {len(aurora)} compounds")
print(f"hERG:    {len(herg)} compounds")

# Canonicalise SMILES
def canonicalise(df, smi_col='smiles'):
    valid = []
    for _, row in df.iterrows():
        mol = Chem.MolFromSmiles(str(row[smi_col]))
        if mol is None:
            continue
        row[smi_col] = Chem.MolToSmiles(mol)
        valid.append(row)
    return pd.DataFrame(valid)

print("\nCanonicalising SMILES...")
mps1   = canonicalise(mps1)
aurora = canonicalise(aurora)
herg   = canonicalise(herg)

# Merge on canonical SMILES — outer join
# Compounds missing a task get NaN
merged = mps1.merge(
    aurora, on='smiles', how='outer'
).merge(
    herg, on='smiles', how='outer'
)

print(f"\nMerged dataset: {len(merged)} unique compounds")
print(f"  Mps1 data:    {merged['mps1_pIC50'].notna().sum()}")
print(f"  Aurora B data: {merged['aurora_pIC50'].notna().sum()}")
print(f"  hERG data:    {merged['herg_pIC50'].notna().sum()}")
print(f"  All 3 tasks:  {merged.dropna().shape[0]}")

# Train/val/test split
train, temp = train_test_split(
    merged, test_size=0.2, random_state=42
)
val, test = train_test_split(
    temp, test_size=0.5, random_state=42
)

print(f"\nSplit:")
print(f"  Train: {len(train)}")
print(f"  Val:   {len(val)}")
print(f"  Test:  {len(test)}")

# Save
train.to_csv(OUTPUT_DIR / "multitask_train.csv", index=False)
val.to_csv(OUTPUT_DIR / "multitask_val.csv",   index=False)
test.to_csv(OUTPUT_DIR / "multitask_test.csv",  index=False)
merged.to_csv(OUTPUT_DIR / "multitask_all.csv", index=False)

print(f"\nSaved to {OUTPUT_DIR}/multitask_*.csv")
print(f"\nTarget columns: mps1_pIC50, aurora_pIC50, herg_pIC50")
print(f"NaN = missing data for that task (masked in training)")
