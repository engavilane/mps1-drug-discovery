"""
Multi-task GNN predictions for Phase 2B candidates.

Predicts simultaneously:
  - Mps1 pIC50 (primary target)
  - Aurora B pIC50 (selectivity)
  - hERG pIC50 (cardiac safety)

Computes:
  - Selectivity index: Aurora B IC50 / Mps1 IC50
  - Safety flag: hERG pIC50 > 6 = cardiac risk
"""

import pandas as pd
import numpy as np
from pathlib import Path
from rdkit import Chem
import subprocess
import tempfile
import os
import argparse
import warnings
warnings.filterwarnings('ignore')

parser = argparse.ArgumentParser(
    description="Multi-task GNN predictions"
)
parser.add_argument("--candidates",
    default="kinase_domain/analysis/phase2b/admet_candidates_clean.csv")
parser.add_argument("--ligands",
    default="kinase_domain/data/phase2b/raw")
parser.add_argument("--model_dir",
    default="kinase_domain/analysis/gnn_model/multitask")
parser.add_argument("--output",
    default="kinase_domain/analysis/phase2b")
args = parser.parse_args()

LIGANDS_DIR = Path(args.ligands)
OUTPUT_DIR  = Path(args.output)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

print("Multi-task GNN Predictions")
print("=" * 55)
print("Tasks: Mps1 pIC50 | Aurora B pIC50 | hERG pIC50")

df = pd.read_csv(args.candidates)
print(f"Candidates: {len(df)}")

# Extract SMILES
print("Extracting SMILES...")
smiles_list = []
for _, row in df.iterrows():
    name  = row['ligand']
    clean = name.replace('_out', '')
    sdf   = LIGANDS_DIR / f"{clean}.sdf"
    if not sdf.exists():
        matches = list(LIGANDS_DIR.glob(f"*{clean}*.sdf"))
        sdf = matches[0] if matches else None
    if sdf is None:
        smiles_list.append(None)
        continue
    mol = Chem.MolFromMolFile(str(sdf))
    smiles_list.append(
        Chem.MolToSmiles(mol) if mol else None
    )

df['smiles'] = smiles_list
valid_df     = df[df['smiles'].notna()].copy()
print(f"Valid SMILES: {len(valid_df)}/{len(df)}")

with tempfile.NamedTemporaryFile(
    suffix='.csv', delete=False, mode='w'
) as tmp:
    tmp_path = tmp.name
    valid_df[['smiles']].to_csv(tmp_path, index=False)

tmp_out = tmp_path.replace('.csv', '_pred.csv')

print("\nRunning multi-task ensemble predictions...")
tasks = ['mps1_pIC50', 'aurora_pIC50', 'herg_pIC50']
all_preds = {task: [] for task in tasks}

for i in range(5):
    model_path = f"{args.model_dir}/model_{i}/best.pt"
    if not os.path.exists(model_path):
        print(f"  Model {i}: not found")
        continue

    result = subprocess.run([
        'chemprop', 'predict',
        '-i',           tmp_path,
        '-o',           tmp_out,
        '--model-path', model_path,
    ], capture_output=True, text=True)

    if result.returncode != 0:
        print(f"  Model {i}: error")
        print(result.stderr[-200:])
        continue

    pred_df  = pd.read_csv(tmp_out)
    pred_cols = [c for c in pred_df.columns if c != 'smiles']
    print(f"  Model {i}: {len(pred_df)} predictions "
          f"({pred_cols})")

    for j, task in enumerate(tasks):
        if j < len(pred_cols):
            all_preds[task].append(
                pred_df[pred_cols[j]].values
            )

print("\nComputing ensemble statistics...")
for task in tasks:
    if all_preds[task]:
        preds = np.array(all_preds[task])
        valid_df[f'{task}_mean'] = preds.mean(axis=0).round(3)
        valid_df[f'{task}_std']  = preds.std(axis=0).round(3)

valid_df['IC50_mps1_nM']   = (
    10 ** (9 - valid_df['mps1_pIC50_mean'])
).round(2)
valid_df['IC50_aurora_nM'] = (
    10 ** (9 - valid_df['aurora_pIC50_mean'])
).round(2)
valid_df['IC50_herg_nM']   = (
    10 ** (9 - valid_df['herg_pIC50_mean'])
).round(2)

valid_df['selectivity_index'] = (
    valid_df['IC50_aurora_nM'] /
    valid_df['IC50_mps1_nM']
).round(2)

valid_df['herg_risk'] = valid_df['herg_pIC50_mean'].apply(
    lambda x: 'high'     if x > 6
    else ('moderate' if x > 5
    else 'low')
)
valid_df['selective'] = valid_df['selectivity_index'] > 100

print(f"\n{'=' * 55}")
print(f"MULTI-TASK PREDICTION RESULTS")
print(f"{'=' * 55}")
print(f"  Compounds predicted: {len(valid_df)}")
print(f"\n  Mps1 pIC50:   "
      f"{valid_df['mps1_pIC50_mean'].min():.2f} - "
      f"{valid_df['mps1_pIC50_mean'].max():.2f}")
print(f"  Aurora pIC50: "
      f"{valid_df['aurora_pIC50_mean'].min():.2f} - "
      f"{valid_df['aurora_pIC50_mean'].max():.2f}")
print(f"  hERG pIC50:   "
      f"{valid_df['herg_pIC50_mean'].min():.2f} - "
      f"{valid_df['herg_pIC50_mean'].max():.2f}")

print(f"\n  Selective (SI > 100): "
      f"{valid_df['selective'].sum()}/{len(valid_df)}")
print(f"  hERG risk high:       "
      f"{(valid_df['herg_risk']=='high').sum()}")
print(f"  hERG risk moderate:   "
      f"{(valid_df['herg_risk']=='moderate').sum()}")
print(f"  hERG risk low:        "
      f"{(valid_df['herg_risk']=='low').sum()}")

print(f"\nTop 15 by Mps1 pIC50:")
print(f"{'Ligand':<33} {'Mps1':>6} {'AurB':>6} "
      f"{'hERG':>6} {'SI':>8} {'Risk':>8} {'Sel':>4}")
print("-" * 75)

for _, r in valid_df.nlargest(
    15, 'mps1_pIC50_mean'
).iterrows():
    name = r['ligand'].replace('phase2b_', '')[:31]
    sel  = "✓" if r['selective'] else "✗"
    print(f"  {name:<31} "
          f"{r['mps1_pIC50_mean']:>6.3f} "
          f"{r['aurora_pIC50_mean']:>6.3f} "
          f"{r['herg_pIC50_mean']:>6.3f} "
          f"{r['selectivity_index']:>8.1f} "
          f"{r['herg_risk']:>8} "
          f"{sel:>4}")

valid_df.to_csv(
    OUTPUT_DIR / "phase2b_multitask_predictions.csv",
    index=False
)
print(f"\nSaved → "
      f"{OUTPUT_DIR}/phase2b_multitask_predictions.csv")

os.unlink(tmp_path)
if os.path.exists(tmp_out):
    os.unlink(tmp_out)
