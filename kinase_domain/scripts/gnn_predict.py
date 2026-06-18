"""
GNN pIC50 prediction for Phase 2 candidates using
the optimised Chemprop D-MPNN ensemble.

Replaces SVR predictions with GNN ensemble predictions
including uncertainty estimates.
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
    description="GNN pIC50 prediction for docking candidates"
)
parser.add_argument("--candidates",
    default="kinase_domain/analysis/phase2/phase2_final_candidates.csv")
parser.add_argument("--ligands",
    default="kinase_domain/data/phase2/raw")
parser.add_argument("--model_dir",
    default="kinase_domain/analysis/gnn_model/optimised")
parser.add_argument("--output",
    default="kinase_domain/analysis/phase2")
args = parser.parse_args()

CANDIDATES_CSV = args.candidates
MODEL_DIR      = args.model_dir
LIGANDS_DIR    = Path(args.ligands)
OUTPUT_DIR     = Path(args.output)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

print("GNN pIC50 Prediction — Optimised Chemprop Ensemble")
print("=" * 55)

# Load candidates
df = pd.read_csv(CANDIDATES_CSV)
print(f"Candidates: {len(df)}")

# Extract SMILES from SDF files
print("Extracting SMILES from SDF files...")
smiles_list = []

for idx, row in df.iterrows():
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
    if mol is None:
        smiles_list.append(None)
        continue

    smiles_list.append(Chem.MolToSmiles(mol))

df['smiles'] = smiles_list
valid_df     = df[df['smiles'].notna()].copy()
print(f"  Valid SMILES: {len(valid_df)}/{len(df)}")

# Write temp input CSV for Chemprop
with tempfile.NamedTemporaryFile(
    suffix='.csv', delete=False, mode='w'
) as tmp_in:
    tmp_in_path = tmp_in.name
    valid_df[['smiles']].to_csv(tmp_in_path, index=False)

tmp_out_path = tmp_in_path.replace('.csv', '_pred.csv')

# Run each model in the ensemble
print("\nRunning GNN ensemble predictions...")
all_preds = []

for i in range(5):
    model_path = f"{MODEL_DIR}/model_{i}/best.pt"
    if not os.path.exists(model_path):
        print(f"  Model {i}: not found")
        continue

    result = subprocess.run([
        "chemprop", "predict",
        "-i",           tmp_in_path,
        "-o",           tmp_out_path,
        "--model-path", model_path,
    ], capture_output=True, text=True)

    if result.returncode != 0:
        print(f"  Model {i}: error")
        print(result.stderr[-200:])
        continue

    pred_df  = pd.read_csv(tmp_out_path)
    pred_col = [c for c in pred_df.columns
                if c != 'smiles'][0]
    all_preds.append(pred_df[pred_col].values)
    print(f"  Model {i}: {len(pred_df)} predictions ✓")

# Compute ensemble mean and uncertainty
print(f"\nComputing ensemble statistics...")
preds_array   = np.array(all_preds)
ensemble_mean = preds_array.mean(axis=0)
ensemble_std  = preds_array.std(axis=0)

valid_df['pIC50_gnn']      = ensemble_mean.round(3)
valid_df['pIC50_gnn_std']  = ensemble_std.round(3)
valid_df['IC50_gnn_nM']    = (10 ** (9 - ensemble_mean)).round(2)
valid_df['gnn_confidence'] = valid_df['pIC50_gnn_std'].apply(
    lambda s: 'high'     if s < 0.1
    else ('moderate' if s < 0.2
    else 'low')
)

# Recompute combined score using GNN pIC50
print("Recomputing combined score with GNN predictions...")
valid_df['vina_norm'] = (
    valid_df['affinity_kcal'].min() -
    valid_df['affinity_kcal']
) / (
    valid_df['affinity_kcal'].min() -
    valid_df['affinity_kcal'].max()
)
valid_df['pIC50_norm'] = (
    valid_df['pIC50_gnn'] -
    valid_df['pIC50_gnn'].min()
) / (
    valid_df['pIC50_gnn'].max() -
    valid_df['pIC50_gnn'].min()
)
valid_df['combined_score_gnn'] = (
    valid_df['vina_norm'] +
    valid_df['pIC50_norm']
) / 2

# Print summary
print(f"\n{'=' * 55}")
print(f"GNN PREDICTION RESULTS")
print(f"{'=' * 55}")
print(f"  Compounds predicted:  {len(valid_df)}")
print(f"  pIC50 range:          "
      f"{valid_df['pIC50_gnn'].min():.2f} - "
      f"{valid_df['pIC50_gnn'].max():.2f}")
print(f"  Best IC50:            "
      f"{valid_df['IC50_gnn_nM'].min():.2f} nM")
print(f"  Mean uncertainty:     "
      f"{valid_df['pIC50_gnn_std'].mean():.3f}")
print(f"\n  Confidence breakdown:")
print(f"    High (σ < 0.1):     "
      f"{(valid_df['gnn_confidence']=='high').sum()}")
print(f"    Moderate (σ < 0.2): "
      f"{(valid_df['gnn_confidence']=='moderate').sum()}")
print(f"    Low (σ ≥ 0.2):      "
      f"{(valid_df['gnn_confidence']=='low').sum()}")

print(f"\nTop 10 by GNN combined score:")
print(f"{'Ligand':<32} {'Vina':>7} {'pIC50':>6} "
      f"{'±':>5} {'IC50(nM)':>9} {'Conf':>8}")
print("-" * 72)

for _, r in valid_df.nlargest(10, 'combined_score_gnn').iterrows():
    name = r['ligand'].replace('phase2_', '')[:30]
    print(f"  {name:<30} "
          f"{r['affinity_kcal']:>7.3f} "
          f"{r['pIC50_gnn']:>6.3f} "
          f"{r['pIC50_gnn_std']:>5.3f} "
          f"{r['IC50_gnn_nM']:>9.2f} "
          f"{r['gnn_confidence']:>8}")

print(f"\nSVR vs GNN — top 10 by GNN pIC50:")
print(f"{'Ligand':<32} {'SVR':>7} {'GNN':>7} {'Δ':>6}")
print("-" * 55)

for _, r in valid_df.nlargest(10, 'pIC50_gnn').iterrows():
    name = r['ligand'].replace('phase2_', '')[:30]
    diff = r['pIC50_gnn'] - r['pIC50_pred']
    print(f"  {name:<30} "
          f"{r['pIC50_pred']:>7.3f} "
          f"{r['pIC50_gnn']:>7.3f} "
          f"{diff:>+6.3f}")

# Save results
valid_df.to_csv(
    OUTPUT_DIR / "phase2_final_candidates_gnn.csv",
    index=False
)
print(f"\nSaved → {OUTPUT_DIR}/phase2_final_candidates_gnn.csv")

# Cleanup temp files
os.unlink(tmp_in_path)
if os.path.exists(tmp_out_path):
    os.unlink(tmp_out_path)
