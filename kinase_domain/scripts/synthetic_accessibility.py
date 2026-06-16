"""
Synthetic Accessibility (SA) scores for Phase 2 candidates.
Uses RDKit SAscore — scale 1-10:
  1-3: easy to synthesise
  4-6: moderate difficulty
  7-10: very hard to synthesise
"""

import pandas as pd
import numpy as np
from pathlib import Path
from rdkit import Chem
import sys
from rdkit.Chem import RDConfig
sys.path.append(RDConfig.RDContribDir)
from SA_Score import sascorer
import warnings
warnings.filterwarnings('ignore')

CANDIDATES_CSV = "kinase_domain/analysis/phase2/phase2_final_candidates.csv"
PHASE2_SDF     = Path("kinase_domain/data/phase2/raw")
OUTPUT_DIR     = Path("kinase_domain/analysis/phase2")

print("Computing Synthetic Accessibility scores...")
print("Scale: 1=easy, 10=very hard to synthesise\n")

df = pd.read_csv(CANDIDATES_CSV)
print(f"Candidates: {len(df)}")

sa_scores = []
failed    = []

for _, row in df.iterrows():
    name = row['ligand']

    # Find SDF
    clean = name.replace('_out', '')
    sdf   = PHASE2_SDF / f"{clean}.sdf"
    if not sdf.exists():
        matches = list(PHASE2_SDF.glob(f"*{clean}*.sdf"))
        if not matches:
            failed.append(name)
            sa_scores.append(None)
            continue
        sdf = matches[0]

    mol = Chem.MolFromMolFile(str(sdf))
    if mol is None:
        failed.append(name)
        sa_scores.append(None)
        continue

    sa = sascorer.calculateScore(mol)
    sa_scores.append(round(sa, 3))

df['sa_score'] = sa_scores

# Classify
def sa_class(score):
    if score is None:
        return 'unknown'
    elif score <= 3:
        return 'easy'
    elif score <= 6:
        return 'moderate'
    else:
        return 'hard'

df['sa_class'] = df['sa_score'].apply(sa_class)

# Summary
print(f"\nSA Score Distribution:")
print(f"  Easy (1-3):     {(df['sa_class']=='easy').sum()}")
print(f"  Moderate (4-6): {(df['sa_class']=='moderate').sum()}")
print(f"  Hard (7-10):    {(df['sa_class']=='hard').sum()}")
print(f"  Failed:         {len(failed)}")
print(f"\n  Mean SA score:  {df['sa_score'].mean():.2f}")
print(f"  Min SA score:   {df['sa_score'].min():.2f}")
print(f"  Max SA score:   {df['sa_score'].max():.2f}")

# Top 10 by combined score with SA info
print(f"\nTop 10 candidates with SA scores:")
print(f"{'Ligand':<35} {'Combined':>8} {'pIC50':>6} "
      f"{'SA':>5} {'Class':>10}")
print("-" * 70)

top10 = df.nlargest(10, 'combined_score')
for _, r in top10.iterrows():
    name = r['ligand'].replace('phase2_','')[:33]
    print(f"  {name:<33} "
          f"{r['combined_score']:>8.3f} "
          f"{r['pIC50_pred']:>6.3f} "
          f"{r['sa_score']:>5.2f} "
          f"{r['sa_class']:>10}")

# Flag hard compounds
hard = df[df['sa_class'] == 'hard']
if len(hard) > 0:
    print(f"\nCompounds flagged as hard to synthesise:")
    for _, r in hard.iterrows():
        print(f"  {r['ligand']}: SA={r['sa_score']}")

df.to_csv(OUTPUT_DIR / "phase2_final_candidates.csv", index=False)
print(f"\nUpdated → {OUTPUT_DIR}/phase2_final_candidates.csv")
