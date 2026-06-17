"""
Phase 2B seed selection for PubChem similarity search.

Selects 5 maximally diverse, highly potent Mps1 inhibitors
from the ChEMBL training set to use as seeds for Phase 2B.

Criteria:
  - pIC50 >= 8.0 (IC50 <= 10 nM)
  - Lipinski Ro5 compliant
  - Max pairwise Tanimoto <= 0.4 (genuine diversity)

Output: kinase_domain/analysis/phase2/phase2b_seeds.csv
"""

import pandas as pd
import numpy as np
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, Descriptors
from rdkit.Chem import rdMolDescriptors
import warnings
warnings.filterwarnings('ignore')

CHEMBL_CSV = "kinase_domain/analysis/ic50/expanded_mps1_activity.csv"
OUTPUT_DIR = Path("kinase_domain/analysis/phase2")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

N_SEEDS    = 5
MIN_PIC50  = 8.0
MAX_SIM    = 0.4

print("Phase 2B Seed Selection")
print("=" * 55)
print(f"Min pIC50:          {MIN_PIC50}")
print(f"Max pairwise Tanimoto: {MAX_SIM}")
print(f"N seeds:            {N_SEEDS}\n")

# Load and filter ChEMBL data
df = pd.read_csv(CHEMBL_CSV)
print(f"ChEMBL compounds: {len(df)}")

results = []
for _, row in df.iterrows():
    mol = Chem.MolFromSmiles(str(row['smiles']))
    if mol is None:
        continue

    mw  = Descriptors.MolWt(mol)
    lp  = Descriptors.MolLogP(mol)
    hbd = rdMolDescriptors.CalcNumHBD(mol)
    hba = rdMolDescriptors.CalcNumHBA(mol)

    if mw > 500 or lp > 5 or hbd > 5 or hba > 10:
        continue
    if row['pIC50'] < MIN_PIC50:
        continue

    results.append({
        'smiles':   row['smiles'],
        'pIC50':    row['pIC50'],
        'IC50_nM':  row['IC50_nM'],
        'MW':       round(mw, 2),
        'LogP':     round(lp, 2),
        'source':   row['source'],
    })

pool = pd.DataFrame(results).drop_duplicates('smiles')
pool = pool.sort_values(
    'pIC50', ascending=False
).reset_index(drop=True)
print(f"Drug-like compounds with pIC50 >= {MIN_PIC50}: {len(pool)}")

# Compute fingerprints
fps       = []
valid_rows = []
for _, row in pool.iterrows():
    mol = Chem.MolFromSmiles(row['smiles'])
    if mol is None:
        continue
    fp = AllChem.GetMorganFingerprintAsBitVect(
        mol, radius=2, nBits=1024
    )
    fps.append(fp)
    valid_rows.append(row)

pool = pd.DataFrame(valid_rows).reset_index(drop=True)

# Greedy diversity selection with max similarity threshold
selected_idx = [0]
selected_fps = [fps[0]]
print(f"\nSeed 0: pIC50={pool.iloc[0]['pIC50']:.3f} "
      f"(highest potency seed)")

for candidate_idx in range(1, len(pool)):
    if len(selected_idx) >= N_SEEDS:
        break

    sims    = DataStructs.BulkTanimotoSimilarity(
        fps[candidate_idx], selected_fps
    )
    max_sim = max(sims)

    if max_sim <= MAX_SIM:
        selected_idx.append(candidate_idx)
        selected_fps.append(fps[candidate_idx])
        print(f"Seed {len(selected_idx)-1}: "
              f"pIC50={pool.iloc[candidate_idx]['pIC50']:.3f} "
              f"max_sim={max_sim:.3f} ✓")

seeds = pool.iloc[selected_idx].reset_index(drop=True)

print(f"\nFinal {len(seeds)} seeds "
      f"(max pairwise Tanimoto <= {MAX_SIM}):")
print(f"{'Seed':>5} {'pIC50':>7} {'IC50_nM':>8} "
      f"{'MW':>6} {'LogP':>5} {'Source':>12}")
print("-" * 50)
for i, row in seeds.iterrows():
    print(f"{i:>5} {row['pIC50']:>7.3f} "
          f"{row['IC50_nM']:>8.4f} "
          f"{row['MW']:>6.1f} "
          f"{row['LogP']:>5.2f} "
          f"{row['source']:>12}")

seeds.to_csv(
    OUTPUT_DIR / "phase2b_seeds.csv", index=False
)
print(f"\nSaved → {OUTPUT_DIR}/phase2b_seeds.csv")
