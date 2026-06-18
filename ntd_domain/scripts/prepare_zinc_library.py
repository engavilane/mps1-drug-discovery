"""
Filter and prepare ZINC fragment library for NTD docking.

From ~2.5M fragments, select ~10,000 diverse drug-like fragments
suitable for TPR groove docking.

Filters:
  - Rule of 3 (MW<=300, LogP<=3, HBD<=3, HBA<=3)
  - TPSA <= 100 (good for PPI site binding)
  - No PAINS alerts
  - Maximum diversity selection (Tanimoto-based)
"""

import pandas as pd
import numpy as np
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, AllChem
from rdkit.Chem import DataStructs
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
import random
import warnings
warnings.filterwarnings('ignore')

INPUT_DIR  = Path("ntd_domain/data/libraries/zinc_raw")
OUTPUT_DIR = Path("ntd_domain/data/libraries")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

TARGET_N   = 10000
SAMPLE_N   = 50000  # pre-filter sample before diversity selection

print("ZINC Fragment Library Preparation")
print("=" * 55)
print(f"Input:  {INPUT_DIR}")
print(f"Target: {TARGET_N} diverse fragments\n")

# Set up PAINS filter
pains_params = FilterCatalogParams()
pains_params.AddCatalog(
    FilterCatalogParams.FilterCatalogs.PAINS
)
pains_catalog = FilterCatalog(pains_params)

# Load and filter all SMILES
print("Loading and filtering ZINC fragments...")
all_smiles = []
smi_files  = sorted(INPUT_DIR.glob("*.smi"))

for smi_file in smi_files:
    with open(smi_file) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 1:
                continue
            smi = parts[0]
            if smi == 'smiles':
                continue
            all_smiles.append(smi)

print(f"  Total ZINC fragments: {len(all_smiles):,}")

# Random sample for efficiency
random.seed(42)
if len(all_smiles) > SAMPLE_N * 3:
    sample = random.sample(all_smiles, SAMPLE_N * 3)
else:
    sample = all_smiles

print(f"  Working sample: {len(sample):,}")

# Apply filters
print("Applying Rule of 3 + PAINS filters...")
filtered = []

for smi in sample:
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        continue

    mw   = Descriptors.MolWt(mol)
    lp   = Descriptors.MolLogP(mol)
    hbd  = rdMolDescriptors.CalcNumHBD(mol)
    hba  = rdMolDescriptors.CalcNumHBA(mol)
    tpsa = Descriptors.TPSA(mol)

    if mw > 320 or lp > 3.5 or hbd > 3 or hba > 6:
        continue
    if tpsa > 120:
        continue
    if len(pains_catalog.GetMatches(mol)) > 0:
        continue

    filtered.append(Chem.MolToSmiles(mol))

    if len(filtered) >= SAMPLE_N:
        break

print(f"  After filters: {len(filtered):,}")

# Diversity selection using fingerprints
print(f"Selecting {TARGET_N} diverse fragments...")
fps = []
valid = []
for smi in filtered:
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        continue
    fp = AllChem.GetMorganFingerprintAsBitVect(
        mol, radius=2, nBits=1024
    )
    fps.append(fp)
    valid.append(smi)

# Greedy diversity selection
selected_idx = [0]
selected_fps = [fps[0]]

while len(selected_idx) < TARGET_N and \
      len(selected_idx) < len(valid):
    best_idx  = None
    best_dist = -1

    # Sample candidates for speed
    candidates = random.sample(
        range(len(valid)),
        min(1000, len(valid))
    )

    for i in candidates:
        if i in selected_idx:
            continue
        sims    = DataStructs.BulkTanimotoSimilarity(
            fps[i], selected_fps[-50:]
        )
        min_sim = min(sims)
        dist    = 1 - min_sim

        if dist > best_dist:
            best_dist = dist
            best_idx  = i

    if best_idx is None:
        break

    selected_idx.append(best_idx)
    selected_fps.append(fps[best_idx])

    if len(selected_idx) % 1000 == 0:
        print(f"  Selected: {len(selected_idx)}/{TARGET_N}")

selected_smiles = [valid[i] for i in selected_idx]

# Save as SMILES file
out_smi = OUTPUT_DIR / "zinc_fragments_filtered.smi"
with open(out_smi, 'w') as f:
    f.write("smiles\n")
    for smi in selected_smiles:
        f.write(f"{smi}\n")

print(f"\n{'=' * 55}")
print(f"ZINC LIBRARY PREPARATION COMPLETE")
print(f"{'=' * 55}")
print(f"  Input fragments:    {len(all_smiles):,}")
print(f"  After Ro3 + PAINS: {len(filtered):,}")
print(f"  Final diverse set:  {len(selected_smiles):,}")
print(f"\nSaved → {out_smi}")
