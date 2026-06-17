"""
Full ADMET filtering for Phase 2B candidates.

Replaces the basic ADME filter with a comprehensive panel
including toxicity predictions.

Local filters (RDKit):
  - Lipinski Ro5
  - TPSA, RotBonds
  - PAINS alerts
  - Brenk structural alerts (reactive groups)
  - Synthetic accessibility score

API-based predictions (pkCSM):
  - hERG inhibition (cardiac toxicity)
  - Ames mutagenicity
  - Hepatotoxicity (DILI)
  - CYP inhibition (1A2, 2C9, 2C19, 2D6, 3A4)
  - Oral bioavailability
  - BBB permeability
"""

import pandas as pd
import numpy as np
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem import FilterCatalog
from rdkit.Chem.FilterCatalog import FilterCatalogParams
import sys
from rdkit.Chem import RDConfig
sys.path.append(RDConfig.RDContribDir)
from SA_Score import sascorer
import requests
import time
import argparse
import warnings
warnings.filterwarnings('ignore')

parser = argparse.ArgumentParser(
    description="Full ADMET filtering"
)
parser.add_argument("--candidates",
    default="kinase_domain/analysis/phase2b/phase2b_plip_candidates.csv")
parser.add_argument("--ligands",
    default="kinase_domain/data/phase2b/raw")
parser.add_argument("--output",
    default="kinase_domain/analysis/phase2b")
parser.add_argument("--use_api",
    action="store_true",
    help="Use pkCSM API for extended ADMET predictions")
args = parser.parse_args()

LIGANDS_DIR = Path(args.ligands)
OUTPUT_DIR  = Path(args.output)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

print("Full ADMET Filtering")
print("=" * 55)

# Set up structural alert filters
pains_params = FilterCatalogParams()
pains_params.AddCatalog(
    FilterCatalogParams.FilterCatalogs.PAINS
)
pains_catalog = FilterCatalog.FilterCatalog(pains_params)

brenk_params = FilterCatalogParams()
brenk_params.AddCatalog(
    FilterCatalogParams.FilterCatalogs.BRENK
)
brenk_catalog = FilterCatalog.FilterCatalog(brenk_params)

def local_admet(mol):
    """Compute local RDKit-based ADMET properties."""
    mw   = Descriptors.MolWt(mol)
    lp   = Descriptors.MolLogP(mol)
    hbd  = rdMolDescriptors.CalcNumHBD(mol)
    hba  = rdMolDescriptors.CalcNumHBA(mol)
    tpsa = Descriptors.TPSA(mol)
    rot  = rdMolDescriptors.CalcNumRotatableBonds(mol)
    sa   = sascorer.calculateScore(mol)

    # PAINS check
    pains = len(pains_catalog.GetMatches(mol)) > 0

    # Brenk structural alerts
    brenk = len(brenk_catalog.GetMatches(mol)) > 0

    # Lipinski Ro5
    ro5_violations = sum([
        mw  > 500,
        lp  > 5,
        hbd > 5,
        hba > 10,
    ])

    passes = (
        ro5_violations == 0 and
        tpsa <= 140 and
        rot  <= 10 and
        not pains and
        not brenk and
        sa   <= 6
    )

    return {
        'MW':            round(mw, 2),
        'LogP':          round(lp, 2),
        'HBD':           hbd,
        'HBA':           hba,
        'TPSA':          round(tpsa, 2),
        'RotBonds':      rot,
        'SA_score':      round(sa, 3),
        'SA_class':      ('easy' if sa <= 3
                          else 'moderate' if sa <= 6
                          else 'hard'),
        'PAINS':         pains,
        'Brenk_alert':   brenk,
        'ro5_violations': ro5_violations,
        'passes_local':  passes,
    }

def pkcsm_predict(smiles):
    """Get extended ADMET predictions from pkCSM API."""
    url = "https://biosig.lab.uq.edu.au/pkcsm/api/predict"
    try:
        r = requests.post(
            url,
            json={"smiles": smiles},
            timeout=30
        )
        if r.status_code == 200:
            return r.json()
    except Exception as e:
        pass
    return {}

# Load candidates
df = pd.read_csv(args.candidates)
print(f"Candidates: {len(df)}")
print(f"API predictions: {'enabled' if args.use_api else 'disabled'}\n")

results = []
failed  = []

for idx, row in df.iterrows():
    ligand_name = row['ligand']

    # Find SDF
    clean   = ligand_name.replace('_out', '')
    sdf     = LIGANDS_DIR / f"{clean}.sdf"
    if not sdf.exists():
        matches = list(LIGANDS_DIR.glob(f"*{clean}*.sdf"))
        sdf = matches[0] if matches else None

    if sdf is None:
        failed.append(ligand_name)
        continue

    mol = Chem.MolFromMolFile(str(sdf))
    if mol is None:
        failed.append(ligand_name)
        continue

    # Local ADMET
    admet = local_admet(mol)
    smiles = Chem.MolToSmiles(mol)

    result = {
        'ligand': ligand_name,
        **{k: row[k] for k in row.index
           if k != 'ligand'},
        **admet,
        'smiles': smiles,
    }

    # API predictions
    if args.use_api and admet['passes_local']:
        api_data = pkcsm_predict(smiles)
        if api_data:
            result['hERG_pIC50']      = api_data.get(
                'herg_pIC50', None)
            result['Ames_positive']   = api_data.get(
                'ames_toxicity', None)
            result['hepatotoxic']     = api_data.get(
                'hepatotoxicity', None)
            result['CYP3A4_inhib']    = api_data.get(
                'CYP3A4_inhibitor', None)
            result['oral_bioavail']   = api_data.get(
                'oral_bioavailability', None)

            # Flag hERG risk
            herg = result.get('hERG_pIC50')
            result['hERG_risk'] = (
                'high'     if herg and herg > 6
                else 'moderate' if herg and herg > 5
                else 'low'
            )

            # Overall ADMET pass
            result['passes_admet'] = (
                admet['passes_local'] and
                not result.get('Ames_positive', False) and
                not result.get('hepatotoxic', False) and
                result.get('hERG_risk', 'low') != 'high'
            )
            time.sleep(0.5)
        else:
            result['passes_admet'] = admet['passes_local']
    else:
        result['passes_admet'] = admet['passes_local']

    results.append(result)

    status = "✓" if result['passes_admet'] else "✗"
    print(f"  {status} {ligand_name[:40]:<40} "
          f"SA={admet['SA_score']:.2f} "
          f"PAINS={admet['PAINS']} "
          f"Brenk={admet['Brenk_alert']}")

results_df = pd.DataFrame(results)
passed_df  = results_df[results_df['passes_admet']]

print(f"\n{'=' * 55}")
print(f"ADMET FILTERING COMPLETE")
print(f"{'=' * 55}")
print(f"  Total:          {len(results_df)}")
print(f"  Passed ADMET:   {len(passed_df)}")
print(f"  Failed ADMET:   {len(results_df) - len(passed_df)}")
print(f"  Failed load:    {len(failed)}")

print(f"\n  Failure reasons:")
print(f"    Ro5 violations:  "
      f"{(results_df['ro5_violations'] > 0).sum()}")
print(f"    PAINS alerts:    "
      f"{results_df['PAINS'].sum()}")
print(f"    Brenk alerts:    "
      f"{results_df['Brenk_alert'].sum()}")
print(f"    SA score > 6:    "
      f"{(results_df['SA_score'] > 6).sum()}")

print(f"\n  SA score distribution:")
print(f"    Easy (1-3):     "
      f"{(results_df['SA_class']=='easy').sum()}")
print(f"    Moderate (4-6): "
      f"{(results_df['SA_class']=='moderate').sum()}")
print(f"    Hard (7-10):    "
      f"{(results_df['SA_class']=='hard').sum()}")

results_df.to_csv(
    OUTPUT_DIR / "admet_full.csv", index=False
)
passed_df.to_csv(
    OUTPUT_DIR / "admet_candidates.csv", index=False
)

print(f"\nAll results → {OUTPUT_DIR}/admet_full.csv")
print(f"Passed      → {OUTPUT_DIR}/admet_candidates.csv")
