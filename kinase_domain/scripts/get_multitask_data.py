"""
Collect multi-task training data from ChEMBL:
  - Mps1/TTK inhibitors (already have)
  - Aurora B (AURKB) inhibitors
  - hERG (KCNH2) inhibitors

For multi-task GNN training to expand applicability domain.
"""

import pandas as pd
import requests
import time
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
import warnings
warnings.filterwarnings('ignore')

OUTPUT_DIR = Path("kinase_domain/analysis/ic50")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

def get_chembl_data(target_id, target_name, max_records=5000):
    print(f"\nFetching {target_name} ({target_id})...")

    url = "https://www.ebi.ac.uk/chembl/api/data/activity"
    params = {
        'target_chembl_id': target_id,
        'standard_type__in': 'IC50,Ki,Kd',
        'standard_relation': '=',
        'standard_units': 'nM',
        'pchembl_value__isnull': False,
        'assay_type': 'B',
        'format': 'json',
        'limit': 1000,
        'offset': 0,
    }

    all_records = []
    while len(all_records) < max_records:
        r = requests.get(url, params=params, timeout=30)
        if r.status_code != 200:
            print(f"  Error: {r.status_code}")
            break

        data = r.json()
        records = data.get('activities', [])
        if not records:
            break

        all_records.extend(records)
        total = data.get('page_meta', {}).get('total_count', 0)
        print(f"  Retrieved: {len(all_records)}/{total}")

        if len(all_records) >= total:
            break
        params['offset'] += 1000
        time.sleep(0.3)

    # Process records
    results = []
    for rec in all_records:
        smiles = rec.get('canonical_smiles', '')
        pic50  = rec.get('pchembl_value')
        if not smiles or not pic50:
            continue

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue

        mw  = Descriptors.MolWt(mol)
        lp  = Descriptors.MolLogP(mol)
        hbd = rdMolDescriptors.CalcNumHBD(mol)
        hba = rdMolDescriptors.CalcNumHBA(mol)

        if mw > 800 or lp > 8:
            continue

        results.append({
            'smiles':  smiles,
            'pIC50':   float(pic50),
            'IC50_nM': 10 ** (9 - float(pic50)),
            'target':  target_name,
            'source':  'ChEMBL',
        })

    df = pd.DataFrame(results).drop_duplicates('smiles')
    print(f"  Valid compounds: {len(df)}")
    return df

# Fetch data for all targets
targets = {
    'CHEMBL2056': 'Aurora_B',
    'CHEMBL240':  'hERG',
}

all_data = {}
for chembl_id, name in targets.items():
    df = get_chembl_data(chembl_id, name)
    all_data[name] = df
    df.to_csv(
        OUTPUT_DIR / f"{name.lower()}_activity.csv",
        index=False
    )
    print(f"  Saved → {name.lower()}_activity.csv")

# Summary
print(f"\n{'='*55}")
print(f"MULTI-TASK DATA COLLECTION COMPLETE")
print(f"{'='*55}")
for name, df in all_data.items():
    print(f"  {name}: {len(df)} compounds")
    print(f"    pIC50 range: {df['pIC50'].min():.2f} - "
          f"{df['pIC50'].max():.2f}")
