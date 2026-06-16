"""
Retrieve Mps1/TTK data from BindingDB via direct download.
UniProt: O14965 (human TTK/Mps1)
"""

import pandas as pd
import numpy as np
import requests
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.inchi import MolToInchiKey
import io
import warnings
warnings.filterwarnings('ignore')

OUTPUT_DIR   = Path("kinase_domain/analysis/ic50")
EXPANDED_CSV = OUTPUT_DIR / "expanded_mps1_activity.csv"

print("Downloading BindingDB TTK data...")

# BindingDB search by UniProt
url = (
    "https://bindingdb.org/axis2/services/"
    "BDBService/getLigandsByUniprots"
    "?uniprot=O14965"
    "&cutoff=10000000"
    "&unit=nM"
    "&response=application/json"
)

bindingdb_rows = []

try:
    r = requests.get(url, timeout=120)
    print(f"  Status: {r.status_code}")

    if r.status_code == 200:
        data  = r.json()
        affinities = data.get("affinities", [])
        print(f"  Records found: {len(affinities)}")

        for entry in affinities:
            smiles = entry.get("smiles", "")
            if not smiles:
                continue

            pic50 = None
            act_type = None

            for atype in ["IC50", "Ki", "Kd", "EC50"]:
                val_str = entry.get(atype, "")
                if not val_str:
                    continue
                try:
                    val_str = str(val_str).strip()
                    val_str = val_str.replace(">","") \
                                     .replace("<","") \
                                     .replace("~","") \
                                     .replace("=","")
                    val  = float(val_str)
                    unit = entry.get(f"{atype}_unit", "nM")
                    if unit == "nM":
                        p = -np.log10(val * 1e-9)
                    elif unit == "uM":
                        p = -np.log10(val * 1e-6)
                    elif unit == "mM":
                        p = -np.log10(val * 1e-3)
                    elif unit == "M":
                        p = -np.log10(val)
                    else:
                        continue
                    if 3 <= p <= 12:
                        pic50    = p
                        act_type = atype
                        break
                except:
                    continue

            if pic50 is None:
                continue

            bindingdb_rows.append({
                "smiles":        smiles,
                "pIC50":         round(pic50, 3),
                "IC50_nM":       round(10**(9-pic50), 3),
                "source":        "BindingDB",
                "activity_type": act_type,
                "chembl_id":     entry.get(
                    "monomerid", ""
                ),
            })

except Exception as e:
    print(f"  Error: {e}")

print(f"  Valid records: {len(bindingdb_rows)}")

if len(bindingdb_rows) == 0:
    print("\n  BindingDB API unavailable.")
    print("  Manual download instructions:")
    print("  1. Go to https://bindingdb.org/rwd/bind/"
          "byUniprot/byUniprotpage.jsp")
    print("  2. Search UniProt: O14965")
    print("  3. Download as TSV")
    print("  4. Save as kinase_domain/data/bindingdb_ttk.tsv")
    print("  Then re-run this script.")
    exit()

# ── Merge with existing expanded data ─────────────────────
existing = pd.read_csv(EXPANDED_CSV)
print(f"\nExisting compounds: {len(existing)}")

bdb_df = pd.DataFrame(bindingdb_rows)

# Canonicalise and deduplicate
from rdkit.Chem.inchi import MolToInchiKey
existing_iks = set()
for smi in existing["smiles"].dropna():
    mol = Chem.MolFromSmiles(str(smi))
    if mol:
        try:
            existing_iks.add(MolToInchiKey(mol))
        except:
            pass

new_rows = []
for _, row in bdb_df.iterrows():
    mol = Chem.MolFromSmiles(str(row["smiles"]))
    if mol is None:
        continue
    try:
        ik = MolToInchiKey(mol)
    except:
        continue
    if ik in existing_iks:
        continue
    existing_iks.add(ik)
    row["smiles"] = Chem.MolToSmiles(mol)
    row["MW"]     = round(Descriptors.MolWt(mol), 2)
    row["LogP"]   = round(Descriptors.MolLogP(mol), 2)
    new_rows.append(row)

new_df = pd.DataFrame(new_rows)
print(f"New from BindingDB: {len(new_df)}")

final_df = pd.concat(
    [existing, new_df], ignore_index=True
)
final_df = final_df.sort_values(
    "pIC50", ascending=False
)
final_df.to_csv(EXPANDED_CSV, index=False)

print(f"Total after BindingDB: {len(final_df)}")
print(f"Saved → {EXPANDED_CSV}")
