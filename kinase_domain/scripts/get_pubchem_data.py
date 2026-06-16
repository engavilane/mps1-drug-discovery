"""
Retrieve Mps1/TTK data from PubChem BioAssay.
Uses broader query including single-concentration screens.
"""

import pandas as pd
import numpy as np
import requests
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.inchi import MolToInchiKey
import time
import warnings
warnings.filterwarnings('ignore')

OUTPUT_DIR   = Path("kinase_domain/analysis/ic50")
EXPANDED_CSV = OUTPUT_DIR / "expanded_mps1_activity.csv"

print("Querying PubChem BioAssay for TTK...")

# Get all TTK assay IDs
aids_url = (
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/"
    "assay/target/genesymbol/TTK/aids/JSON"
)
r = requests.get(aids_url, timeout=30)
aids = r.json().get("IdentifierList", {}).get("AID", [])
print(f"  Found {len(aids)} TTK assays")

pubchem_rows = []
processed    = 0

for aid in aids:
    time.sleep(0.3)

    # Get assay description
    desc_url = (
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/"
        f"assay/aid/{aid}/description/JSON"
    )
    try:
        dr = requests.get(desc_url, timeout=20)
        if dr.status_code != 200:
            continue
        desc = dr.json().get(
            "PC_AssayContainer", [{}]
        )[0].get("assay", {}).get("descr", {})

        assay_type = desc.get("aid_type", 0)
        # Include both confirmatory (1) and dose-response (2)
        if assay_type not in [1, 2]:
            continue

        # Get active compounds
        active_url = (
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/"
            f"assay/aid/{aid}/concise/JSON"
        )
        ar = requests.get(active_url, timeout=30)
        if ar.status_code != 200:
            continue

        table = ar.json().get("Table", {})
        cols  = table.get("Columns", {}).get(
            "Column", []
        )
        rows  = table.get("Row", [])

        if not rows:
            continue

        # Find AC50 or activity score column
        ac50_idx   = None
        score_idx  = None
        active_idx = None

        for i, col in enumerate(cols):
            col_lower = col.lower()
            if "ac50" in col_lower or "ic50" in col_lower:
                ac50_idx = i
            if "activity" in col_lower:
                active_idx = i

        for row in rows:
            cells = row.get("Cell", [])
            if not cells:
                continue

            try:
                cid = int(cells[0])
            except:
                continue

            # Check activity
            if active_idx and active_idx < len(cells):
                if cells[active_idx].lower() != "active":
                    continue

            # Get AC50 if available
            pic50 = None
            if ac50_idx and ac50_idx < len(cells):
                try:
                    ac50  = float(cells[ac50_idx])
                    pic50 = -np.log10(ac50 * 1e-6)
                except:
                    pass

            # Default for confirmed actives without AC50
            if pic50 is None:
                pic50 = 6.0  # ~1 uM — conservative

            if not 3 <= pic50 <= 12:
                continue

            # Get SMILES
            time.sleep(0.1)
            smi_url = (
                f"https://pubchem.ncbi.nlm.nih.gov/rest/"
                f"pug/compound/cid/{cid}/"
                f"property/IsomericSMILES/JSON"
            )
            sr = requests.get(smi_url, timeout=10)
            if sr.status_code != 200:
                continue

            props  = sr.json().get(
                "PropertyTable", {}
            ).get("Properties", [{}])[0]
            smiles = props.get("IsomericSMILES", "")
            if not smiles:
                continue

            pubchem_rows.append({
                "smiles":        smiles,
                "pIC50":         round(pic50, 3),
                "IC50_nM":       round(10**(9-pic50), 3),
                "source":        "PubChem",
                "activity_type": "AC50",
                "chembl_id":     f"CID_{cid}",
            })

        processed += 1
        print(f"  Processed AID {aid} "
              f"({processed} assays, "
              f"{len(pubchem_rows)} compounds so far)")

        if processed >= 30:
            print("  Stopping at 30 assays "
                  "(rate limit protection)")
            break

    except Exception as e:
        continue

print(f"\nTotal PubChem records: {len(pubchem_rows)}")

if len(pubchem_rows) == 0:
    print("No PubChem records retrieved")
    exit()

# ── Merge with existing ───────────────────────────────────
existing = pd.read_csv(EXPANDED_CSV)
print(f"Existing compounds: {len(existing)}")

existing_iks = set()
for smi in existing["smiles"].dropna():
    mol = Chem.MolFromSmiles(str(smi))
    if mol:
        try:
            existing_iks.add(MolToInchiKey(mol))
        except:
            pass

new_rows = []
for row in pubchem_rows:
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

new_df   = pd.DataFrame(new_rows)
final_df = pd.concat(
    [existing, pd.DataFrame(new_rows)],
    ignore_index=True
)
final_df = final_df.sort_values(
    "pIC50", ascending=False
)
final_df.to_csv(EXPANDED_CSV, index=False)

print(f"New from PubChem:   {len(new_rows)}")
print(f"Total:              {len(final_df)}")
print(f"Saved → {EXPANDED_CSV}")
