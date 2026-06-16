"""
Training data expansion for Mps1/TTK QSAR model.

Sources:
  1. ChEMBL — broader query including Ki, Kd, percent inhibition
  2. BindingDB — Ki and IC50 values
  3. PubChem BioAssay — HTS data

Target: 5,000-10,000 unique compounds with pIC50 equivalent
"""

import pandas as pd
import numpy as np
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from chembl_webresource_client.new_client import new_client
import requests
import time
import warnings
warnings.filterwarnings('ignore')

# ── Paths ─────────────────────────────────────────────────
OUTPUT_DIR = Path("kinase_domain/analysis/ic50")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

EXISTING_CSV = OUTPUT_DIR / "chembl_mps1_ic50.csv"
EXPANDED_CSV = OUTPUT_DIR / "expanded_mps1_activity.csv"

# ── Constants ─────────────────────────────────────────────
# ChEMBL target ID for Mps1/TTK
CHEMBL_TARGET = "CHEMBL4523"

# Conversion factors
# pIC50 = -log10(IC50 in M)
# pKi   = -log10(Ki in M)
# These are comparable for competitive inhibitors

# ── Source 1 — ChEMBL broader query ──────────────────────
print("=" * 60)
print("SOURCE 1 — ChEMBL (expanded query)")
print("=" * 60)

activity  = new_client.activity
molecule  = new_client.molecule

print("Querying ChEMBL for all activity types...")

# Query all activity types not just IC50
all_activities = activity.filter(
    target_chembl_id=CHEMBL_TARGET,
    assay_type="B"  # Binding assays
).only([
    "molecule_chembl_id",
    "canonical_smiles",
    "standard_type",
    "standard_value",
    "standard_units",
    "standard_relation",
    "pchembl_value",
    "assay_chembl_id",
    "assay_description",
])

chembl_rows = []
for act in all_activities:
    # Skip if no SMILES or value
    if not act.get("canonical_smiles"):
        continue
    if not act.get("standard_value"):
        continue

    std_type  = act.get("standard_type", "")
    std_value = act.get("standard_value")
    std_units = act.get("standard_units", "")
    pchembl   = act.get("pchembl_value")

    # Convert to pIC50 equivalent
    pic50 = None

    # Use pChEMBL value directly if available
    if pchembl:
        try:
            pic50 = float(pchembl)
        except:
            pass

    # Otherwise convert from standard value
    if pic50 is None and std_value and std_units:
        try:
            val = float(std_value)
            # Convert nM to pIC50
            if std_units == "nM":
                pic50 = -np.log10(val * 1e-9)
            elif std_units == "uM":
                pic50 = -np.log10(val * 1e-6)
            elif std_units == "mM":
                pic50 = -np.log10(val * 1e-3)
            elif std_units == "M":
                pic50 = -np.log10(val)
        except:
            pass

    if pic50 is None:
        continue

    # Filter reasonable range
    if pic50 < 3 or pic50 > 12:
        continue

    chembl_rows.append({
        "smiles":      act["canonical_smiles"],
        "pIC50":       round(pic50, 3),
        "IC50_nM":     10**(9 - pic50),
        "source":      "ChEMBL",
        "activity_type": std_type,
        "chembl_id":   act.get("molecule_chembl_id", ""),
    })

chembl_df = pd.DataFrame(chembl_rows)
print(f"  Raw ChEMBL records: {len(chembl_df)}")
print(f"  Activity types: {chembl_df['activity_type'].value_counts().to_dict()}")

# ── Source 2 — BindingDB ──────────────────────────────────
print(f"\n{'=' * 60}")
print("SOURCE 2 — BindingDB")
print("=" * 60)

print("Querying BindingDB for Mps1/TTK...")

# BindingDB REST API
# UniProt ID for human Mps1/TTK: O14965
UNIPROT_ID = "O14965"

bindingdb_url = (
    f"https://bindingdb.org/axis2/services/BDBService/"
    f"getLigandsByUniprots"
    f"?uniprot={UNIPROT_ID}"
    f"&response=application/json"
)

bindingdb_rows = []
try:
    r = requests.get(bindingdb_url, timeout=60)
    if r.status_code == 200:
        data = r.json()
        affinities = data.get(
            "affinities", []
        )
        print(f"  Raw BindingDB records: {len(affinities)}")

        for entry in affinities:
            smiles = entry.get("smiles", "")
            if not smiles:
                continue

            # Try IC50, Ki, Kd in order
            pic50 = None
            for affinity_type in ["IC50", "Ki", "Kd"]:
                val_str = entry.get(affinity_type, "")
                if not val_str or val_str in ["", "None"]:
                    continue
                try:
                    # Handle ranges like ">10000" or "<1"
                    val_str = val_str.replace(">", "")
                    val_str = val_str.replace("<", "")
                    val_str = val_str.replace("~", "")
                    val     = float(val_str)
                    units   = entry.get(
                        f"{affinity_type}_unit", "nM"
                    )
                    if units == "nM":
                        pic50 = -np.log10(val * 1e-9)
                    elif units == "uM":
                        pic50 = -np.log10(val * 1e-6)
                    elif units == "mM":
                        pic50 = -np.log10(val * 1e-3)
                    elif units == "M":
                        pic50 = -np.log10(val)
                    if pic50 and 3 <= pic50 <= 12:
                        break
                except:
                    continue

            if pic50 is None:
                continue

            bindingdb_rows.append({
                "smiles":        smiles,
                "pIC50":         round(pic50, 3),
                "IC50_nM":       10**(9 - pic50),
                "source":        "BindingDB",
                "activity_type": affinity_type,
                "chembl_id":     "",
            })

    else:
        print(f"  BindingDB API error: {r.status_code}")

except Exception as e:
    print(f"  BindingDB query failed: {e}")
    print("  Trying alternative BindingDB endpoint...")

    # Alternative: download TSV directly
    alt_url = (
        "https://bindingdb.org/bind/downloads/"
        "BindingDB_TTK_2024m9.tsv.zip"
    )
    print(f"  Manual download required from:")
    print(f"  https://bindingdb.org/rwd/bind/byUniprot/"
          f"byUniprotpage.jsp?tab=protein"
          f"&uniprot={UNIPROT_ID}")

bindingdb_df = pd.DataFrame(bindingdb_rows)
print(f"  Valid BindingDB records: {len(bindingdb_df)}")

# ── Source 3 — PubChem BioAssay ───────────────────────────
print(f"\n{'=' * 60}")
print("SOURCE 3 — PubChem BioAssay")
print("=" * 60)

print("Querying PubChem BioAssay for TTK/Mps1...")

pubchem_rows = []
try:
    # Search for TTK assays
    search_url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/"
        "assay/target/genesymbol/TTK/aids/JSON"
    )
    r = requests.get(search_url, timeout=30)

    if r.status_code == 200:
        aids = r.json().get("IdentifierList", {}).get(
            "AID", []
        )
        print(f"  Found {len(aids)} TTK assays")

        # Get data from most relevant assays
        # Focus on dose-response assays (type = 2)
        for aid in aids[:20]:  # Limit to first 20
            time.sleep(0.2)  # Rate limiting
            assay_url = (
                f"https://pubchem.ncbi.nlm.nih.gov/"
                f"rest/pug/assay/aid/{aid}/JSON"
            )
            ar = requests.get(assay_url, timeout=30)
            if ar.status_code != 200:
                continue

            assay_data = ar.json()
            assay_info = assay_data.get(
                "PC_AssayContainer", [{}]
            )[0].get("assay", {}).get("descr", {})

            assay_type = assay_info.get("aid_type", 0)

            # Only dose-response assays
            if assay_type != 2:
                continue

            # Get active compounds
            active_url = (
                f"https://pubchem.ncbi.nlm.nih.gov/"
                f"rest/pug/assay/aid/{aid}/"
                f"concise/JSON"
                f"?activity=active"
            )
            cr = requests.get(active_url, timeout=30)
            if cr.status_code != 200:
                continue

            compounds = cr.json().get(
                "Table", {}
            ).get("Row", [])

            for comp in compounds[:100]:
                cells = comp.get("Cell", [])
                if len(cells) < 4:
                    continue
                try:
                    cid       = cells[0]
                    ac50_str  = cells[3]
                    ac50      = float(ac50_str)
                    pic50     = -np.log10(ac50 * 1e-6)

                    if not 3 <= pic50 <= 12:
                        continue

                    # Get SMILES from CID
                    smiles_url = (
                        f"https://pubchem.ncbi.nlm.nih.gov/"
                        f"rest/pug/compound/cid/{cid}/"
                        f"property/IsomericSMILES/JSON"
                    )
                    sr = requests.get(
                        smiles_url, timeout=10
                    )
                    if sr.status_code != 200:
                        continue
                    props = sr.json().get(
                        "PropertyTable", {}
                    ).get("Properties", [{}])[0]
                    smiles = props.get(
                        "IsomericSMILES", ""
                    )
                    if not smiles:
                        continue

                    pubchem_rows.append({
                        "smiles":        smiles,
                        "pIC50":         round(pic50, 3),
                        "IC50_nM":       10**(9-pic50),
                        "source":        "PubChem",
                        "activity_type": "AC50",
                        "chembl_id":     f"CID_{cid}",
                    })
                except:
                    continue

except Exception as e:
    print(f"  PubChem query failed: {e}")

pubchem_df = pd.DataFrame(pubchem_rows)
print(f"  Valid PubChem records: {len(pubchem_df)}")

# ── Combine all sources ───────────────────────────────────
print(f"\n{'=' * 60}")
print("COMBINING AND DEDUPLICATING")
print("=" * 60)

# Load existing ChEMBL IC50 data
existing_df = pd.read_csv(EXISTING_CSV)
existing_df["source"]        = "ChEMBL_IC50"
existing_df["activity_type"] = "IC50"
if "chembl_id" not in existing_df.columns:
    existing_df["chembl_id"] = ""

# Standardise columns
keep_cols = [
    "smiles", "pIC50", "IC50_nM",
    "source", "activity_type", "chembl_id"
]

# Use canonical_smiles if smiles not present
if "canonical_smiles" in existing_df.columns:
    existing_df = existing_df.rename(
        columns={"canonical_smiles": "smiles"}
    )

all_dfs = [existing_df[keep_cols]]
if len(chembl_df) > 0:
    all_dfs.append(chembl_df[keep_cols])
if len(bindingdb_df) > 0:
    all_dfs.append(bindingdb_df[keep_cols])
if len(pubchem_df) > 0:
    all_dfs.append(pubchem_df[keep_cols])

combined_df = pd.concat(all_dfs, ignore_index=True)
print(f"Total before deduplication: {len(combined_df)}")

# ── Standardise SMILES and deduplicate ────────────────────
print("Standardising SMILES...")
valid_rows = []
inchikeys  = set()

for _, row in combined_df.iterrows():
    smiles = row["smiles"]
    if not smiles or pd.isna(smiles):
        continue

    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        continue

    # Canonicalise
    canon_smiles = Chem.MolToSmiles(mol)

    # Deduplicate by InChIKey
    from rdkit.Chem.inchi import MolToInchiKey
    try:
        inchikey = MolToInchiKey(mol)
    except:
        inchikey = canon_smiles

    if inchikey in inchikeys:
        continue
    inchikeys.add(inchikey)

    # Basic filters
    mw = Descriptors.MolWt(mol)
    if mw < 100 or mw > 1000:
        continue

    valid_rows.append({
        "smiles":        canon_smiles,
        "pIC50":         row["pIC50"],
        "IC50_nM":       row["IC50_nM"],
        "source":        row["source"],
        "activity_type": row["activity_type"],
        "chembl_id":     row.get("chembl_id", ""),
        "MW":            round(mw, 2),
        "LogP":          round(Descriptors.MolLogP(mol), 2),
        "HBD":           rdMolDescriptors.CalcNumHBD(mol),
        "HBA":           rdMolDescriptors.CalcNumHBA(mol),
    })

final_df = pd.DataFrame(valid_rows)
final_df = final_df.sort_values("pIC50", ascending=False)
final_df.to_csv(EXPANDED_CSV, index=False)

# ── Summary ───────────────────────────────────────────────
print(f"\n{'=' * 60}")
print("EXPANSION SUMMARY")
print("=" * 60)
print(f"Original ChEMBL IC50:     {len(existing_df)}")
print(f"New ChEMBL (expanded):    {len(chembl_df)}")
print(f"BindingDB:                {len(bindingdb_df)}")
print(f"PubChem BioAssay:         {len(pubchem_df)}")
print(f"After deduplication:      {len(final_df)}")
print(f"\nSource breakdown:")
print(final_df["source"].value_counts().to_string())
print(f"\npIC50 distribution:")
print(f"  > 9  (< 1 nM):       {(final_df['pIC50']>9).sum()}")
print(f"  8-9  (1-10 nM):      {((final_df['pIC50']>=8)&(final_df['pIC50']<9)).sum()}")
print(f"  7-8  (10-100 nM):    {((final_df['pIC50']>=7)&(final_df['pIC50']<8)).sum()}")
print(f"  6-7  (0.1-1 uM):     {((final_df['pIC50']>=6)&(final_df['pIC50']<7)).sum()}")
print(f"  5-6  (1-10 uM):      {((final_df['pIC50']>=5)&(final_df['pIC50']<6)).sum()}")
print(f"  < 5  (> 10 uM):      {(final_df['pIC50']<5).sum()}")
print(f"\nSaved → {EXPANDED_CSV}")
