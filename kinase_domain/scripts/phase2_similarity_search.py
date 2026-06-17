import pubchempy as pcp
import pandas as pd
import numpy as np
import requests
import os
import time
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors, DataStructs


# Paths
LIGANDS_DIR    = Path("kinase_domain/data/ligands/raw")
PHASE2_DIR     = Path("kinase_domain/data/phase2b")
PHASE2_SDF     = PHASE2_DIR / "raw"
PHASE2_PDBQT   = PHASE2_DIR / "pdbqt"
OUTPUT_DIR     = Path("kinase_domain/analysis/phase2b")

for d in [PHASE2_SDF, PHASE2_PDBQT, OUTPUT_DIR]:
    d.mkdir(parents=True, exist_ok=True)


# Pugh et al. 2022 filters (Section 4.1) 
FILTERS = {
    "MW":        (270, 500),
    "HBD":       (1,   7),
    "HBA":       (3,   11),
    "RotBonds":  (1,   9),
    "HeavyAtoms":(19,  37),
    "LogP":      (1,   7),
}
TANIMOTO_THRESHOLD = 0.7


# Load seeds from CSV

seeds_df = pd.read_csv(
    "kinase_domain/analysis/phase2/phase2b_seeds.csv"
)
SEEDS = seeds_df['smiles'].tolist()


# Load existing ligand fingerprints to avoid duplicates 
print("Loading existing ligand fingerprints...")
existing_fps  = []
existing_smiles = set()

for sdf_file in LIGANDS_DIR.glob("*.sdf"):
    mol = Chem.MolFromMolFile(str(sdf_file))
    if mol:
        fp = AllChem.GetMorganFingerprintAsBitVect(
            mol, radius=2, nBits=1024
        )
        existing_fps.append(fp)
        smi = Chem.MolToSmiles(mol)
        existing_smiles.add(smi)

print(f"  {len(existing_fps)} existing compounds loaded\n")


# Helper: check filters 
def passes_filters(mol):
    """Check Pugh et al. physicochemical filters."""
    mw    = Descriptors.MolWt(mol)
    hbd   = rdMolDescriptors.CalcNumHBD(mol)
    hba   = rdMolDescriptors.CalcNumHBA(mol)
    rotb  = rdMolDescriptors.CalcNumRotatableBonds(mol)
    heavy = mol.GetNumHeavyAtoms()
    logp  = Descriptors.MolLogP(mol)

    return (
        FILTERS["MW"][0]        <= mw    <= FILTERS["MW"][1]        and
        FILTERS["HBD"][0]       <= hbd   <= FILTERS["HBD"][1]       and
        FILTERS["HBA"][0]       <= hba   <= FILTERS["HBA"][1]       and
        FILTERS["RotBonds"][0]  <= rotb  <= FILTERS["RotBonds"][1]  and
        FILTERS["HeavyAtoms"][0]<= heavy <= FILTERS["HeavyAtoms"][1]and
        FILTERS["LogP"][0]      <= logp  <= FILTERS["LogP"][1]
    )

def is_novel(mol, existing_fps, threshold=0.9):
    """Check compound is not too similar to existing ones."""
    fp = AllChem.GetMorganFingerprintAsBitVect(
        mol, radius=2, nBits=1024
    )
    if existing_fps:
        sims = DataStructs.BulkTanimotoSimilarity(fp, existing_fps)
        if max(sims) >= threshold:
            return False, max(sims)
    return True, 0.0


# Main search loop 
all_candidates = []
seen_cids      = set()

print("=" * 60)
print("PHASE 2 — PUBCHEM SIMILARITY SEARCH")
print(f"Tanimoto threshold: {TANIMOTO_THRESHOLD}")
print(f"Seeds: {', '.join(SEEDS)}")
print("=" * 60)

for seed_idx, seed_smiles in enumerate(SEEDS):
    seed_name = f"seed_{seed_idx}"
    print(f"\nSearching similar compounds for: {seed_name} "
          f"(pIC50={seeds_df.iloc[seed_idx]['pIC50']:.3f})")

    try:
        seed_results = pcp.get_compounds(
            seed_smiles, "smiles"
        )
        if not seed_results:
            print(f"   Could not find CID for {seed_name}")
            continue

        seed_cid = seed_results[0].cid
        print(f"  Seed CID: {seed_cid}")

        # PubChem similarity search
        similar = pcp.get_compounds(
            seed_cid,
            searchtype="similarity",
            Threshold=int(TANIMOTO_THRESHOLD * 100),
            MaxRecords=200
        )
        print(f"  Found {len(similar)} similar compounds")

    except Exception as e:
        print(f"   PubChem search error: {e}")
        continue


    # Filter candidates 
    accepted  = 0
    rejected  = 0

    for compound in similar:
        cid = compound.cid

        # Skip if already seen
        if cid in seen_cids:
            continue
        seen_cids.add(cid)

        # Skip seed itself
        if cid == seed_cid:
            continue

        # Get SMILES
        smiles = compound.connectivity_smiles
        if not smiles:
            continue

        # Parse molecule
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue

        # Skip if same as existing ligand
        if smiles in existing_smiles:
            rejected += 1
            continue

        # Apply physicochemical filters
        if not passes_filters(mol):
            rejected += 1
            continue

        # Check novelty vs existing dataset
        novel, max_sim = is_novel(mol, existing_fps)
        if not novel:
            rejected += 1
            continue


        # Download 3D SDF 
        try:
            url = (f"https://pubchem.ncbi.nlm.nih.gov/rest/pug"
                   f"/compound/cid/{cid}/SDF?record_type=3d")
            r = requests.get(url, timeout=10)

            if r.status_code != 200:
                # Try 2D + note it
                url2d = (f"https://pubchem.ncbi.nlm.nih.gov/rest"
                         f"/pug/compound/cid/{cid}/SDF"
                         f"?record_type=2d")
                r = requests.get(url2d, timeout=10)
                if r.status_code != 200:
                    continue
                conformer = "2d"
            else:
                conformer = "3d"

            filename = f"phase2b_{seed_name}_{cid}.sdf"
            filepath = PHASE2_SDF / filename
            with open(filepath, "w") as f:
                f.write(r.text)

            # Compute descriptors
            mw    = Descriptors.MolWt(mol)
            logp  = Descriptors.MolLogP(mol)
            hbd   = rdMolDescriptors.CalcNumHBD(mol)
            hba   = rdMolDescriptors.CalcNumHBA(mol)
            tpsa  = Descriptors.TPSA(mol)
            rotb  = rdMolDescriptors.CalcNumRotatableBonds(mol)
            heavy = mol.GetNumHeavyAtoms()
            rings = rdMolDescriptors.CalcNumRings(mol)
            arom  = rdMolDescriptors.CalcNumAromaticRings(mol)

            all_candidates.append({
                "cid":          cid,
                "seed":         seed_name,
                "smiles":       smiles,
                "sdf_file":     str(filepath),
                "conformer":    conformer,
                "MW":           round(mw, 2),
                "LogP":         round(logp, 2),
                "HBD":          hbd,
                "HBA":          hba,
                "TPSA":         round(tpsa, 2),
                "RotBonds":     rotb,
                "HeavyAtoms":   heavy,
                "Rings":        rings,
                "AromaticRings":arom,
            })
            accepted += 1

        except Exception as e:
            continue

        time.sleep(0.3)

    print(f"  ✓ Accepted: {accepted} | Rejected: {rejected}")


# Save candidate list 
candidates_df = pd.DataFrame(all_candidates)

if candidates_df.empty:
    print("\n No candidates found — try lowering threshold")
else:
    # Remove duplicates across seeds
    candidates_df = candidates_df.drop_duplicates(
        subset="cid"
    ).reset_index(drop=True)

    candidates_df.to_csv(
        OUTPUT_DIR / "phase2_candidates.csv", index=False
    )

    print(f"\n{'=' * 60}")
    print(f"PHASE 2 SIMILARITY SEARCH COMPLETE")
    print(f"{'=' * 60}")
    print(f"Total novel candidates: {len(candidates_df)}")
    print(f"SDF files saved to:     {PHASE2_SDF}")
    print(f"Candidate list:         {OUTPUT_DIR}/phase2_candidates.csv")
    print(f"\nDescriptor summary:")
    print(f"  MW:    {candidates_df['MW'].mean():.1f} ± "
          f"{candidates_df['MW'].std():.1f}")
    print(f"  LogP:  {candidates_df['LogP'].mean():.2f} ± "
          f"{candidates_df['LogP'].std():.2f}")
    print(f"  TPSA:  {candidates_df['TPSA'].mean():.1f} ± "
          f"{candidates_df['TPSA'].std():.1f}")
    print(f"\nCandidates per seed:")
    print(candidates_df.groupby("seed").size().to_string())
