"""
Ensemble receptor preparation for Mps1/TTK.

Uses all 75 KLIFS-curated structures (all DFG-in confirmed).
Pre-filters to top 30 by quality/resolution before expensive
alignment and RMSD calculations.

Method:
  1. Get KLIFS quality scores (API)
  2. Pre-filter to top 30 (quality >= 8.0, resolution <= 3.0 A)
  3. Align filtered structures to 5LJJ via DSSP rigid core
  4. Compute pairwise binding site Ca RMSD matrix
  5. Greedy maximum diversity selection (quality-weighted)
  6. Prepare 7 receptor PDBQTs
"""

from Bio.PDB import PDBParser, PDBIO, Select, Superimposer
from Bio.PDB.DSSP import DSSP
import numpy as np
import pandas as pd
from pathlib import Path
import subprocess
import requests
import warnings
warnings.filterwarnings('ignore')

# ── Configuration ─────────────────────────────────────────
RAW_DIR      = Path("kinase_domain/data/raw")
RECEPTOR_DIR = Path("kinase_domain/data/receptor/ensemble")
RECEPTOR_DIR.mkdir(parents=True, exist_ok=True)

REF_PDB    = "kinase_domain/data/raw/5LJJ.pdb"
N_ENSEMBLE = 7

# Quality pre-filter thresholds
MIN_QUALITY    = 8.0
MAX_RESOLUTION = 3.0
TOP_N_PREFILTER = 30

STANDARD_AA = [
    'ALA','ARG','ASN','ASP','CYS','GLN','GLU',
    'GLY','HIS','ILE','LEU','LYS','MET','PHE',
    'PRO','SER','THR','TRP','TYR','VAL'
]
EXCLUDE_RES = set(range(672, 691))  # activation + P+1 loop
EXCLUDE_RES.update(range(603, 606)) # hinge loop
RIGID_SS    = {'H', 'G', 'I', 'E'}

BINDING_SITE_RES = set(range(520, 700))

parser = PDBParser(QUIET=True)

# ── Step 0 — Get KLIFS metadata and pre-filter ────────────
print("=" * 60)
print("STEP 0 — KLIFS quality pre-filter")
print("=" * 60)

print("Fetching KLIFS metadata for all TTK structures...")
r = requests.get(
    'https://klifs.net/api/structures_list',
    params={'kinase_ID': 326},
    timeout=30
)
klifs_df = pd.DataFrame(r.json())
klifs_df['pdb'] = klifs_df['pdb'].str.upper()

# Keep best entry per PDB (highest quality score)
klifs_best = klifs_df.sort_values(
    'quality_score', ascending=False
).drop_duplicates('pdb').reset_index(drop=True)

print(f"  Total KLIFS structures: {len(klifs_best)}")
print(f"  All DFG-in: {(klifs_df['DFG'] == 'in').all()}")
print(f"  All aC-in:  {(klifs_df['aC_helix'] == 'in').all()}")

# Pre-filter by quality and resolution
pre_filtered = klifs_best[
    (klifs_best['resolution'].astype(float)
     <= MAX_RESOLUTION) &
    (klifs_best['quality_score'].astype(float)
     >= MIN_QUALITY)
].head(TOP_N_PREFILTER).copy()

print(f"\nPre-filter criteria:")
print(f"  Quality score >= {MIN_QUALITY}")
print(f"  Resolution    <= {MAX_RESOLUTION} A")
print(f"  Top N:           {TOP_N_PREFILTER}")
print(f"  Passing:         {len(pre_filtered)} structures")

# Always include 5LJJ as reference
if '5LJJ' not in pre_filtered['pdb'].values:
    ljj_row      = klifs_best[klifs_best['pdb'] == '5LJJ']
    pre_filtered = pd.concat(
        [ljj_row, pre_filtered]
    ).drop_duplicates('pdb').reset_index(drop=True)
    print(f"  Added 5LJJ (reference) -> {len(pre_filtered)} total")

candidate_pdbs = set(pre_filtered['pdb'].values)

print(f"\nCandidates:")
for _, row in pre_filtered.sort_values(
    'quality_score', ascending=False
).iterrows():
    print(f"  {row['pdb']:<6} quality={row['quality_score']} "
          f"res={row['resolution']} A "
          f"ligand={row['ligand']}")

# ── Step 1 — Find available PDB files ────────────────────
print(f"\n{'=' * 60}")
print("STEP 1 — Loading candidate structures")
print("=" * 60)

all_pdbs = sorted(RAW_DIR.glob("*.pdb"))

exclude_names = {
    'receptor_clean', 'receptor',
    '5N7V_aligned', '4JS8_aligned',
    '5NAD_aligned', '7LQD_aligned',
    '7LQD_aligned_trimmed'
}

receptor_pdbs = [
    p for p in all_pdbs
    if p.stem.upper() in candidate_pdbs
    and not any(ex in p.stem for ex in exclude_names)
    and len(p.stem) == 4
]
print(f"  {len(receptor_pdbs)} candidate PDB files found")

# ── Helper functions ──────────────────────────────────────
def get_rigid_ca(structure, pdb_file):
    """Get Ca atoms in rigid SS elements via DSSP."""
    model = list(structure.get_models())[0]
    try:
        dssp      = DSSP(model, str(pdb_file), dssp='mkdssp')
        dssp_dict = {
            (k[0], k[1][1]): v[2]
            for k, v in dssp.property_dict.items()
        }
    except Exception:
        dssp_dict = {}

    residues = {}
    for mod in structure:
        for chain in mod:
            for residue in chain:
                resid    = residue.get_id()[1]
                resname  = residue.get_resname()
                chain_id = chain.get_id()
                if resname not in STANDARD_AA:
                    continue
                if resid in EXCLUDE_RES:
                    continue
                if 'CA' not in residue:
                    continue
                ss = dssp_dict.get((chain_id, resid), 'C')
                if ss not in RIGID_SS:
                    continue
                residues[resid] = residue['CA']
    return residues

def get_binding_site_coords(structure):
    """Get Ca coordinates for binding site residues."""
    coords = []
    resids = []
    for model in structure:
        for chain in model:
            for residue in chain:
                resid   = residue.get_id()[1]
                resname = residue.get_resname()
                if resname not in STANDARD_AA:
                    continue
                if resid not in BINDING_SITE_RES:
                    continue
                if 'CA' not in residue:
                    continue
                coords.append(
                    list(residue['CA'].get_vector())
                )
                resids.append(resid)
    return np.array(coords), resids

class ProteinSelect(Select):
    def accept_residue(self, residue):
        return residue.get_resname() in STANDARD_AA

# ── Step 2 — Align structures to 5LJJ ────────────────────
print(f"\n{'=' * 60}")
print("STEP 2 — Aligning structures to 5LJJ")
print("=" * 60)

ref_structure = parser.get_structure('5LJJ', REF_PDB)
ref_rigid     = get_rigid_ca(ref_structure, REF_PDB)
print(f"  Reference rigid Ca: {len(ref_rigid)}")

aligned_structures = {}
sup_rmsds          = {}
quality_scores     = {}

for pdb_path in receptor_pdbs:
    pdb_id = pdb_path.stem.upper()
    try:
        mob_structure = parser.get_structure(
            pdb_id, str(pdb_path)
        )
        mob_rigid = get_rigid_ca(mob_structure, pdb_path)

        common = sorted(
            set(ref_rigid.keys()) & set(mob_rigid.keys())
        )
        if len(common) < 10:
            print(f"  ✗ {pdb_id}: too few common Ca "
                  f"({len(common)})")
            continue

        ref_ca = [ref_rigid[r] for r in common]
        mob_ca = [mob_rigid[r] for r in common]

        sup = Superimposer()
        sup.set_atoms(ref_ca, mob_ca)
        sup.apply(mob_structure.get_atoms())

        aligned_structures[pdb_id] = mob_structure
        sup_rmsds[pdb_id]          = round(sup.rms, 3)

        # KLIFS quality score
        row = pre_filtered[pre_filtered['pdb'] == pdb_id]
        quality_scores[pdb_id] = float(
            row['quality_score'].values[0]
        ) if len(row) > 0 else 0.0

        print(f"  ✓ {pdb_id}: sup={sup.rms:.3f} A "
              f"quality={quality_scores[pdb_id]} "
              f"({len(common)} Ca)")

    except Exception as e:
        print(f"  ✗ {pdb_id}: {e}")

print(f"\n  Successfully aligned: {len(aligned_structures)}")

# ── Step 3 — Compute pairwise RMSD matrix ─────────────────
print(f"\n{'=' * 60}")
print("STEP 3 — Pairwise binding site RMSD matrix")
print("=" * 60)

pdb_ids = list(aligned_structures.keys())
n       = len(pdb_ids)
print(f"  Computing {n}x{n} = {n*n} pairwise RMSDs...")

all_coords = {}
all_resids = {}
for pdb_id in pdb_ids:
    coords, resids = get_binding_site_coords(
        aligned_structures[pdb_id]
    )
    all_coords[pdb_id] = coords
    all_resids[pdb_id] = resids

rmsd_matrix = np.zeros((n, n))
for i, id1 in enumerate(pdb_ids):
    for j, id2 in enumerate(pdb_ids):
        if i >= j:
            continue
        res1       = set(all_resids[id1])
        res2       = set(all_resids[id2])
        common_res = sorted(res1 & res2)
        if len(common_res) < 10:
            rmsd_matrix[i, j] = 999
            rmsd_matrix[j, i] = 999
            continue
        idx1 = [all_resids[id1].index(r)
                for r in common_res]
        idx2 = [all_resids[id2].index(r)
                for r in common_res]
        c1   = all_coords[id1][idx1]
        c2   = all_coords[id2][idx2]
        n_c  = min(len(c1), len(c2))
        diff = c1[:n_c] - c2[:n_c]
        rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
        rmsd_matrix[i, j] = round(rmsd, 3)
        rmsd_matrix[j, i] = round(rmsd, 3)

rmsd_df = pd.DataFrame(
    rmsd_matrix, index=pdb_ids, columns=pdb_ids
)
rmsd_df.to_csv(RECEPTOR_DIR / "pairwise_rmsd_matrix.csv")

valid = rmsd_matrix[(rmsd_matrix > 0) & (rmsd_matrix < 999)]
print(f"  Mean pairwise RMSD: {valid.mean():.3f} A")
print(f"  Max pairwise RMSD:  {valid.max():.3f} A")
print(f"  Min pairwise RMSD:  {valid.min():.3f} A")

# ── Step 4 — Quality-weighted diversity selection ─────────
print(f"\n{'=' * 60}")
print(f"STEP 4 — Selecting {N_ENSEMBLE} diverse structures")
print("=" * 60)
print("Method: greedy maximum diversity, quality-weighted")
print("        score = min_dist * (1 + 0.1 * quality)")
print()

selected  = ['5LJJ']
remaining = [p for p in pdb_ids if p != '5LJJ']

while len(selected) < N_ENSEMBLE and remaining:
    best_pdb   = None
    best_score = -1

    for cand in remaining:
        cand_idx = pdb_ids.index(cand)
        dists    = [
            rmsd_matrix[cand_idx, pdb_ids.index(sel)]
            for sel in selected
            if rmsd_matrix[
                cand_idx, pdb_ids.index(sel)
            ] < 999
        ]
        if not dists:
            continue
        min_dist = min(dists)
        quality  = quality_scores.get(cand, 0)
        score    = min_dist * (1 + 0.1 * quality)

        if score > best_score:
            best_score = score
            best_pdb   = cand

    if best_pdb:
        selected.append(best_pdb)
        remaining.remove(best_pdb)
        q   = quality_scores.get(best_pdb, 0)
        row = pre_filtered[
            pre_filtered['pdb'] == best_pdb
        ]
        res = row['resolution'].values[0] \
            if len(row) > 0 else 'N/A'
        lig = row['ligand'].values[0] \
            if len(row) > 0 else 'N/A'
        print(f"  Round {len(selected)-1}: {best_pdb} "
              f"selected (score={best_score:.3f}, "
              f"quality={q}, res={res}A, ligand={lig})")

# ── Step 5 — Final ensemble summary ──────────────────────
print(f"\n{'=' * 60}")
print(f"FINAL ENSEMBLE ({len(selected)} structures)")
print("=" * 60)
print(f"{'PDB':<6} {'Sup RMSD':>9} {'Quality':>8} "
      f"{'Res':>6} {'Ligand':>8} {'Mean dist':>10}")
print("-" * 55)

for pdb_id in selected:
    row  = pre_filtered[pre_filtered['pdb'] == pdb_id]
    res  = row['resolution'].values[0] \
        if len(row) > 0 else 'N/A'
    lig  = row['ligand'].values[0] \
        if len(row) > 0 else 'N/A'
    qual = quality_scores.get(pdb_id, 0)
    sup  = sup_rmsds.get(pdb_id, 0)
    idx  = pdb_ids.index(pdb_id)
    dists = [
        rmsd_matrix[idx, pdb_ids.index(s)]
        for s in selected if s != pdb_id
        and rmsd_matrix[idx, pdb_ids.index(s)] < 999
    ]
    mean_d = np.mean(dists) if dists else 0
    print(f"{pdb_id:<6} {sup:>9.3f} A {qual:>8.1f} "
          f"{str(res):>6} A {str(lig):>8} "
          f"{mean_d:>10.3f} A")

# ── Step 6 — Save aligned receptors + prepare PDBQTs ─────
print(f"\n{'=' * 60}")
print("STEP 6 — Preparing receptor PDBQTs")
print("=" * 60)

# Clean old ensemble files
for old in RECEPTOR_DIR.glob("*_ensemble*"):
    old.unlink()

io            = PDBIO()
ensemble_info = []

for pdb_id in selected:
    structure = (
        ref_structure
        if pdb_id == '5LJJ'
        else aligned_structures[pdb_id]
    )

    aligned_pdb   = RECEPTOR_DIR / f"{pdb_id}_ensemble.pdb"
    aligned_pdbqt = RECEPTOR_DIR / \
        f"{pdb_id}_ensemble_receptor"

    io.set_structure(structure)
    io.save(str(aligned_pdb), ProteinSelect())

    result = subprocess.run([
        'mk_prepare_receptor.py',
        '-i', str(aligned_pdb),
        '-o', str(aligned_pdbqt),
        '--allow_bad_res', '-p'
    ], capture_output=True, text=True)

    pdbqt_path = Path(f"{aligned_pdbqt}.pdbqt")
    status     = "✓" if pdbqt_path.exists() else "✗"
    print(f"  {status} {pdb_id} → "
          f"{pdbqt_path.name}")

    row = pre_filtered[pre_filtered['pdb'] == pdb_id]
    ensemble_info.append({
        "pdb_id":        pdb_id,
        "sup_rmsd":      sup_rmsds.get(pdb_id, 0),
        "quality_score": quality_scores.get(pdb_id, 0),
        "resolution":    row['resolution'].values[0]
                         if len(row) > 0 else None,
        "ligand":        row['ligand'].values[0]
                         if len(row) > 0 else None,
        "dfg":           'in',
        "grich_distance": row['Grich_distance'].values[0]
                          if len(row) > 0 else None,
        "aligned_pdb":   str(aligned_pdb),
        "pdbqt":         str(pdbqt_path),
        "prepared":      pdbqt_path.exists(),
    })

ensemble_df = pd.DataFrame(ensemble_info)
ensemble_df.to_csv(
    RECEPTOR_DIR / "ensemble_info.csv", index=False
)

print(f"\n{'=' * 60}")
print("ENSEMBLE RECEPTOR PREPARATION COMPLETE")
print("=" * 60)
print(f"  Structures in KLIFS:    75 (all DFG-in)")
print(f"  After quality filter:   {len(pre_filtered)}")
print(f"  Successfully aligned:   {len(aligned_structures)}")
print(f"  Final ensemble:         {len(selected)}")
print(f"  Prepared PDBQTs:        "
      f"{ensemble_df['prepared'].sum()}/{len(selected)}")
print(f"\n  Ensemble: {selected}")
print(f"\nResults → {RECEPTOR_DIR}/")
