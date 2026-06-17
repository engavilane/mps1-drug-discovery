"""
NTD TPR domain receptor preparation for docking.

Structure: 4H7Y (1.8 A, apo TPR domain)
Binding site: Pocket 42 (fpocket) = Hec1/Ndc80 interface
Reference: Screpanti et al. 2011

Steps:
  1. Run fpocket on 4H7Y
  2. Identify pocket closest to Hec1 interface
  3. Define docking grid
  4. Prepare receptor PDBQT
"""

import subprocess
import numpy as np
import pandas as pd
from pathlib import Path
from Bio.PDB import PDBParser
import requests
import warnings
warnings.filterwarnings('ignore')

RAW_DIR      = Path("ntd_domain/data/raw")
RECEPTOR_DIR = Path("ntd_domain/data/receptor")
ANALYSIS_DIR = Path("ntd_domain/analysis/binding_site")

RECEPTOR_DIR.mkdir(parents=True, exist_ok=True)
ANALYSIS_DIR.mkdir(parents=True, exist_ok=True)

PDB_ID   = "4H7Y"
PDB_PATH = RAW_DIR / f"{PDB_ID}.pdb"

# Hec1 interface residues from Screpanti et al. 2011
HEC1_INTERFACE = {
    104, 107, 108, 111, 140, 141,
    144, 145, 173, 174, 177, 178
}

print("NTD TPR Domain Receptor Preparation")
print("=" * 55)

# Download 4H7Y if not present
if not PDB_PATH.exists():
    print(f"Downloading {PDB_ID}...")
    r = requests.get(
        f"https://files.rcsb.org/download/{PDB_ID}.pdb",
        timeout=30
    )
    PDB_PATH.write_text(r.text)
    print(f"  Downloaded: {PDB_PATH.stat().st_size} bytes")
else:
    print(f"  {PDB_ID}.pdb already present")

# Run fpocket
print(f"\nRunning fpocket on {PDB_ID}...")
fpocket_out = ANALYSIS_DIR / f"{PDB_ID}_out"

result = subprocess.run(
    ["fpocket", "-f", str(PDB_PATH)],
    capture_output=True, text=True,
    cwd=str(ANALYSIS_DIR)
)

# Move output to correct location
raw_out = RAW_DIR / f"{PDB_ID}_out"
if raw_out.exists() and not fpocket_out.exists():
    import shutil
    shutil.move(str(raw_out), str(fpocket_out))

print(f"  fpocket output: {fpocket_out}")

# Parse fpocket pockets
print("\nParsing fpocket pockets...")
pqr_file = fpocket_out / f"{PDB_ID}_pockets.pqr"

pockets = {}
with open(pqr_file) as f:
    for line in f:
        if not line.startswith('ATOM'):
            continue
        parts = line.split()
        try:
            pocket_num = int(parts[4])
            x = float(parts[5])
            y = float(parts[6])
            z = float(parts[7])
            if pocket_num not in pockets:
                pockets[pocket_num] = []
            pockets[pocket_num].append([x, y, z])
        except:
            continue

print(f"  {len(pockets)} pockets detected")

# Get Hec1 interface centre from crystal structure
parser = PDBParser(QUIET=True)
struct = parser.get_structure(PDB_ID, str(PDB_PATH))

interface_coords = []
for model in struct:
    for chain in model:
        for residue in chain:
            resid = residue.get_id()[1]
            if resid in HEC1_INTERFACE:
                if 'CA' in residue:
                    interface_coords.append(
                        residue['CA'].get_vector().get_array()
                    )

interface_centre = np.mean(interface_coords, axis=0)
print(f"\nHec1 interface centre: {interface_centre.round(2)}")
print(f"Interface residues found: {len(interface_coords)}")

# Find closest pocket to interface
print("\nPocket distances to Hec1 interface:")
pocket_distances = []
for pocket_num, coords in pockets.items():
    centre = np.mean(coords, axis=0)
    dist   = np.linalg.norm(centre - interface_centre)
    pocket_distances.append((dist, pocket_num, centre))

pocket_distances.sort()
for dist, num, centre in pocket_distances[:5]:
    print(f"  Pocket {num}: {dist:.2f} A from interface "
          f"[{centre[0]:.1f}, {centre[1]:.1f}, {centre[2]:.1f}]")

# Select best pocket
best_dist, best_pocket, pocket_centre = pocket_distances[0]
print(f"\nSelected: Pocket {best_pocket} "
      f"({best_dist:.2f} A from Hec1 interface)")

# Define grid centre as midpoint
grid_centre = (pocket_centre + interface_centre) / 2
print(f"Grid centre: {grid_centre.round(2)}")

# Save pocket analysis
analysis_df = pd.DataFrame([
    {
        'pocket_id':       num,
        'dist_to_hec1_A':  round(dist, 2),
        'centre_x':        round(centre[0], 2),
        'centre_y':        round(centre[1], 2),
        'centre_z':        round(centre[2], 2),
        'selected':        (num == best_pocket),
    }
    for dist, num, centre in pocket_distances[:10]
])
analysis_df.to_csv(
    ANALYSIS_DIR / "pocket_analysis.csv", index=False
)

# Save grid config
grid_config = f"""# NTD TPR domain docking grid — {PDB_ID}
# Binding site: Pocket {best_pocket} (fpocket)
#               = Hec1/Ndc80 interface
# Reference: Screpanti et al. 2011
# Distance pocket to interface: {best_dist:.2f} A
# Grid size: 30x30x30 A (PPI interface, conservative)
# AlphaFold2-Multimer: TPR pLDDT=96.1 (confident)
#                      Hec1 pLDDT=33.3 (disordered — expected)

center_x = {grid_centre[0]:.2f}
center_y = {grid_centre[1]:.2f}
center_z = {grid_centre[2]:.2f}

size_x = 30
size_y = 30
size_z = 30

exhaustiveness = 16
num_modes      = 9
energy_range   = 3
"""

grid_path = RECEPTOR_DIR / f"{PDB_ID}_grid.txt"
grid_path.write_text(grid_config)
print(f"Grid config saved → {grid_path}")

# Prepare PDBQT
print(f"\nPreparing {PDB_ID} receptor PDBQT...")
pdbqt_path = RECEPTOR_DIR / f"{PDB_ID}_receptor"

result = subprocess.run([
    "mk_prepare_receptor.py",
    "-i", str(PDB_PATH),
    "-o", str(pdbqt_path),
    "--allow_bad_res",
    "--default_altloc", "A",
    "-p"
], capture_output=True, text=True)

final_pdbqt = Path(f"{pdbqt_path}.pdbqt")
if final_pdbqt.exists():
    print(f"  ✓ {final_pdbqt.name} "
          f"({final_pdbqt.stat().st_size} bytes)")
else:
    print(f"  ✗ PDBQT preparation failed")
    print(result.stderr[-200:])

print(f"\n{'=' * 55}")
print(f"NTD RECEPTOR PREPARATION COMPLETE")
print(f"{'=' * 55}")
print(f"  Structure:    {PDB_ID} (1.8 A, apo TPR domain)")
print(f"  Best pocket:  Pocket {best_pocket} "
      f"({best_dist:.2f} A from Hec1 interface)")
print(f"  Grid centre:  {grid_centre.round(2)}")
print(f"  Grid size:    30 x 30 x 30 A")
print(f"  PDBQT:        {final_pdbqt.name}")
