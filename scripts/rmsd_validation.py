import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign
from pathlib import Path
from itertools import permutations
from scipy.spatial.distance import cdist
from scipy.optimize import linear_sum_assignment


# Paths 
CRYSTAL_SDF  = "data/ligands/raw/5LJJ_AD5.sdf"
REDOCK_PDBQT = "docking/results/5LJJ_AD5_redock_out.pdbqt"
OUTPUT_DIR   = Path("analysis/validation")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


import numpy as np

# Load crystal coordinates directly from PDB
print("Loading crystal coordinates from 5LJJ.pdb...")
crystal_coords = np.load(
    "data/receptor/reversine_crystal_coords.npy"
)
print(f"  Crystal pose: {len(crystal_coords)} heavy atoms")


# Parse docked poses from PDBQT 
print("\nParsing docked poses...")

def parse_pdbqt_poses(filepath):
    """Extract all poses from a multi-model PDBQT file."""
    poses = []
    current = []
    with open(filepath) as f:
        for line in f:
            if line.startswith("MODEL"):
                current = []
            elif line.startswith("ENDMDL"):
                poses.append(current)
            elif line.startswith(("ATOM", "HETATM")):
                current.append(line)
    return poses



# Convert PDBQT pose to RDKit mol 
def pdbqt_to_coords(pdbqt_lines):
    """Extract heavy atom coordinates from PDBQT lines."""
    coords = []
    for line in pdbqt_lines:
        atom_type = line[77:79].strip() if len(line) > 77 else ""
        # Skip hydrogens
        if atom_type.startswith("H"):
            continue
        try:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            coords.append([x, y, z])
        except ValueError:
            continue
    return np.array(coords)

poses = parse_pdbqt_poses(REDOCK_PDBQT)
print(f"  Found {len(poses)} docked poses")


# Calculate RMSD for each pose 
print("\nCalculating RMSD for each pose vs crystal structure...")
print("=" * 50)


results = []

for i, pose_lines in enumerate(poses):
    docked_coords = pdbqt_to_coords(pose_lines)

    if len(docked_coords) == 0:
        continue

    n = min(len(crystal_coords), len(docked_coords))
    if n == 0:
        continue

    # Distance matrix between all pairs of atoms
    dist_matrix = cdist(
        crystal_coords[:n], docked_coords[:n]
    )

    # Optimal assignment
    row_ind, col_ind = linear_sum_assignment(dist_matrix)
    matched_crystal = crystal_coords[row_ind]
    matched_docked  = docked_coords[col_ind]

    diff = matched_crystal - matched_docked
    rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))

    results.append({"pose": i+1, "rmsd": round(rmsd, 3)})
    print(f"  Pose {i+1}: RMSD = {rmsd:.3f} Å"
          f"  {'✓ PASS' if rmsd < 2.0 else '✗ FAIL'}")


# Summary 
print("\n" + "=" * 50)
best = min(results, key=lambda x: x["rmsd"])
print(f"Best RMSD: {best['rmsd']:.3f} Å (Pose {best['pose']})")

if best["rmsd"] < 2.0:
    print("✓ VALIDATION PASSED — docking correctly")
elif best["rmsd"] < 3.0:
    print(" BORDERLINE — pose is close but above 2 Å threshold")
else:
    print(" VALIDATION FAILED — docking does not")


# Save results 
df = pd.DataFrame(results)
df.to_csv(OUTPUT_DIR / "rmsd_validation.csv", index=False)
print(f"\nResults saved → {OUTPUT_DIR}/rmsd_validation.csv")


