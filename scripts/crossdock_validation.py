"""
Cross-docking validation of reversine across multiple
Mps1 crystal structures representing different scaffolds.

Structures:
  5LJJ — Reversine (reference, self-docking)
  5N7V — MPI-0479605 (purine scaffold)        [done]
  7LQD — RMS-07 (pyrrolopyrimidine scaffold)
  5NAD — BAY-1217389 (methylbenzamide scaffold)
  4JS8 — Compound 401348 (indazole scaffold)
"""

from Bio.PDB import PDBParser, PDBIO, Select, Superimposer
from Bio.PDB.DSSP import DSSP
from scipy.spatial.distance import cdist
from scipy.optimize import linear_sum_assignment
from vina import Vina
import numpy as np
import subprocess
import pandas as pd
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# ── Configuration ─────────────────────────────────────────
STANDARD_AA = [
    'ALA','ARG','ASN','ASP','CYS','GLN','GLU',
    'GLY','HIS','ILE','LEU','LYS','MET','PHE',
    'PRO','SER','THR','TRP','TYR','VAL'
]
EXCLUDE_RESIDUES = set(range(672, 691))
EXCLUDE_RESIDUES.update(range(603, 606))
RIGID_SS = {'H', 'G', 'I', 'E'}

OUTPUT_DIR = Path("analysis/validation")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Target structures ─────────────────────────────────────
TARGETS = [
    {
        "pdb_id":  "5N7V",
        "ligand":  "MPI-0479605",
        "scaffold":"Purine",
        "done":    True,   # already done
        "best_rmsd": 1.332
    },
    {
        "pdb_id":  "7LQD",
        "ligand":  "RMS-07",
        "scaffold":"Pyrrolopyrimidine",
        "done":    False
    },
    {
        "pdb_id":  "5NAD",
        "ligand":  "BAY-1217389",
        "scaffold":"Methylbenzamide",
        "done":    False
    },
    {
        "pdb_id":  "4JS8",
        "ligand":  "Compound 401348",
        "scaffold":"Indazole",
        "done":    False
    },
]

# ── Helper functions ──────────────────────────────────────
parser = PDBParser(QUIET=True)

def get_rigid_ca(structure, pdb_file):
    """Get Cα atoms in rigid SS elements via DSSP."""
    model = list(structure.get_models())[0]
    try:
        dssp = DSSP(model, pdb_file, dssp='mkdssp')
        dssp_dict = {
            (k[0], k[1][1]): v[2]
            for k, v in dssp.property_dict.items()
        }
    except Exception as e:
        print(f"    DSSP warning: {e}")
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
                if resid in EXCLUDE_RESIDUES:
                    continue
                if 'CA' not in residue:
                    continue
                ss = dssp_dict.get((chain_id, resid), 'C')
                if ss not in RIGID_SS:
                    continue
                residues[resid] = residue['CA']
    return residues

def pdbqt_to_coords(lines):
    """Extract heavy atom coordinates from PDBQT lines."""
    coords = []
    for line in lines:
        atom_type = line[77:79].strip() if len(line) > 77 else ''
        if atom_type.startswith('H'):
            continue
        try:
            coords.append([
                float(line[30:38]),
                float(line[38:46]),
                float(line[46:54])
            ])
        except ValueError:
            continue
    return np.array(coords)

def parse_pdbqt_poses(filepath):
    """Parse all poses from multi-model PDBQT."""
    poses, current = [], []
    with open(filepath) as f:
        for line in f:
            if line.startswith('MODEL'):
                current = []
            elif line.startswith('ENDMDL'):
                poses.append(current)
            elif line.startswith(('ATOM', 'HETATM')):
                current.append(line)
    return poses

def calc_best_rmsd(poses, crystal_coords):
    """Calculate best RMSD across all poses."""
    results = []
    for pose_lines in poses:
        docked = pdbqt_to_coords(pose_lines)
        if len(docked) == 0:
            continue
        n = min(len(crystal_coords), len(docked))
        dist_matrix      = cdist(
            crystal_coords[:n], docked[:n]
        )
        row_ind, col_ind = linear_sum_assignment(dist_matrix)
        diff = crystal_coords[row_ind] - docked[col_ind]
        rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
        results.append(round(rmsd, 3))
    return min(results), results

# ── Load reference ────────────────────────────────────────
print("Loading reference structure (5LJJ)...")
ref_structure  = parser.get_structure(
    '5LJJ', 'data/raw/5LJJ.pdb'
)
ref_rigid      = get_rigid_ca(ref_structure, 'data/raw/5LJJ.pdb')
crystal_coords = np.load(
    'data/receptor/reversine_crystal_coords.npy'
)
print(f"  Reference rigid Cα: {len(ref_rigid)}")

# ── Run cross-docking for each target ────────────────────
all_results = []

# Add already-done 5N7V result
all_results.append({
    "pdb_id":          "5N7V",
    "ligand":          "MPI-0479605",
    "scaffold":        "Purine",
    "sup_rmsd":        0.423,
    "best_cross_rmsd": 1.332,
    "pass":            True,
})

print("\n" + "=" * 60)
print("CROSS-DOCKING VALIDATION")
print("=" * 60)

for target in TARGETS:
    if target["done"]:
        continue

    pdb_id   = target["pdb_id"]
    pdb_file = f'data/raw/{pdb_id}.pdb'

    print(f"\n── {pdb_id} ({target['ligand']}, "
          f"{target['scaffold']} scaffold) ──")

    # Load mobile structure
    mob_structure = parser.get_structure(pdb_id, pdb_file)

    # Get rigid Cα
    mob_rigid = get_rigid_ca(mob_structure, pdb_file)
    print(f"  Mobile rigid Cα: {len(mob_rigid)}")

    # Find common rigid residues
    common = sorted(
        set(ref_rigid.keys()) & set(mob_rigid.keys())
    )
    print(f"  Common rigid residues: {len(common)}")

    if len(common) < 10:
        print(f"  ⚠ Too few common residues — skipping {pdb_id}")
        continue

    # Superimpose
    ref_ca_list = [ref_rigid[r] for r in common]
    mob_ca_list = [mob_rigid[r] for r in common]

    sup = Superimposer()
    sup.set_atoms(ref_ca_list, mob_ca_list)
    sup.apply(mob_structure.get_atoms())
    print(f"  Superimposition RMSD: {sup.rms:.3f} Å")

    # Save aligned receptor
    class ProteinSelect(Select):
        def accept_residue(self, residue):
            return residue.get_resname() in STANDARD_AA

    aligned_pdb  = f'data/receptor/{pdb_id}_aligned.pdb'
    aligned_pdbqt = f'data/receptor/{pdb_id}_aligned_receptor'

    io = PDBIO()
    io.set_structure(mob_structure)
    io.save(aligned_pdb, ProteinSelect())

    # Prepare PDBQT
    result = subprocess.run([
        'mk_prepare_receptor.py',
        '-i', aligned_pdb,
        '-o', aligned_pdbqt,
        '-p'
    ], capture_output=True, text=True)

    if result.returncode != 0:
        print(f"  ✗ PDBQT preparation failed")
        continue
    print(f"  ✓ Receptor prepared")

    # Dock reversine
    v = Vina(sf_name='vina')
    v.set_receptor(f'{aligned_pdbqt}.pdbqt')
    v.compute_vina_maps(
        center=[-34.48, -15.66, -10.38],
        box_size=[20, 33, 21]
    )
    v.set_ligand_from_file(
        'data/ligands/pdbqt/5LJJ_AD5_redock.pdbqt'
    )
    v.dock(exhaustiveness=32, n_poses=10)

    out_pdbqt = (f'docking/results/'
                 f'reversine_crossdock_{pdb_id}_out.pdbqt')
    v.write_poses(out_pdbqt, n_poses=10, overwrite=True)

    score = v.energies(n_poses=1)[0][0]
    print(f"  Best docking score: {score:.3f} kcal/mol")

    # Calculate RMSD
    poses              = parse_pdbqt_poses(out_pdbqt)
    best_rmsd, all_rmsds = calc_best_rmsd(
        poses, crystal_coords
    )

    passed = best_rmsd < 2.0
    print(f"  Pose RMSDs: "
          f"{[f'{r:.2f}' for r in all_rmsds]}")
    print(f"  Best RMSD: {best_rmsd:.3f} Å  "
          f"{'✓ PASS' if passed else '✗ FAIL'}")

    all_results.append({
        "pdb_id":          pdb_id,
        "ligand":          target["ligand"],
        "scaffold":        target["scaffold"],
        "sup_rmsd":        round(sup.rms, 3),
        "best_cross_rmsd": best_rmsd,
        "pass":            passed,
    })

# ── Summary ───────────────────────────────────────────────
print("\n" + "=" * 60)
print("CROSS-DOCKING VALIDATION SUMMARY")
print("=" * 60)
print(f"{'PDB':<6} {'Scaffold':<20} {'Sup RMSD':>9} "
      f"{'Cross RMSD':>11} {'Result':>8}")
print("-" * 60)

rmsds_passed = []
for r in all_results:
    result_str = "✓ PASS" if r["pass"] else "✗ FAIL"
    print(f"{r['pdb_id']:<6} {r['scaffold']:<20} "
          f"{r['sup_rmsd']:>9.3f} "
          f"{r['best_cross_rmsd']:>11.3f} "
          f"{result_str:>8}")
    if r["pass"]:
        rmsds_passed.append(r["best_cross_rmsd"])

all_rmsds_val = [r["best_cross_rmsd"] for r in all_results]
n_pass        = sum(r["pass"] for r in all_results)

print("-" * 60)
print(f"\nSelf-docking RMSD:         0.666 Å  ✓ PASS")
print(f"Cross-docking mean RMSD:   "
      f"{np.mean(all_rmsds_val):.3f} ± "
      f"{np.std(all_rmsds_val):.3f} Å")
print(f"Cross-docking passed:      {n_pass}/{len(all_results)}")

# ── Save results ──────────────────────────────────────────
df = pd.DataFrame(all_results)
df.to_csv(
    OUTPUT_DIR / "crossdocking_validation.csv", index=False
)
print(f"\nResults saved → "
      f"{OUTPUT_DIR}/crossdocking_validation.csv")

# ── Report text ───────────────────────────────────────────
print(f"\nMethods text:")
print(f"  'Cross-docking validation of reversine across")
print(f"   {len(all_results)} structurally diverse Mps1 crystal")
print(f"   structures (representing {len(all_results)} different")
print(f"   inhibitor scaffolds) yielded a mean RMSD of")
print(f"   {np.mean(all_rmsds_val):.3f} ± "
      f"{np.std(all_rmsds_val):.3f} Å,")
print(f"   with {n_pass}/{len(all_results)} structures passing")
print(f"   the 2.0 Å threshold.'")
