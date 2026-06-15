# scripts/crossdock_5N7V.py
from Bio.PDB import PDBParser, PDBIO, Select, Superimposer
import numpy as np
import subprocess
import warnings
warnings.filterwarnings('ignore')

# ── Load both structures ──────────────────────────────────
parser = PDBParser(QUIET=True)
ref_structure = parser.get_structure(
    '5LJJ', 'data/raw/5LJJ.pdb'
)
mob_structure = parser.get_structure(
    '5N7V', 'data/raw/5N7V.pdb'
)

# ── Get Cα atoms for superimposition ─────────────────────
# Use hinge region residues 516-808 (common to both structures)
def get_ca_atoms(structure, res_range=(516, 808)):
    atoms = []
    for model in structure:
        for chain in model:
            for residue in chain:
                resid = residue.get_id()[1]
                if res_range[0] <= resid <= res_range[1]:
                    if 'CA' in residue:
                        atoms.append(residue['CA'])
    return atoms

print("Getting Cα atoms for superimposition...")
ref_atoms = get_ca_atoms(ref_structure)
mob_atoms = get_ca_atoms(mob_structure)

# Match by residue number
ref_dict = {a.get_parent().get_id()[1]: a for a in ref_atoms}
mob_dict = {a.get_parent().get_id()[1]: a for a in mob_atoms}
common   = sorted(set(ref_dict.keys()) & set(mob_dict.keys()))

ref_list = [ref_dict[r] for r in common]
mob_list = [mob_dict[r] for r in common]
print(f"  Common Cα atoms: {len(common)}")

# ── Superimpose ───────────────────────────────────────────
sup = Superimposer()
sup.set_atoms(ref_list, mob_list)
sup.apply(mob_structure.get_atoms())
print(f"  RMSD after superimposition: {sup.rms:.3f} Å")

# ── Save aligned 5N7V (protein only) ─────────────────────
class ProteinSelect(Select):
    def accept_residue(self, residue):
        standard = [
            'ALA','ARG','ASN','ASP','CYS','GLN','GLU',
            'GLY','HIS','ILE','LEU','LYS','MET','PHE',
            'PRO','SER','THR','TRP','TYR','VAL'
        ]
        return residue.get_resname() in standard

io = PDBIO()
io.set_structure(mob_structure)
io.save('data/receptor/5N7V_aligned.pdb', ProteinSelect())
print("  ✓ Aligned 5N7V saved")

# ── Prepare PDBQT ─────────────────────────────────────────
result = subprocess.run([
    'mk_prepare_receptor.py',
    '-i', 'data/receptor/5N7V_aligned.pdb',
    '-o', 'data/receptor/5N7V_aligned_receptor',
    '-p'
], capture_output=True, text=True)

if result.returncode == 0:
    print("  ✓ 5N7V aligned receptor.pdbqt prepared\n")
else:
    print(f"  ✗ Error: {result.stderr}")
    exit()

# ── Cross-dock reversine ──────────────────────────────────
print("Cross-docking reversine into aligned 5N7V...")
from vina import Vina

v = Vina(sf_name='vina')
v.set_receptor('data/receptor/5N7V_aligned_receptor.pdbqt')
v.compute_vina_maps(
    center=[-34.48, -15.66, -10.38],  # now valid after alignment
    box_size=[20, 33, 21]
)
v.set_ligand_from_file(
    'data/ligands/pdbqt/5LJJ_AD5_redock.pdbqt'
)
v.dock(exhaustiveness=32, n_poses=10)
v.write_poses(
    'docking/results/5LJJ_AD5_crossdock_5N7V_aligned_out.pdbqt',
    n_poses=10, overwrite=True
)

score = v.energies(n_poses=1)[0][0]
print(f"Best cross-docking score: {score:.3f} kcal/mol\n")

# ── RMSD calculation ──────────────────────────────────────
from scipy.spatial.distance import cdist
from scipy.optimize import linear_sum_assignment

crystal_coords = np.load(
    'data/receptor/reversine_crystal_coords.npy'
)

def pdbqt_to_coords(lines):
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

poses = []
current = []
pdbqt = 'docking/results/5LJJ_AD5_crossdock_5N7V_aligned_out.pdbqt'
with open(pdbqt) as f:
    for line in f:
        if line.startswith('MODEL'):
            current = []
        elif line.startswith('ENDMDL'):
            poses.append(current)
        elif line.startswith(('ATOM','HETATM')):
            current.append(line)

print("CROSS-DOCKING RMSD (reversine into aligned 5N7V):")
print("=" * 55)

results = []
for i, pose_lines in enumerate(poses):
    docked = pdbqt_to_coords(pose_lines)
    if len(docked) == 0:
        continue
    n = min(len(crystal_coords), len(docked))
    dist_matrix = cdist(crystal_coords[:n], docked[:n])
    row_ind, col_ind = linear_sum_assignment(dist_matrix)
    diff = crystal_coords[row_ind] - docked[col_ind]
    rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
    results.append(rmsd)
    star = "⭐" if rmsd < 1.0 else ""
    print(f"  Pose {i+1}: RMSD = {rmsd:.3f} Å"
          f"  {'✓ PASS' if rmsd < 2.0 else '✗ FAIL'} {star}")

best = min(results)
print(f"\nStructural superimposition RMSD: {sup.rms:.3f} Å")
print(f"Best cross-docking RMSD:         {best:.3f} Å")

if best < 2.0:
    print("✓ CROSS-DOCKING VALIDATION PASSED")
elif best < 3.0:
    print("⚠ BORDERLINE — good for rigid cross-docking")
else:
    print("✗ CROSS-DOCKING FAILED")
    print(f"  Note: structural difference between 5LJJ and 5N7V")
    print(f"  is {sup.rms:.3f} Å — receptor flexibility may explain result")
