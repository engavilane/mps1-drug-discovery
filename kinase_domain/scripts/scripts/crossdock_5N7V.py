from Bio.PDB import PDBParser, PDBIO, Select, Superimposer
from Bio.PDB.DSSP import DSSP
import numpy as np
import subprocess
import warnings
warnings.filterwarnings('ignore')

# ── Configuration ─────────────────────────────────────────
STANDARD_AA = [
    'ALA','ARG','ASN','ASP','CYS','GLN','GLU',
    'GLY','HIS','ILE','LEU','LYS','MET','PHE',
    'PRO','SER','THR','TRP','TYR','VAL'
]

# Flexible regions to exclude — known from Mps1 literature
# Bolanos-Garcia 2025, Pugh et al. 2022
EXCLUDE_RESIDUES = set(range(672, 691))  # activation + P+1 loop
EXCLUDE_RESIDUES.update(range(603, 606)) # hinge loop

# Secondary structure types considered rigid
# H=alpha helix, G=3-10 helix, I=pi helix, E=beta strand
RIGID_SS = {'H', 'G', 'I', 'E'}

# ── Load structures ───────────────────────────────────────
print("Loading structures...")
parser = PDBParser(QUIET=True)
ref_structure = parser.get_structure(
    '5LJJ', 'data/raw/5LJJ.pdb'
)
mob_structure = parser.get_structure(
    '5N7V', 'data/raw/5N7V.pdb'
)
print("  ✓ 5LJJ and 5N7V loaded")

# ── Run DSSP on both structures ───────────────────────────
print("\nRunning DSSP secondary structure assignment...")

def run_dssp(structure, pdb_file):
    model = list(structure.get_models())[0]
    try:
        dssp = DSSP(model, pdb_file, dssp='mkdssp')
        return {
            (k[0], k[1][1]): v[2]   # (chain, resid) → ss_type
            for k, v in dssp.property_dict.items()
        }
    except Exception as e:
        print(f"  DSSP error: {e}")
        return {}

ref_dssp = run_dssp(ref_structure, 'data/raw/5LJJ.pdb')
mob_dssp = run_dssp(mob_structure, 'data/raw/5N7V.pdb')

print(f"  5LJJ: {len(ref_dssp)} residues assigned")
print(f"  5N7V: {len(mob_dssp)} residues assigned")

# ── Get Cα atoms filtered by DSSP ────────────────────────
def get_rigid_ca(structure, dssp_dict):
    """
    Returns dict of resid → CA atom
    Only rigid secondary structure elements,
    excluding known flexible regions
    """
    residues = {}
    for model in structure:
        for chain in model:
            for residue in chain:
                resid   = residue.get_id()[1]
                resname = residue.get_resname()
                chain_id = chain.get_id()

                # Standard amino acids only
                if resname not in STANDARD_AA:
                    continue

                # Exclude known flexible regions
                if resid in EXCLUDE_RESIDUES:
                    continue

                # Must have Cα
                if 'CA' not in residue:
                    continue

                # Check DSSP secondary structure
                ss = dssp_dict.get((chain_id, resid), 'C')

                # Keep only rigid secondary structure
                if ss not in RIGID_SS:
                    continue

                residues[resid] = {
                    'CA': residue['CA'],
                    'ss': ss,
                    'resname': resname
                }
    return residues

print("\nSelecting rigid Cα atoms...")
ref_rigid = get_rigid_ca(ref_structure, ref_dssp)
mob_rigid = get_rigid_ca(mob_structure, mob_dssp)

print(f"  5LJJ rigid residues: {len(ref_rigid)}")
print(f"  5N7V rigid residues: {len(mob_rigid)}")

# ── Find common rigid residues ────────────────────────────
common = sorted(
    set(ref_rigid.keys()) & set(mob_rigid.keys())
)

# Keep only residues where both agree on SS type
common_same_ss = [
    r for r in common
    if ref_rigid[r]['ss'] == mob_rigid[r]['ss']
]

print(f"  Common rigid residues:           {len(common)}")
print(f"  Common with same SS assignment:  {len(common_same_ss)}")

if len(common_same_ss) < 10:
    print("  ⚠ Too few — falling back to all common Cα")
    common_same_ss = sorted(
        set(r for r in ref_rigid) &
        set(r for r in mob_rigid)
    )

# ── Superimpose on rigid core ─────────────────────────────
print("\nSuperimposing on rigid secondary structure core...")
ref_ca_list = [ref_rigid[r]['CA'] for r in common_same_ss
               if r in ref_rigid]
mob_ca_list = [mob_rigid[r]['CA'] for r in common_same_ss
               if r in mob_rigid]

sup = Superimposer()
sup.set_atoms(ref_ca_list, mob_ca_list)
sup.apply(mob_structure.get_atoms())

print(f"  Residues used:           {len(common_same_ss)}")
print(f"  Superimposition RMSD:    {sup.rms:.3f} Å")
print(f"  (helices + sheets only,")
print(f"   excluding activation loop 672-690")
print(f"   and hinge loop 603-605)")

# ── Save aligned receptor ─────────────────────────────────
class ProteinSelect(Select):
    def accept_residue(self, residue):
        return residue.get_resname() in STANDARD_AA

io = PDBIO()
io.set_structure(mob_structure)
io.save('data/receptor/5N7V_aligned_dssp.pdb', ProteinSelect())
print("\n  ✓ Aligned 5N7V saved (DSSP rigid core)")

# ── Prepare PDBQT ─────────────────────────────────────────
result = subprocess.run([
    'mk_prepare_receptor.py',
    '-i', 'data/receptor/5N7V_aligned_dssp.pdb',
    '-o', 'data/receptor/5N7V_aligned_dssp_receptor',
    '-p'
], capture_output=True, text=True)

if result.returncode == 0:
    print("  ✓ Receptor PDBQT prepared")
else:
    print(f"  ✗ Error: {result.stderr[:200]}")
    exit()

# ── Cross-dock reversine ──────────────────────────────────
print("\nCross-docking reversine into aligned 5N7V...")
from vina import Vina

v = Vina(sf_name='vina')
v.set_receptor(
    'data/receptor/5N7V_aligned_dssp_receptor.pdbqt'
)
v.compute_vina_maps(
    center=[-34.48, -15.66, -10.38],
    box_size=[20, 33, 21]
)
v.set_ligand_from_file(
    'data/ligands/pdbqt/5LJJ_AD5_redock.pdbqt'
)
v.dock(exhaustiveness=32, n_poses=10)
v.write_poses(
    'docking/results/5LJJ_AD5_crossdock_5N7V_dssp_out.pdbqt',
    n_poses=10, overwrite=True
)
score = v.energies(n_poses=1)[0][0]
print(f"  Best cross-docking score: {score:.3f} kcal/mol")

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

poses   = []
current = []
pdbqt   = ('docking/results/'
            '5LJJ_AD5_crossdock_5N7V_dssp_out.pdbqt')

with open(pdbqt) as f:
    for line in f:
        if line.startswith('MODEL'):
            current = []
        elif line.startswith('ENDMDL'):
            poses.append(current)
        elif line.startswith(('ATOM', 'HETATM')):
            current.append(line)

print(f"\nCROSS-DOCKING RMSD")
print(f"(reversine into DSSP-aligned 5N7V")
print(f" vs 5LJJ crystal pose)")
print("=" * 55)

results = []
for i, pose_lines in enumerate(poses):
    docked = pdbqt_to_coords(pose_lines)
    if len(docked) == 0:
        continue
    n = min(len(crystal_coords), len(docked))
    dist_matrix      = cdist(crystal_coords[:n], docked[:n])
    row_ind, col_ind = linear_sum_assignment(dist_matrix)
    diff = crystal_coords[row_ind] - docked[col_ind]
    rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
    results.append(rmsd)
    star = "⭐" if rmsd < 1.0 else ""
    print(f"  Pose {i+1:2d}: RMSD = {rmsd:.3f} Å"
          f"  {'✓ PASS' if rmsd < 2.0 else '✗ FAIL'} {star}")

best = min(results)
print(f"\n{'=' * 55}")
print(f"Structural superimposition RMSD: {sup.rms:.3f} Å")
print(f"Best cross-docking RMSD:         {best:.3f} Å")
print(f"Self-docking RMSD (reference):   0.666 Å")

if best < 2.0:
    print("\n✓ CROSS-DOCKING VALIDATION PASSED")
    print("  Vina correctly places reversine even in a")
    print("  receptor from a different crystal structure")
elif best < 3.0:
    print(f"\n⚠ BORDERLINE ({best:.3f} Å)")
    print(f"  Structural difference between receptors:")
    print(f"  {sup.rms:.3f} Å — some deviation expected")
    print(f"  with rigid docking across conformations")
else:
    print(f"\n✗ CROSS-DOCKING FAILED ({best:.3f} Å)")
    print(f"  Structural difference between receptors:")
    print(f"  {sup.rms:.3f} Å")
    print(f"  Rigid receptor limitation — expected result")
    print(f"  Self-docking (0.666 Å) remains valid")
