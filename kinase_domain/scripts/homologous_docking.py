"""
Homologous docking validation — top 10 Phase 1 candidates
docked into their own co-crystallised receptor conformation.

Tests whether the 5LJJ-based ranking introduces systematic
single-receptor conformation bias.

If Spearman correlation > 0.7 → ranking is reliable
If Spearman correlation < 0.5 → significant bias detected
"""

from Bio.PDB import PDBParser, PDBIO, Select, Superimposer
from Bio.PDB.DSSP import DSSP
from scipy.spatial.distance import cdist
from scipy.optimize import linear_sum_assignment
from scipy.stats import spearmanr
from vina import Vina
import numpy as np
import subprocess
import pandas as pd
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Configuration 
STANDARD_AA = [
    'ALA','ARG','ASN','ASP','CYS','GLN','GLU',
    'GLY','HIS','ILE','LEU','LYS','MET','PHE',
    'PRO','SER','THR','TRP','TYR','VAL'
]
EXCLUDE_RESIDUES = set(range(672, 691))
EXCLUDE_RESIDUES.update(range(603, 606))
RIGID_SS  = {'H', 'G', 'I', 'E'}
OUTPUT_DIR = Path("kinase_domain/analysis/validation")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# Top 10 compounds — ligand file → native PDB
TOP_10 = {
    "7CHN_LIG": {"pdb": "7CHN", "vina_score": -10.324},
    "4JS8_LIG": {"pdb": "4JS8", "vina_score": -10.316},
    "7CHM_LIG": {"pdb": "7CHM", "vina_score": -10.288},
    "7CJA_LIG": {"pdb": "7CJA", "vina_score": -10.281},
    "7CHT_LIG": {"pdb": "7CHT", "vina_score": -9.970},
    "5EI8_LIG": {"pdb": "5EI8", "vina_score": -9.867},
    "5EI6_LIG": {"pdb": "5EI6", "vina_score": -9.752},
    "4D2S_LIG": {"pdb": "4D2S", "vina_score": -9.692},
    "7CIL_LIG": {"pdb": "7CIL", "vina_score": -9.482},
    "4CV9_LIG": {"pdb": "4CV9", "vina_score": -9.478},
}


# Helper functions
parser = PDBParser(QUIET=True)

def get_rigid_ca(structure, pdb_file):
    model = list(structure.get_models())[0]
    try:
        dssp = DSSP(model, pdb_file, dssp='mkdssp')
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
                if resid in EXCLUDE_RESIDUES:
                    continue
                if 'CA' not in residue:
                    continue
                ss = dssp_dict.get((chain_id, resid), 'C')
                if ss not in RIGID_SS:
                    continue
                residues[resid] = residue['CA']
    return residues

class ProteinSelect(Select):
    def accept_residue(self, residue):
        return residue.get_resname() in STANDARD_AA


# Load reference (5LJJ) 
print("Loading reference structure (5LJJ)...")
ref_structure = parser.get_structure(
    '5LJJ', 'kinase_domain/data/raw/5LJJ.pdb'
)
ref_rigid = get_rigid_ca(
    ref_structure, 'kinase_domain/data/raw/5LJJ.pdb'
)
print(f"  Reference rigid Cα: {len(ref_rigid)}\n")


# Main loop 
print("=" * 65)
print("HOMOLOGOUS DOCKING — TOP 10 CANDIDATES")
print("Each compound docked into its own crystal structure")
print("=" * 65)

results = []

for ligand_name, info in TOP_10.items():
    pdb_id     = info["pdb"]
    vina_5ljj  = info["vina_score"]
    pdb_file   = f'kinase_domain/data/raw/{pdb_id}.pdb'
    ligand_sdf = f'kinase_domain/data/ligands/raw/{ligand_name}.sdf'
    ligand_pdbqt = f'kinase_domain/data/ligands/pdbqt/{ligand_name}.pdbqt'

    print(f"\n── {ligand_name} → {pdb_id} ──")
    print(f"   5LJJ score: {vina_5ljj} kcal/mol")

    # Check files exist
    if not Path(pdb_file).exists():
        print(f"   ✗ PDB not found: {pdb_file}")
        continue
    if not Path(ligand_pdbqt).exists():
        print(f"   ✗ Ligand PDBQT not found: {ligand_pdbqt}")
        continue

    # Load and align native receptor 
    mob_structure = parser.get_structure(pdb_id, pdb_file)
    mob_rigid     = get_rigid_ca(mob_structure, pdb_file)

    common = sorted(
        set(ref_rigid.keys()) & set(mob_rigid.keys())
    )

    if len(common) < 10:
        print(f"   ✗ Too few common Cα: {len(common)}")
        continue

    ref_ca_list = [ref_rigid[r] for r in common]
    mob_ca_list = [mob_rigid[r] for r in common]

    sup = Superimposer()
    sup.set_atoms(ref_ca_list, mob_ca_list)
    sup.apply(mob_structure.get_atoms())
    print(f"   Superimposition RMSD: {sup.rms:.3f} Å "
          f"({len(common)} Cα)")


    # Save aligned receptor 
    aligned_pdb  = f'kinase_domain/data/receptor/{pdb_id}_homo_aligned.pdb'
    aligned_pdbqt = f'kinase_domain/data/receptor/{pdb_id}_homo_receptor'

    io = PDBIO()
    io.set_structure(mob_structure)
    io.save(aligned_pdb, ProteinSelect())


    # Prepare PDBQT 
    result = subprocess.run([
        'mk_prepare_receptor.py',
        '-i', aligned_pdb,
        '-o', aligned_pdbqt,
        '--allow_bad_res',
        '-p'
    ], capture_output=True, text=True)

    pdbqt_path = Path(f'{aligned_pdbqt}.pdbqt')
    if not pdbqt_path.exists():
        print(f"   ✗ Receptor preparation failed")
        print(f"   {result.stderr[-200:]}")
        continue
    print(f"   ✓ Receptor prepared")


    # Dock into native conformation 
    out_pdbqt = (f'kinase_domain/docking/results/'
                 f'{ligand_name}_homo_{pdb_id}_out.pdbqt')

    v = Vina(sf_name='vina')
    v.set_receptor(str(pdbqt_path))
    v.compute_vina_maps(
        center=[-34.48, -15.66, -10.38],
        box_size=[20, 33, 21]
    )
    v.set_ligand_from_file(ligand_pdbqt)
    v.dock(exhaustiveness=16, n_poses=5)
    v.write_poses(out_pdbqt, n_poses=5, overwrite=True)

    native_score = v.energies(n_poses=1)[0][0]
    print(f"   Native score:  {native_score:.3f} kcal/mol")
    print(f"   5LJJ score:    {vina_5ljj:.3f} kcal/mol")
    print(f"   Difference:    "
          f"{native_score - vina_5ljj:+.3f} kcal/mol")

    results.append({
        "ligand":        ligand_name,
        "native_pdb":    pdb_id,
        "sup_rmsd":      round(sup.rms, 3),
        "score_5ljj":    vina_5ljj,
        "score_native":  round(native_score, 3),
        "score_diff":    round(native_score - vina_5ljj, 3),
    })


# Summary and Spearman correlation 
print("\n" + "=" * 65)
print("HOMOLOGOUS DOCKING SUMMARY")
print("=" * 65)
print(f"{'Ligand':<12} {'PDB':>5} {'SupRMSD':>8} "
      f"{'5LJJ':>8} {'Native':>8} {'Diff':>7}")
print("-" * 65)

for r in results:
    print(f"{r['ligand']:<12} {r['native_pdb']:>5} "
          f"{r['sup_rmsd']:>8.3f} "
          f"{r['score_5ljj']:>8.3f} "
          f"{r['score_native']:>8.3f} "
          f"{r['score_diff']:>+7.3f}")

if len(results) >= 3:
    scores_5ljj   = [r["score_5ljj"]   for r in results]
    scores_native = [r["score_native"]  for r in results]
    score_diffs   = [r["score_diff"]    for r in results]

    rho, pval = spearmanr(scores_5ljj, scores_native)

    print("-" * 65)
    print(f"\nSpearman rank correlation: ρ = {rho:.3f} "
          f"(p = {pval:.4f})")
    print(f"Mean score difference:     "
          f"{np.mean(score_diffs):.3f} ± "
          f"{np.std(score_diffs):.3f} kcal/mol")
    print(f"Max score difference:      "
          f"{max(score_diffs, key=abs):.3f} kcal/mol")

    if rho >= 0.7:
        print("\n✓ RANKING VALIDATED — 5LJJ-based ranking")
        print("  is consistent with native conformation ranking")
        print(f"  (ρ = {rho:.3f} ≥ 0.7 threshold)")
    elif rho >= 0.5:
        print(f"\n⚠ MODERATE CORRELATION (ρ = {rho:.3f})")
        print("  Ranking broadly consistent but some")
        print("  compounds may be re-ordered")
    else:
        print(f"\n✗ WEAK CORRELATION (ρ = {rho:.3f})")
        print("  Single receptor conformation introduces")
        print("  significant ranking bias")


    # Save results 
    df = pd.DataFrame(results)
    df.to_csv(
        OUTPUT_DIR / "homologous_docking.csv", index=False
    )
    print(f"\nResults saved → "
          f"{OUTPUT_DIR}/homologous_docking.csv")


    # Methods text 
    print(f"\nMethods text:")
    print(f"  'To assess single-receptor conformation bias,")
    print(f"   the top {len(results)} Phase 1 candidates were")
    print(f"   re-docked into their own co-crystallised Mps1")
    print(f"   receptor conformations, aligned to 5LJJ via")
    print(f"   DSSP rigid core superimposition. Spearman rank")
    print(f"   correlation between 5LJJ-based and native")
    print(f"   conformation scores was ρ = {rho:.3f}")
    print(f"   (p = {pval:.4f}), indicating")
    if rho >= 0.7:
        print(f"   that the 5LJJ-based ranking is reliable.'")
    elif rho >= 0.5:
        print(f"   moderate ranking consistency.'")
    else:
        print(f"   potential ranking bias from single receptor use.'")
