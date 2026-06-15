# Usage :
#   Phase 1 :
#   Phase 2 : python scripts/run_docking.py \
#               --ligands        data/phase2/pdbqt \
#               --results        docking/phase2_results \
#               --exhaustiveness 8


from pathlib import Path
from vina import Vina
import pandas as pd
import argparse


# Argument parser
parser = argparse.ArgumentParser(
    description="Run AutoDock Vina docking"
)
parser.add_argument("--ligands",      
                    default="kinase_domain/data/ligands/pdbqt")
parser.add_argument("--results",      
                    default="kinase_domain/docking/results")
parser.add_argument("--exhaustiveness",
                    type=int, default=16)
parser.add_argument("--scores_csv",   
                    default=None,
                    help="Output CSV (auto if not set)")
args = parser.parse_args()


# PATHS
RECEPTOR_PATH = "kinase_domain/data/receptor/receptor.pdbqt"
LIGANDS_DIR = Path(args.ligands)
RESULTS_DIR = Path(args.results)
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

SCORES_CSV  = (args.scores_csv if args.scores_csv
               else str(RESULTS_DIR / "docking_scores.csv"))

# Results collector
results = []

# Initialise Vina with receptor
v = Vina (sf_name='vina')
v.set_receptor(RECEPTOR_PATH)

# Coordinates obtained through PyMOL using centerofmass and visual verification for grid boxing
v.compute_vina_maps(
    center=[-34.48, -15.66, -10.38],
    box_size=[20, 33, 21]
)

# Loop over all 45 ligands
pdbqt_files = sorted(LIGANDS_DIR.glob("*.pdbqt"))
print(f"Found {len(pdbqt_files)} ligands to dock\n")

for ligand_path in pdbqt_files:
    ligand_name = ligand_path.stem
    output_path = RESULTS_DIR / f"{ligand_name}_out.pdbqt"

    print(f"Docking: {ligand_name}")
    try:
        v.set_ligand_from_file(str(ligand_path))
        v.dock(exhaustiveness=args.exhaustiveness, n_poses=5)
        v.write_poses(str(output_path), n_poses=5, overwrite=True)

        # Extract best score (pose 1)
        energies = v.energies(n_poses=1)
        best_score = energies[0][0]

        print(f"   Best affinity: {best_score} kcal/mol")
        results.append({
            "ligand":          ligand_name,
            "affinity_kcal":   best_score,
            "output_file":     str(output_path)
        })

    except Exception as e:
        print(f"   Failed: {e}")
        results.append({
            "ligand":        ligand_name,
            "affinity_kcal": None,
            "output_file":   None
        })

# Save results to CSV 
results_df = pd.DataFrame(results)
results_df = results_df.sort_values("affinity_kcal")  # best (most negative) first
results_df.to_csv(SCORES_CSV, index=False)

print(f"\n{'='*50}")
print(f" Docking complete — {len(results_df)} ligands processed")
print(f"Best binder: {results_df.iloc[0]['ligand']} "
      f"({results_df.iloc[0]['affinity_kcal']} kcal/mol)")
print(f"Results saved to docking/results/docking_scores.csv")

