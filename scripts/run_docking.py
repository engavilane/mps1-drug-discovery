from pathlib import Path
from vina import Vina
import pandas as pd


# PATHS
receptor_path = 'data/receptor/receptor.pdbqt'
ligands_dir = Path('data/ligands/pdbqt')
results_dir = Path ('docking/results')
results_dir.mkdir(parents=True, exist_ok=True)

# Results collector
results = []

# Initialise Vina with receptor
v = Vina (sf_name='vina')
v.set_receptor(receptor_path)

# Coordinates obtained through PyMOL using centerofmass and visual verification for grid boxing
v.compute_vina_maps(
    center=[-34.48, -15.66, -10.38],
    box_size=[20, 33, 21]
)

# Loop over all 45 ligands
pdbqt_files = sorted(ligands_dir.glob("*.pdbqt"))
print(f"Found {len(pdbqt_files)} ligands to dock\n")

for ligand_path in pdbqt_files:
    ligand_name = ligand_path.stem
    output_path = results_dir / f"{ligand_name}_out.pdbqt"

    print(f"Docking: {ligand_name}")
    try:
        v.set_ligand_from_file(str(ligand_path))
        v.dock(exhaustiveness=16, n_poses=5)
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
results_df.to_csv("docking/results/docking_scores.csv", index=False)

print(f"\n{'='*50}")
print(f" Docking complete — {len(results_df)} ligands processed")
print(f"Best binder: {results_df.iloc[0]['ligand']} "
      f"({results_df.iloc[0]['affinity_kcal']} kcal/mol)")
print(f"Results saved to docking/results/docking_scores.csv")

