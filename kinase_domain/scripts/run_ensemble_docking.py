"""
Ensemble docking for Phase 2B candidates.
Docks each ligand into all 7 ensemble receptors and
computes consensus scores.
"""

from pathlib import Path
from vina import Vina
import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser(
    description="Ensemble docking with AutoDock Vina"
)
parser.add_argument("--ligands",
    default="kinase_domain/data/phase2b/pdbqt")
parser.add_argument("--receptors",
    default="kinase_domain/data/receptor/ensemble")
parser.add_argument("--results",
    default="kinase_domain/docking/phase2b_results")
parser.add_argument("--exhaustiveness",
    type=int, default=8)
parser.add_argument("--scores_csv",
    default=None)
args = parser.parse_args()

LIGANDS_DIR  = Path(args.ligands)
RECEPTOR_DIR = Path(args.receptors)
RESULTS_DIR  = Path(args.results)
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

SCORES_CSV = (
    args.scores_csv if args.scores_csv
    else str(RESULTS_DIR / "docking_scores.csv")
)

# Grid box — same as Phase 1 kinase domain
CENTER   = [-34.48, -15.66, -10.38]
BOX_SIZE = [20, 33, 21]

# Load ensemble receptor info
ensemble_csv = RECEPTOR_DIR / "ensemble_info.csv"
ensemble_df  = pd.read_csv(ensemble_csv)
ensemble_df  = ensemble_df[ensemble_df['prepared'] == True]

print("Ensemble Docking — Phase 2B")
print("=" * 55)
print(f"Ligands:    {LIGANDS_DIR}")
print(f"Receptors:  {len(ensemble_df)} ensemble members")
print(f"Results:    {RESULTS_DIR}")
print(f"Exhaustiveness: {args.exhaustiveness}\n")

print("Ensemble members:")
for _, row in ensemble_df.iterrows():
    print(f"  {row['pdb_id']}: quality={row['quality_score']} "
          f"res={row['resolution']} A")

pdbqt_files = sorted(LIGANDS_DIR.glob("*.pdbqt"))
print(f"\nFound {len(pdbqt_files)} ligands to dock\n")
print("=" * 55)

all_results = []

for lig_idx, ligand_path in enumerate(pdbqt_files):
    ligand_name = ligand_path.stem
    print(f"\n[{lig_idx+1}/{len(pdbqt_files)}] {ligand_name}")

    receptor_scores = {}

    for _, rec_row in ensemble_df.iterrows():
        pdb_id   = rec_row['pdb_id']
        rec_path = rec_row['pdbqt']

        if not Path(rec_path).exists():
            continue

        try:
            v = Vina(sf_name='vina', verbosity=0)
            v.set_receptor(rec_path)
            v.compute_vina_maps(
                center=CENTER,
                box_size=BOX_SIZE
            )
            v.set_ligand_from_file(str(ligand_path))
            v.dock(
                exhaustiveness=args.exhaustiveness,
                n_poses=5
            )

            # Save pose
            out_path = (
                RESULTS_DIR /
                f"{ligand_name}_{pdb_id}_out.pdbqt"
            )
            v.write_poses(
                str(out_path), n_poses=5, overwrite=True
            )

            score = v.energies(n_poses=1)[0][0]
            receptor_scores[pdb_id] = round(score, 3)
            print(f"  {pdb_id}: {score:.3f} kcal/mol")

        except Exception as e:
            print(f"  {pdb_id}: failed — {e}")
            receptor_scores[pdb_id] = None

    # Compute consensus scores
    valid_scores = [
        s for s in receptor_scores.values()
        if s is not None
    ]

    if not valid_scores:
        print(f"  All receptors failed — skipping")
        continue

    best_score  = min(valid_scores)
    mean_score  = np.mean(valid_scores)
    std_score   = np.std(valid_scores)
    n_receptors = len(valid_scores)

    # Consensus: how many receptors score < -7.0
    n_good = sum(1 for s in valid_scores if s < -7.0)

    result = {
        "ligand":         ligand_name,
        "best_score":     round(best_score, 3),
        "mean_score":     round(mean_score, 3),
        "std_score":      round(std_score, 3),
        "n_receptors":    n_receptors,
        "n_good":         n_good,
        "consensus_good": n_good >= 4,
    }

    for pdb_id, score in receptor_scores.items():
        result[f"score_{pdb_id}"] = score

    all_results.append(result)
    print(f"  Best: {best_score:.3f} | "
          f"Mean: {mean_score:.3f} | "
          f"Good receptors: {n_good}/{n_receptors}")

# Save results
results_df = pd.DataFrame(all_results)
results_df = results_df.sort_values(
    "best_score", ascending=True
)
results_df.to_csv(SCORES_CSV, index=False)

print(f"\n{'=' * 55}")
print(f"ENSEMBLE DOCKING COMPLETE")
print(f"{'=' * 55}")
print(f"  Ligands docked:     {len(results_df)}")
print(f"  Consensus binders:  "
      f"{results_df['consensus_good'].sum()}")
print(f"\n  Top 10 by best score:")
print(f"  {'Ligand':<35} {'Best':>7} {'Mean':>7} "
      f"{'Std':>6} {'Good':>5}")
print("  " + "-" * 65)

for _, r in results_df.head(10).iterrows():
    name = r['ligand'][:33]
    print(f"  {name:<35} "
          f"{r['best_score']:>7.3f} "
          f"{r['mean_score']:>7.3f} "
          f"{r['std_score']:>6.3f} "
          f"{r['n_good']:>5}")

print(f"\nResults → {SCORES_CSV}")
