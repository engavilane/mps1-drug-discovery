"""
NTD TPR domain docking campaign.

Docks a fragment/lead-like library into the Hec1 interface
groove of Mps1 TPR domain (4H7Y).

Binding site: Pocket 42 (fpocket), validated by Screpanti 2011
Grid: center=[45.64, 27.59, 88.97], size=30x30x30 A
"""

from pathlib import Path
from vina import Vina
import pandas as pd
import numpy as np
import argparse
import warnings
warnings.filterwarnings('ignore')

parser = argparse.ArgumentParser(
    description="NTD TPR domain docking"
)
parser.add_argument("--ligands",
    default="ntd_domain/data/libraries/pdbqt")
parser.add_argument("--receptor",
    default="ntd_domain/data/receptor/4H7Y_receptor.pdbqt")
parser.add_argument("--results",
    default="ntd_domain/docking/results")
parser.add_argument("--exhaustiveness",
    type=int, default=16)
parser.add_argument("--scores_csv",
    default=None)
args = parser.parse_args()

LIGANDS_DIR  = Path(args.ligands)
RESULTS_DIR  = Path(args.results)
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

SCORES_CSV = (
    args.scores_csv if args.scores_csv
    else str(RESULTS_DIR / "ntd_docking_scores.csv")
)

# NTD grid — Pocket 42, Hec1 interface
CENTER   = [45.64, 27.59, 88.97]
BOX_SIZE = [30, 30, 30]

print("NTD TPR Domain Docking Campaign")
print("=" * 55)
print(f"Receptor:  {args.receptor}")
print(f"Ligands:   {LIGANDS_DIR}")
print(f"Grid:      centre={CENTER}")
print(f"           box={BOX_SIZE} A")
print(f"Exhaustiveness: {args.exhaustiveness}\n")

# Initialise Vina
v = Vina(sf_name='vina', verbosity=0)
v.set_receptor(args.receptor)
v.compute_vina_maps(
    center=CENTER,
    box_size=BOX_SIZE
)

pdbqt_files = sorted(LIGANDS_DIR.glob("*.pdbqt"))
print(f"Found {len(pdbqt_files)} ligands to dock\n")
print("=" * 55)

results = []

for idx, lig_path in enumerate(pdbqt_files):
    lig_name = lig_path.stem
    out_path = RESULTS_DIR / f"{lig_name}_out.pdbqt"

    try:
        v.set_ligand_from_file(str(lig_path))
        v.dock(
            exhaustiveness=args.exhaustiveness,
            n_poses=5
        )
        v.write_poses(
            str(out_path), n_poses=5, overwrite=True
        )

        score = v.energies(n_poses=1)[0][0]

        results.append({
            'ligand':        lig_name,
            'affinity_kcal': round(score, 3),
            'output_file':   str(out_path),
        })

        print(f"  [{idx+1}/{len(pdbqt_files)}] "
              f"{lig_name[:40]:<40} "
              f"{score:.3f} kcal/mol")

    except Exception as e:
        print(f"  [{idx+1}/{len(pdbqt_files)}] "
              f"{lig_name}: failed — {e}")
        results.append({
            'ligand':        lig_name,
            'affinity_kcal': None,
            'output_file':   None,
        })

results_df = pd.DataFrame(results)
results_df = results_df.dropna(subset=['affinity_kcal'])
results_df = results_df.sort_values(
    'affinity_kcal', ascending=True
)
results_df.to_csv(SCORES_CSV, index=False)

print(f"\n{'=' * 55}")
print(f"NTD DOCKING COMPLETE")
print(f"{'=' * 55}")
print(f"  Ligands docked: {len(results_df)}")
print(f"  Best score:     "
      f"{results_df['affinity_kcal'].min():.3f} kcal/mol")
print(f"  Best binder:    "
      f"{results_df.iloc[0]['ligand']}")

print(f"\nTop 10 NTD binders:")
print(f"{'Ligand':<40} {'Score':>8}")
print("-" * 50)
for _, r in results_df.head(10).iterrows():
    print(f"  {r['ligand']:<38} "
          f"{r['affinity_kcal']:>8.3f}")

print(f"\nResults → {SCORES_CSV}")
