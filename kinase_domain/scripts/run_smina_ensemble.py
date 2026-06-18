"""
Smina ensemble docking for Phase 2B candidates.

Uses 3 most representative receptor conformations:
  5AP7 — highest quality (9.6), holo
  5AP6 — second highest quality (9.4), holo
  3CEK — best apo structure (9.0)

Smina uses Vina scoring — fast CPU docking.
Top 50 candidates will be rescored with GNINA CNN on Colab.
"""

import subprocess
import pandas as pd
import numpy as np
from pathlib import Path
import argparse
import time
import warnings
warnings.filterwarnings('ignore')

parser = argparse.ArgumentParser(
    description="Smina ensemble docking"
)
parser.add_argument("--ligands",
    default="kinase_domain/data/phase2b/pdbqt")
parser.add_argument("--receptors",
    default="kinase_domain/data/receptor/ensemble")
parser.add_argument("--results",
    default="kinase_domain/docking/phase2b_results")
parser.add_argument("--exhaustiveness",
    type=int, default=8)
args = parser.parse_args()

LIGANDS_DIR  = Path(args.ligands)
RECEPTOR_DIR = Path(args.receptors)
RESULTS_DIR  = Path(args.results)
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

CENTER   = [-34.48, -15.66, -10.38]
BOX_SIZE = [20, 33, 21]

SELECTED_RECEPTORS = ['5AP7']

receptors = [
    r for r in sorted(RECEPTOR_DIR.glob('*.pdbqt'))
    if any(s in r.stem for s in SELECTED_RECEPTORS)
]

ligands = sorted(LIGANDS_DIR.glob('*.pdbqt'))

print("Smina Ensemble Docking — Phase 2B")
print("=" * 55)
print(f"Ligands:        {len(ligands)}")
print(f"Receptors:      {len(receptors)}")
for r in receptors:
    print(f"  {r.stem}")
print(f"Total runs:     {len(ligands) * len(receptors)}")
print(f"Exhaustiveness: {args.exhaustiveness}\n")

results = []
start   = time.time()

for lig_idx, lig_path in enumerate(ligands):
    lig_name   = lig_path.stem
    rec_scores = {}

    for rec_path in receptors:
        rec_name = rec_path.stem.replace(
            '_ensemble_receptor', ''
        )
        out_path = (
            RESULTS_DIR /
            f"{lig_name}_{rec_name}_out.pdbqt"
        )

        cmd = [
            'smina',
            '--receptor',       str(rec_path),
            '--ligand',         str(lig_path),
            '--out',            str(out_path),
            '--center_x',       str(CENTER[0]),
            '--center_y',       str(CENTER[1]),
            '--center_z',       str(CENTER[2]),
            '--size_x',         str(BOX_SIZE[0]),
            '--size_y',         str(BOX_SIZE[1]),
            '--size_z',         str(BOX_SIZE[2]),
            '--exhaustiveness', str(args.exhaustiveness),
            '--num_modes',      '5',
        ]

        try:
            r = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=120
            )

            score = None
            for line in r.stdout.split('\n'):
                parts = line.split()
                if len(parts) >= 2 and parts[0] == '1':
                    try:
                        score = float(parts[1])
                        break
                    except:
                        pass

            rec_scores[rec_name] = score

        except Exception as e:
            rec_scores[rec_name] = None

    valid = [s for s in rec_scores.values()
             if s is not None]

    if valid:
        result = {
            'ligand':      lig_name,
            'best_score':  round(min(valid), 3),
            'mean_score':  round(np.mean(valid), 3),
            'std_score':   round(np.std(valid), 3),
            'n_receptors': len(valid),
            'n_good':      sum(1 for s in valid
                               if s < -7.0),
            'consensus':   sum(1 for s in valid
                               if s < -7.0) >= 2,
        }
        for rec in SELECTED_RECEPTORS:
            result[f'score_{rec}'] = rec_scores.get(rec)

        results.append(result)

        elapsed = time.time() - start
        eta     = (elapsed/(lig_idx+1)) * (
            len(ligands)-lig_idx-1
        )
        print(f"  [{lig_idx+1}/{len(ligands)}] "
              f"{lig_name[:35]:<35} "
              f"best={result['best_score']:.3f} "
              f"ETA: {eta/60:.0f} min")

    else:
        print(f"  [{lig_idx+1}/{len(ligands)}] "
              f"{lig_name}: all receptors failed")

results_df = pd.DataFrame(results)
results_df = results_df.sort_values(
    'best_score', ascending=True
)

scores_path = RESULTS_DIR / "docking_scores.csv"
results_df.to_csv(scores_path, index=False)

print(f"\n{'=' * 55}")
print(f"SMINA ENSEMBLE DOCKING COMPLETE")
print(f"{'=' * 55}")
print(f"  Ligands docked:     {len(results_df)}")
print(f"  Consensus binders:  "
      f"{results_df['consensus'].sum()}")
print(f"\n  Top 10 by best score:")
print(f"  {'Ligand':<35} {'Best':>7} "
      f"{'Mean':>7} {'Good':>5}")
print("  " + "-" * 58)

for _, r in results_df.head(10).iterrows():
    name = r['ligand'][:33]
    print(f"  {name:<33} "
          f"{r['best_score']:>7.3f} "
          f"{r['mean_score']:>7.3f} "
          f"{r['n_good']:>5}")

print(f"\nResults → {scores_path}")
