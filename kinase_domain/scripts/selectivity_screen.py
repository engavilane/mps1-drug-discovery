"""
Selectivity counter-screen for top Mps1 candidates.

Docks top candidates into Aurora B (4TND) and computes
selectivity index:
  SI = IC50(Aurora B) / IC50(Mps1)
  SI > 100 = highly selective for Mps1

Uses Vina scores as proxy for IC50 (relative comparison only).
"""

from pathlib import Path
import subprocess
import pandas as pd
import numpy as np
import argparse
import warnings
warnings.filterwarnings('ignore')

parser = argparse.ArgumentParser(
    description="Selectivity counter-screen"
)
parser.add_argument("--candidates",
    default="kinase_domain/analysis/phase2b/phase2b_multitask_predictions.csv")
parser.add_argument("--ligands",
    default="kinase_domain/data/phase2b/pdbqt")
parser.add_argument("--receptor",
    default="kinase_domain/data/receptor/selectivity/4TND_receptor.pdbqt")
parser.add_argument("--results",
    default="kinase_domain/docking/selectivity_results")
parser.add_argument("--top_n",
    type=int, default=10)
parser.add_argument("--exhaustiveness",
    type=int, default=16)
args = parser.parse_args()

LIGANDS_DIR = Path(args.ligands)
RESULTS_DIR = Path(args.results)
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

CENTER   = [28.42, -38.19, 26.10]
BOX_SIZE = [22, 22, 22]

print("Selectivity Counter-Screen — Aurora B (4TND)")
print("=" * 55)
print(f"Top N:          {args.top_n}")
print(f"Exhaustiveness: {args.exhaustiveness}\n")

cand_df = pd.read_csv(args.candidates)
print(f"Candidates loaded: {len(cand_df)}")

score_col = (
    'mps1_pIC50_mean' if 'mps1_pIC50_mean' in cand_df.columns
    else 'pIC50_gnn'  if 'pIC50_gnn'       in cand_df.columns
    else 'affinity_kcal'
)

top = cand_df.nlargest(
    args.top_n, score_col
).reset_index(drop=True)

print(f"Top {args.top_n} by {score_col}")
print("=" * 55)

results = []

for idx, row in top.iterrows():
    ligand_name = row['ligand']

    lig_path = LIGANDS_DIR / f"{ligand_name}.pdbqt"
    if not lig_path.exists():
        clean = ligand_name.replace('_out', '')
        lig_path = LIGANDS_DIR / f"{clean}.pdbqt"
    if not lig_path.exists():
        print(f"  [{idx+1}/{args.top_n}] {ligand_name}: "
              f"PDBQT not found")
        continue

    mps1_score = row.get(
        'affinity_kcal',
        row.get('best_score', 0)
    )

    cmd = [
        'vina',
        '--receptor',       args.receptor,
        '--ligand',         str(lig_path),
        '--out',            str(RESULTS_DIR /
                            f"{ligand_name}_aurB_out.pdbqt"),
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
            cmd, capture_output=True,
            text=True, timeout=180
        )
        aurora_score = None
        for line in r.stdout.split('\n'):
            parts = line.split()
            if len(parts) >= 2 and parts[0] == '1':
                try:
                    aurora_score = float(parts[1])
                    break
                except:
                    pass

        if aurora_score is None:
            print(f"  [{idx+1}/{args.top_n}] "
                  f"{ligand_name}: docking failed")
            continue

        score_diff = aurora_score - mps1_score
        selective  = score_diff > 2.0

        result = {
            'ligand':       ligand_name,
            'mps1_score':   round(mps1_score, 3),
            'aurora_score': round(aurora_score, 3),
            'score_diff':   round(score_diff, 3),
            'selective':    selective,
        }

        if 'mps1_pIC50_mean' in row:
            result['mps1_pIC50']    = row['mps1_pIC50_mean']
            result['aurora_pIC50']  = row['aurora_pIC50_mean']
            result['herg_pIC50']    = row['herg_pIC50_mean']
            result['SI_multitask']  = row['selectivity_index']
            result['herg_risk']     = row['herg_risk']

        results.append(result)

        sel_flag = "✓ selective" if selective else "~ non-selective"
        print(f"  [{idx+1}/{args.top_n}] "
              f"{ligand_name[:35]}")
        print(f"    Mps1: {mps1_score:.3f} | "
              f"AurB: {aurora_score:.3f} | "
              f"Diff: {score_diff:+.3f} | {sel_flag}")

    except Exception as e:
        print(f"  [{idx+1}/{args.top_n}] "
              f"{ligand_name}: error — {e}")

results_df = pd.DataFrame(results)
results_df = results_df.sort_values(
    'score_diff', ascending=False
)
results_df.to_csv(
    RESULTS_DIR / "selectivity_results.csv", index=False
)

print(f"\n{'=' * 55}")
print(f"SELECTIVITY SCREEN COMPLETE")
print(f"{'=' * 55}")
print(f"  Screened:            {len(results_df)}")
print(f"  Selective (Δ>2.0):  "
      f"{results_df['selective'].sum()}/{len(results_df)}")

print(f"\nRanked by selectivity:")
print(f"{'Ligand':<35} {'Mps1':>7} {'AurB':>7} "
      f"{'Diff':>7} {'SI_MT':>8} {'hERG':>8}")
print("-" * 75)

for _, r in results_df.iterrows():
    name = r['ligand'].replace('phase2b_', '')[:33]
    si   = f"{r.get('SI_multitask', 0):.0f}"
    herg = r.get('herg_risk', 'N/A')
    print(f"  {name:<33} "
          f"{r['mps1_score']:>7.3f} "
          f"{r['aurora_score']:>7.3f} "
          f"{r['score_diff']:>+7.3f} "
          f"{si:>8} "
          f"{herg:>8}")

print(f"\nResults → "
      f"{RESULTS_DIR}/selectivity_results.csv")
