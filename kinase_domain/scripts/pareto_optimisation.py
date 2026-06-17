"""
Pareto optimisation for Phase 2B final candidate ranking.

Replaces the combined score (arbitrary equal weighting) with
multi-objective Pareto optimisation across:
  - Ensemble docking score (minimise — more negative = better)
  - GNN predicted pIC50 (maximise)
  - GNN uncertainty (minimise — lower = more reliable)
  - SA score (minimise — lower = easier to synthesise)
  - TPSA (minimise — lower = better absorption)

A compound is Pareto-optimal if no other compound is better
on ALL objectives simultaneously.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import argparse
import warnings
warnings.filterwarnings('ignore')

parser = argparse.ArgumentParser(
    description="Pareto optimisation for candidate ranking"
)
parser.add_argument("--candidates",
    default="kinase_domain/analysis/phase2b/admet_candidates.csv")
parser.add_argument("--output",
    default="kinase_domain/analysis/phase2b")
args = parser.parse_args()

OUTPUT_DIR = Path(args.output)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

print("Pareto Optimisation — Phase 2B Candidates")
print("=" * 55)

df = pd.read_csv(args.candidates)
print(f"Candidates: {len(df)}")

# Define objectives
# All converted so that LOWER = BETTER (minimisation)
objectives = {}

if 'best_score' in df.columns:
    objectives['docking'] = df['best_score'].values
elif 'affinity_kcal' in df.columns:
    objectives['docking'] = df['affinity_kcal'].values

if 'pIC50_gnn' in df.columns:
    objectives['potency'] = -df['pIC50_gnn'].values
elif 'pIC50_pred' in df.columns:
    objectives['potency'] = -df['pIC50_pred'].values

if 'pIC50_gnn_std' in df.columns:
    objectives['uncertainty'] = df['pIC50_gnn_std'].values

if 'SA_score' in df.columns:
    objectives['sa_score'] = df['SA_score'].values

if 'TPSA' in df.columns:
    objectives['tpsa'] = df['TPSA'].values

print(f"Objectives: {list(objectives.keys())}")

# Build objective matrix
obj_matrix = np.column_stack(
    [objectives[k] for k in objectives]
)
print(f"Objective matrix: {obj_matrix.shape}\n")

# Pareto dominance check
def is_dominated(a, b):
    """Return True if solution a is dominated by solution b.
    b dominates a if b is <= a on all objectives
    and strictly < on at least one."""
    return (
        np.all(b <= a) and
        np.any(b < a)
    )

def get_pareto_front(matrix):
    """Get indices of Pareto-optimal solutions."""
    n = len(matrix)
    pareto_mask = np.ones(n, dtype=bool)

    for i in range(n):
        if not pareto_mask[i]:
            continue
        for j in range(n):
            if i == j or not pareto_mask[j]:
                continue
            if is_dominated(matrix[i], matrix[j]):
                pareto_mask[i] = False
                break

    return np.where(pareto_mask)[0]

# Get Pareto fronts iteratively
print("Computing Pareto fronts...")
remaining_idx = np.arange(len(df))
all_fronts    = []
front_num     = 1

while len(remaining_idx) > 0:
    sub_matrix = obj_matrix[remaining_idx]
    pareto_idx = get_pareto_front(sub_matrix)
    front_idx  = remaining_idx[pareto_idx]

    all_fronts.append(front_idx)
    print(f"  Front {front_num}: {len(front_idx)} candidates")

    remaining_idx = np.setdiff1d(remaining_idx, front_idx)
    front_num += 1

    if front_num > 10:
        print(f"  ... {len(remaining_idx)} remaining candidates "
              f"in lower fronts")
        break

# Assign front numbers
df['pareto_front'] = 0
for i, front_idx in enumerate(all_fronts, 1):
    df.loc[front_idx, 'pareto_front'] = i

# Sort by Pareto front then by docking score
sort_col = (
    'best_score' if 'best_score' in df.columns
    else 'affinity_kcal'
)
df = df.sort_values(
    ['pareto_front', sort_col],
    ascending=[True, True]
)

pareto_front1 = df[df['pareto_front'] == 1]

print(f"\n{'=' * 55}")
print(f"PARETO OPTIMISATION RESULTS")
print(f"{'=' * 55}")
print(f"  Total candidates:     {len(df)}")
print(f"  Pareto front 1:       {len(pareto_front1)}")
print(f"  Number of fronts:     {len(all_fronts)}")

print(f"\nPareto Front 1 (optimal candidates):")
show_cols = (
    ['ligand', 'pareto_front', sort_col] +
    (['pIC50_gnn', 'pIC50_gnn_std']
     if 'pIC50_gnn' in df.columns
     else ['pIC50_pred']) +
    (['SA_score'] if 'SA_score' in df.columns else [])
)
show_cols = [c for c in show_cols if c in df.columns]

print(f"{'Ligand':<35} {'Front':>6} "
      f"{'Docking':>8} {'pIC50':>6} "
      f"{'±':>5} {'SA':>5}")
print("-" * 70)

for _, r in pareto_front1.head(15).iterrows():
    name  = r['ligand'].replace('phase2b_', '')[:33]
    score = r.get(sort_col, 0)
    pic50 = r.get('pIC50_gnn', r.get('pIC50_pred', 0))
    std   = r.get('pIC50_gnn_std', 0)
    sa    = r.get('SA_score', 0)
    print(f"  {name:<33} "
          f"{r['pareto_front']:>6} "
          f"{score:>8.3f} "
          f"{pic50:>6.3f} "
          f"{std:>5.3f} "
          f"{sa:>5.2f}")

df.to_csv(
    OUTPUT_DIR / "phase2b_pareto_candidates.csv",
    index=False
)
pareto_front1.to_csv(
    OUTPUT_DIR / "phase2b_pareto_front1.csv",
    index=False
)

print(f"\nAll candidates → "
      f"{OUTPUT_DIR}/phase2b_pareto_candidates.csv")
print(f"Pareto front 1 → "
      f"{OUTPUT_DIR}/phase2b_pareto_front1.csv")
