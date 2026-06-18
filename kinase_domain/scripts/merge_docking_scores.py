"""
Merge Smina ensemble scores (all 354 candidates) with
GNINA CNN rescoring (top 50 candidates).

Final ranking uses:
  - CNN score where available (top 50)
  - Smina score for remaining candidates
"""

import pandas as pd
import numpy as np
from pathlib import Path

SMINA_CSV = "kinase_domain/docking/phase2b_results/docking_scores.csv"
CNN_CSV   = "kinase_domain/docking/phase2b_results/top50_gnina_cnn_scores.csv"
OUTPUT    = "kinase_domain/docking/phase2b_results/final_docking_scores.csv"

print("Merging Smina + GNINA CNN scores")
print("=" * 55)

smina_df = pd.read_csv(SMINA_CSV)
print(f"Smina scores:    {len(smina_df)} candidates")

if not Path(CNN_CSV).exists():
    print("CNN scores not found — using Smina only")
    smina_df['final_score'] = smina_df['best_score']
    smina_df['score_source'] = 'smina'
    smina_df.to_csv(OUTPUT, index=False)
else:
    cnn_df = pd.read_csv(CNN_CSV)
    print(f"CNN scores:      {len(cnn_df)} candidates")

    merged = smina_df.merge(
        cnn_df[['ligand', 'best_cnn', 'mean_cnn',
                'best_vina', 'consensus']],
        on='ligand', how='left'
    )

    # Use CNN score where available, Smina otherwise
    merged['final_score'] = np.where(
        merged['best_cnn'].notna(),
        merged['best_cnn'],
        merged['best_score']
    )
    merged['score_source'] = np.where(
        merged['best_cnn'].notna(),
        'gnina_cnn', 'smina'
    )

    merged = merged.sort_values(
        'final_score', ascending=True
    )
    merged.to_csv(OUTPUT, index=False)

    print(f"\nScore sources:")
    print(f"  GNINA CNN: "
          f"{(merged['score_source']=='gnina_cnn').sum()}")
    print(f"  Smina:     "
          f"{(merged['score_source']=='smina').sum()}")

    print(f"\nTop 10 by final score:")
    print(f"{'Ligand':<35} {'Score':>7} {'Source':>10}")
    print("-" * 55)
    for _, r in merged.head(10).iterrows():
        name = r['ligand'].replace('phase2b_', '')[:33]
        print(f"  {name:<33} "
              f"{r['final_score']:>7.3f} "
              f"{r['score_source']:>10}")

print(f"\nSaved → {OUTPUT}")
