#!/bin/bash
# Phase 2B post-docking pipeline
# Run after GNINA CNN scores are downloaded from Colab

echo "Phase 2B Post-Docking Pipeline"
echo "================================"

# Step 1 — Merge Smina + CNN scores
echo "Step 1: Merging docking scores..."
conda run -n docking python \
    kinase_domain/scripts/merge_docking_scores.py

# Step 2 — PLIP interaction analysis
echo "Step 2: PLIP analysis..."
conda run -n docking python \
    kinase_domain/scripts/plip_analysis.py \
    --scores   kinase_domain/docking/phase2b_results/final_docking_scores.csv \
    --results  kinase_domain/docking/phase2b_results \
    --output   kinase_domain/analysis/phase2b/plip \
    --receptor kinase_domain/data/receptor/ensemble/5AP7_ensemble_receptor.pdbqt

# Step 3a - Filter hinge binders
echo "Step 3a: Filtering hinge binders..."
conda run -n docking python -c "
import pandas as pd
df = pd.read_csv (
    'kinase_domain/analysis/phase2b/plip/plip_interactions.csv'
)
hinge = df[df['is_hinge_binder'] == True]
print(f'Hinge binders: {len(hinge)}/{len(df)}')
hinge.to_csv(
    'kinase_domain/analysis/phase2b/plip/hinge_binders.csv',
    index = False
)

# Step 3b — ADMET filtering
echo "Step 3b: ADMET filtering..."
conda run -n docking python \
    kinase_domain/scripts/admet_filter.py \
    --candidates kinase_domain/analysis/phase2b/plip/hinge_binders.csv \
    --ligands    kinase_domain/data/phase2b/raw \
    --output     kinase_domain/analysis/phase2b

# Step 4 — GNN predictions
echo "Step 4: GNN pIC50 predictions..."
conda run -n docking python \
    kinase_domain/scripts/gnn_predict.py \
    --candidates kinase_domain/analysis/phase2b/admet_candidates.csv \
    --ligands    kinase_domain/data/phase2b/raw \
    --model_dir  kinase_domain/analysis/gnn_model/optimised \
    --output     kinase_domain/analysis/phase2b

# Step 5 — Pareto optimisation
echo "Step 5: Pareto optimisation..."
conda run -n docking python \
    kinase_domain/scripts/pareto_optimisation.py \
    --candidates kinase_domain/analysis/phase2b/admet_candidates.csv \
    --output     kinase_domain/analysis/phase2b

echo "================================"
echo "Phase 2B pipeline complete"
echo "Results: kinase_domain/analysis/phase2b/phase2b_pareto_front1.csv"
