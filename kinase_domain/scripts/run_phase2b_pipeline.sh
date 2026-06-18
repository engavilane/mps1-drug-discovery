#!/bin/bash
# Phase 2B post-docking pipeline

set -e
echo "Phase 2B Post-Docking Pipeline"
echo "================================"

conda run -n docking python \
    kinase_domain/scripts/merge_docking_scores.py
echo "Step 1 done: scores merged"

conda run -n docking python \
    kinase_domain/scripts/plip_analysis.py \
    --scores   kinase_domain/docking/phase2b_results/final_docking_scores.csv \
    --results  kinase_domain/docking/phase2b_results \
    --output   kinase_domain/analysis/phase2b/plip \
    --receptor kinase_domain/data/receptor/ensemble/5AP7_ensemble_receptor.pdbqt
echo "Step 2 done: PLIP analysis"

conda run -n docking python \
    kinase_domain/scripts/filter_hinge_binders.py \
    kinase_domain/analysis/phase2b/plip/plip_interactions.csv \
    kinase_domain/analysis/phase2b/plip/hinge_binders.csv
echo "Step 3a done: hinge binders filtered"

conda run -n docking python \
    kinase_domain/scripts/admet_filter.py \
    --candidates kinase_domain/analysis/phase2b/plip/hinge_binders.csv \
    --ligands    kinase_domain/data/phase2b/raw \
    --output     kinase_domain/analysis/phase2b
echo "Step 3b done: ADMET filtering"

conda run -n docking python \
    kinase_domain/scripts/gnn_predict.py \
    --candidates kinase_domain/analysis/phase2b/admet_candidates.csv \
    --ligands    kinase_domain/data/phase2b/raw \
    --model_dir  kinase_domain/analysis/gnn_model/optimised \
    --output     kinase_domain/analysis/phase2b
echo "Step 4 done: GNN predictions"

conda run -n docking python \
    kinase_domain/scripts/pareto_optimisation.py \
    --candidates kinase_domain/analysis/phase2b/admet_candidates.csv \
    --output     kinase_domain/analysis/phase2b
echo "Step 5 done: Pareto optimisation"

echo "================================"
echo "Phase 2B pipeline complete"
echo "Results: kinase_domain/analysis/phase2b/phase2b_pareto_front1.csv"
