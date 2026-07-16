# Mps1/TTK In Silico Drug Discovery Pipeline

> Ensemble docking, graph neural network QSAR, and first computational allosteric targeting of the Mps1 NTD TPR domain — MSc Bioinformatics internship, Oxford Brookes University, 2026.

**Author:** Enya Gavilan Esquitino  
**Supervisor:** Prof. Victor M. Bolanos-Garcia  
**Institution:** Oxford Brookes University  

---

## Biological Context

**Mps1/TTK** (Monopolar Spindle 1 / Threonine Tyrosine Kinase) is the master regulator of the Spindle Assembly Checkpoint (SAC), overexpressed in triple-negative breast cancer, glioblastoma, pancreatic ductal adenocarcinoma, and colorectal cancer. This pipeline targets two complementary sites:

- **Kinase domain** — ATP-binding hinge region (Gly605/Glu603 hydrogen bond pharmacophore)
- **NTD TPR domain** — Hec1/Ndc80 protein-protein interaction interface (novel allosteric strategy)

---

## Pipeline Overview

### Campaign 1 — Kinase Domain: Validation & Baseline (Phase 1)
1. 45 co-crystallised Mps1 inhibitors retrieved from RCSB/KLIFS
2. Ensemble receptor preparation (75 KLIFS structures → 7 diverse conformations)
3. AutoDock Vina docking (exhaustiveness=16)
4. PLIP interaction analysis (Gly605/Glu603 H-bond detection, geometric criteria)
5. ADME filtering (Lipinski Ro5, TPSA, PAINS, Brenk, SA score)
6. ML Model 1 — Ridge Regression (n=45, LOOCV, R²=0.706)
7. ML Model 2 — SVR (n=2,352 ChEMBL, 10-fold CV, R²=0.729)
8. ML Model 3 — D-MPNN GNN ensemble (n=3,607, R²=0.816, σ=0.186)

### Campaign 2 — Kinase Domain: Novel Candidate Discovery (Phase 2B)
1. 5 diverse ChEMBL seeds (pIC50 ≥ 10.4, max pairwise Tanimoto ≤ 0.4)
2. PubChem similarity search (Tanimoto ≥ 0.7) → 395 candidates
3. Smina ensemble docking (3 receptors) + GNINA CNN rescoring (top 50)
4. PLIP hinge interaction analysis → 233/354 hinge binders
5. Full ADMET filtering (SwissADME) → 163/233 pass
6. D-MPNN GNN pIC50 prediction + uncertainty quantification
7. Multi-task GNN (Mps1 + Aurora B + hERG) → selectivity and safety
8. Aurora B counter-screen (docking into 4TND)
9. Multi-objective Pareto optimisation (5 objectives) → 25 front-1 candidates
10. Resistance mutation docking (C604Y mutant, 5NTT)

### Campaign 3 — NTD Domain: First Allosteric Campaign
1. 4H7Y TPR domain (1.8 Å, apo) receptor preparation
2. fpocket binding site detection → Hec1 interface (Pocket 42)
3. AlphaFold2-Multimer complex prediction (TPR pLDDT=96.1)
4. 988 ZINC fragments screened (MW ≤ 320 Da, LogP ≤ 3.5)
5. GNINA docking (exhaustiveness=16) → 964 scored
6. PLIP interface analysis (residues 104,107,108,111,140,141,144,145,173,174,177,178)
7. Top 10 interface binders identified (Lys107/Arg108 contacts)

---

## Key Results

| Metric | Value |
|--------|-------|
| Self-docking RMSD (5LJJ) | 0.666 Å |
| Cross-docking mean RMSD | 1.903 ± 0.633 Å |
| Homologous docking ρ | 0.833 (p=0.005) |
| PLIP Phase 1 hinge binders | 39/45 (86.7%) |
| PLIP Phase 2B hinge binders | 233/354 (65.8%) |
| Ridge Regression R² | 0.706 (LOOCV) |
| SVR R² | 0.729 (10-fold CV) |
| GNN R² | 0.816 (held-out test) |
| Y-randomisation Z-score | 15.88 (p < 0.01) |
| Phase 2B ADMET pass rate | 163/233 (70.0%) |
| Pareto front 1 candidates | 25 |
| **Top candidate pIC50** | **9.228 (IC50 = 0.59 nM)** |
| **Top candidate SI vs Aurora B** | **1353** |
| NTD top fragment score | -8.943 kcal/mol (zinc_div_0277) |
| NTD interface binders | 10/20 (Lys107 or Arg108) |
| C604Y resistance Δ | +2.56 to +2.70 kcal/mol (all susceptible) |

---

## Repository Structure

```
mps1-drug-discovery/
│
├── kinase_domain/
│   ├── data/
│   │   ├── raw/                        # PDB structures (5LJJ, 5NTT, etc.)
│   │   ├── receptor/
│   │   │   └── ensemble/               # 7 KLIFS-validated receptor PDBQTs
│   │   ├── ligands/
│   │   │   ├── raw/                    # 45 Phase 1 SDF files
│   │   │   └── pdbqt/                  # 45 Phase 1 PDBQT files
│   │   ├── phase2/                     # Phase 2A candidates (PDB-seeded)
│   │   ├── phase2b/                    # Phase 2B candidates (ChEMBL-seeded)
│   │   │   ├── raw/                    # 395 SDF files
│   │   │   └── pdbqt/                  # 385 PDBQT files
│   │   ├── libraries/                  # Chemical libraries
│   │   └── bindingdb_ttk.tsv           # BindingDB TTK activity data
│   │
│   ├── docking/
│   │   ├── results/                    # Phase 1 docked poses + scores
│   │   ├── phase2_results/             # Phase 2A docked poses + scores
│   │   ├── phase2b_results/            # Phase 2B docked poses + scores
│   │   ├── selectivity_results/        # Aurora B counter-screen (4TND)
│   │   └── resistance_results/         # C604Y mutant docking (5NTT)
│   │
│   ├── analysis/
│   │   ├── validation/                 # RMSD, cross-docking, homologous
│   │   ├── interactions/               # Legacy distance-based (deprecated)
│   │   ├── interactions_plip/
│   │   │   ├── phase1/                 # PLIP results (39/45 hinge binders)
│   │   │   └── phase2/                 # PLIP results (Phase 2A)
│   │   ├── adme/                       # Phase 1 ADME results
│   │   ├── adme_plip/                  # ADME + PLIP combined (Phase 1+2A)
│   │   ├── admet/                      # Full ADMET results
│   │   ├── ic50/                       # ChEMBL + BindingDB training data
│   │   ├── ml_model/
│   │   │   ├── best_model.pkl          # Ridge Regression (Model 1)
│   │   │   ├── ridge_coefficients.csv
│   │   │   └── chembl/                 # SVR + validation (Model 2)
│   │   ├── gnn_model/
│   │   │   ├── optimised/              # D-MPNN GNN ensemble (5 models)
│   │   │   └── multitask/              # Multi-task GNN (Mps1+AuroraB+hERG)
│   │   ├── phase2/                     # Phase 2A results
│   │   ├── phase2b/
│   │   │   ├── plip/                   # PLIP hinge binders (233/354)
│   │   │   ├── admet_candidates_clean.csv  # 163 ADMET-passing candidates
│   │   │   ├── phase2b_multitask_predictions.csv
│   │   │   └── swissadme_top5.csv
│   │   ├── pareto/                     # Pareto optimisation results
│   │   ├── selectivity/                # Selectivity analysis
│   │   ├── resistance/                 # Resistance mutation analysis
│   │   ├── consensus_docking/          # Consensus docking results
│   │   └── md_simulations/             # MD analysis results
│   │
│   ├── notebooks/
│   │   ├── Mps1_ChemProp_Training.ipynb   # D-MPNN hyperparameter optimisation
│   │   ├── GNINA_top50.ipynb              # Phase 2B GNINA CNN rescoring
│   │   └── Multi_task_GNN.ipynb           # Multi-task GNN training
│   │
│   └── scripts/
│       ├── download_ligands.py
│       ├── cleanup_ligands.py
│       ├── prepare_ligands.py
│       ├── ensemble_receptor_prep.py
│       ├── run_docking.py
│       ├── run_ensemble_docking.py
│       ├── run_smina_ensemble.py
│       ├── plip_analysis.py
│       ├── filter_hinge_binders.py
│       ├── admet_filter.py
│       ├── ml_model.py
│       ├── ml_chembl.py
│       ├── expand_training_data.py
│       ├── gnn_predict.py
│       ├── multitask_predict.py
│       ├── get_multitask_data.py
│       ├── prepare_multitask_data.py
│       ├── y_randomisation.py
│       ├── leverage_ad.py
│       ├── applicability_domain.py
│       ├── synthetic_accessibility.py
│       ├── phase2_similarity_search.py
│       ├── phase2b_seed_selection.py
│       ├── phase2_predict.py
│       ├── merge_docking_scores.py
│       ├── pareto_optimisation.py
│       ├── selectivity_screen.py
│       ├── rmsd_validation.py
│       ├── crossdock_5N7V.py
│       ├── crossdock_validation.py
│       ├── homologous_docking.py
│       ├── prepare_crossdock_receptor.py
│       ├── get_ic50.py
│       └── run_phase2b_pipeline.sh
│
├── ntd_domain/
│   ├── data/
│   │   ├── raw/                        # 4H7Y PDB + fixed structure
│   │   ├── receptor/                   # 4H7Y PDBQT receptor
│   │   └── libraries/
│   │       ├── zinc_1000_diverse.smi   # 988 ZINC fragment library
│   │       └── pdbqt/                  # Prepared fragment PDBQTs
│   ├── docking/
│   │   ├── results/                    # All GNINA docking results
│   │   └── top20_results/              # Top 20 docked poses
│   ├── analysis/
│   │   ├── ntd_docking_scores.csv      # All 964 fragment scores
│   │   ├── ntd_top_candidates.csv      # Top 10 fragments
│   │   ├── alphafold_model.png         # AlphaFold2-Multimer figure
│   │   ├── admet/                      # NTD fragment ADMET results
│   │   ├── binding_site/               # fpocket binding site analysis
│   │   ├── fragments/                  # Fragment analysis
│   │   ├── leads/                      # Lead fragment results
│   │   └── plip_top20/
│   │       └── ntd_plip_interactions.csv
│   └── scripts/
│       ├── ntd_receptor_prep.py
│       ├── prepare_zinc_library.py
│       ├── ntd_docking.py
│       └── ntd_plip_analysis.py
│
├── md_simulations/
│   ├── kinase_domain/
│   │   ├── systems/                    # 5 holo systems (EM-minimised)
│   │   ├── topol/                      # GROMACS topologies (GAFF2/AMBER99sb-ildn)
│   │   ├── mdp/                        # MDP files (EM, NVT, NPT, MD)
│   │   ├── analysis/                   # MD trajectory analysis
│   │   └── scripts/                    # MD preparation scripts
│   └── ntd_domain/
│       ├── systems/                    # NTD apo system
│       ├── topol/                      # NTD topology
│       ├── mdp/                        # NTD MDP files
│       ├── results/                    # NTD MD results (future work)
│       └── analysis/                   # NTD trajectory analysis
│
├── shared/
│   ├── figures/                        # Shared publication figures
│   └── models/                         # Shared ML model files
│
├── report_figures/
│   ├── poses/                          # PDB files for PyMOL visualisation
│   ├── kinase_top3_structures.png      # Chemical structures (top 3 kinase)
│   ├── ntd_top3_structures.png         # Chemical structures (top 3 NTD)
│   ├── kinase_docking_pose.png         # PyMOL — kinase binding pose
│   ├── ntd_docking_pose.png            # PyMOL — NTD binding pose
│   ├── resistance_WT_pose.png          # PyMOL — WT binding pose
│   ├── resistance_C604Y_pose.png       # PyMOL — C604Y mutant pose
│   ├── resistance_superimposed.png     # PyMOL — WT vs C604Y overlay
│   └── pkcsm_input.smi                 # SMILES for ADMET submission
│
├── README.md
└── environment.yml                     # Conda reproducible environment
```

---

## Reproducing the Pipeline

### Setup

```bash
git clone https://github.com/engavilane/mps1-drug-discovery
cd mps1-drug-discovery
conda env create -f environment.yml
conda activate docking
```

### Campaign 1 — Kinase Domain Phase 1

```bash
python kinase_domain/scripts/download_ligands.py
python kinase_domain/scripts/cleanup_ligands.py
python kinase_domain/scripts/prepare_ligands.py
python kinase_domain/scripts/ensemble_receptor_prep.py
python kinase_domain/scripts/run_docking.py
python kinase_domain/scripts/rmsd_validation.py
python kinase_domain/scripts/crossdock_validation.py
python kinase_domain/scripts/homologous_docking.py
python kinase_domain/scripts/plip_analysis.py
python kinase_domain/scripts/admet_filter.py
python kinase_domain/scripts/ml_model.py
python kinase_domain/scripts/ml_chembl.py
python kinase_domain/scripts/y_randomisation.py
python kinase_domain/scripts/leverage_ad.py
# GNN training: run Mps1_ChemProp_Training.ipynb on Google Colab
```

### Campaign 2 — Kinase Domain Phase 2B

```bash
python kinase_domain/scripts/phase2b_seed_selection.py
python kinase_domain/scripts/phase2_similarity_search.py
python kinase_domain/scripts/prepare_ligands.py \
    --raw kinase_domain/data/phase2b/raw \
    --pdbqt kinase_domain/data/phase2b/pdbqt
# GNINA docking: run GNINA_top50.ipynb on Colab GPU
bash kinase_domain/scripts/run_phase2b_pipeline.sh
# Multi-task GNN: run Multi_task_GNN.ipynb on Colab
python kinase_domain/scripts/selectivity_screen.py
```

### Campaign 3 — NTD Domain

```bash
python ntd_domain/scripts/ntd_receptor_prep.py
python ntd_domain/scripts/prepare_zinc_library.py
# NTD docking: run on Colab GPU (see ntd_domain/scripts/ntd_docking.py)
python ntd_domain/scripts/ntd_plip_analysis.py
```

---

## Top Candidates

### Kinase Domain (Phase 2B)

| Candidate | pIC50 | IC50 (nM) | SI vs Aurora B | hERG risk |
|-----------|-------|-----------|----------------|-----------|
| phase2b_seed_4_44199470 ★ | 9.228 | 0.59 | 1353 | moderate |
| phase2b_seed_4_66555847 | 8.870 | 1.35 | 270 | high |
| phase2b_seed_4_90072314 | 8.808 | 1.56 | 1117 | moderate |
| phase2b_seed_4_86713879 | 8.587 | 2.59 | 290 | moderate |
| phase2b_seed_4_44217530 | 8.431 | 3.71 | 150 | moderate |

★ Best overall — sub-nM potency, SI=1353, moderate safety profile.

### NTD Domain

| Fragment | Score (kcal/mol) | Contact | MW (Da) | LogP |
|----------|-----------------|---------|---------|------|
| zinc_div_0277 ★ | -8.943 | Arg108 | 318.3 | -1.24 |
| zinc_div_0276 | -8.862 | Lys107 | 318.3 | -1.39 |
| zinc_div_0300 | -8.422 | Lys107 | 318.4 | -1.06 |

★ Unique Arg108 contact — distinct binding mode within TPR groove.

---

## Software

| Tool | Version | Purpose |
|------|---------|---------|
| AutoDock Vina | 1.2.5 | Phase 1 docking |
| Smina | latest | Phase 2B ensemble docking |
| GNINA | 1.0.3 | CNN rescoring + NTD docking |
| Meeko | latest | Receptor/ligand preparation |
| PLIP | 2.x | Protein-ligand interaction analysis |
| Chemprop | 2.x | D-MPNN GNN (Models 3 + multi-task) |
| RDKit | 2024.x | Cheminformatics + ADMET filtering |
| scikit-learn | 1.x | Ridge + SVR models |
| GROMACS | 2026.0 | MD simulations (future work) |
| PyTorch | 2.x | GNN training |
| AlphaFold2-Multimer | v3 | NTD-Hec1 complex prediction |
| SwissADME | web | Extended ADMET panel |
| PyMOL | 2.x | Docking pose visualisation |

---

## References

Biology & Target

1. Musacchio A, Salmon ED. The spindle-assembly checkpoint in space and time. Nat Rev Mol Cell Biol. 2007. DOI: 10.1038/nrm2163
2. Pugh KM et al. In silico identification of potential Mps1 inhibitors from natural isoflavones. Int J Mol Sci. 2022. DOI: 10.3390/ijms232214228
3. Nijenhuis W et al. A TPR domain–containing N-terminal module of MPS1 is required for its kinetochore localization by Aurora B. J Cell Biol. 2013. DOI: 10.1083/jcb.201210117
4. Bolanos-Garcia VM. Trends Biochem Sci. 2025. (in press — ask supervisor for DOI)

Docking & Scoring

5. Trott O, Olson AJ. AutoDock Vina: improving the speed and accuracy of docking. J Comput Chem. 2010. DOI: 10.1002/jcc.21334
6. McNutt AT et al. GNINA 1.0: molecular docking with deep learning. J Cheminform. 2021. DOI: 10.1186/s13321-021-00522-2
7. Quiroga R, Villarreal MA. Vinardo: A Scoring Function Based on Autodock Vina Improves Scoring, Docking, and Virtual Screening. PLoS ONE. 2016. DOI: 10.1371/journal.pone.0155183

Structure Databases

8. Berman HM et al. The Protein Data Bank. Nucleic Acids Res. 2000. DOI: 10.1093/nar/28.1.235
9. Kanev GK et al. KLIFS: an overhaul of the kinase-ligand interaction fingerprints and structures database. Nucleic Acids Res. 2021. DOI: 10.1093/nar/gkaa895

Interaction Analysis

10. Salentin S et al. PLIP: fully automated protein-ligand interaction profiler. Nucleic Acids Res. 2015. DOI: 10.1093/nar/gkv315

Machine Learning & QSAR

11. Heid E et al. Chemprop: A machine learning package for chemical property prediction. J Chem Inf Model. 2024. DOI: 10.1021/acs.jcim.3c01250
12. Yang K et al. Analyzing learned molecular representations for property prediction. J Chem Inf Model. 2019. DOI: 10.1021/acs.jcim.9b00237
13. Liaw A, Wiener M. Classification and Regression by randomForest. R News. 2002. (scikit-learn implementation)
14. Cortes C, Vapnik V. Support-vector networks. Mach Learn. 1995. DOI: 10.1007/BF00994018

Cheminformatics & ADMET

15. RDKit: Open-source cheminformatics. rdkit.org
16.Ertl P, Schuffenhauer A. Estimation of synthetic accessibility score. J Cheminform. 2009. DOI: 10.1186/1758-2946-1-8
17. Daina A et al. SwissADME: a free web tool to evaluate pharmacokinetics and drug-likeness. Sci Rep. 2017. DOI: 10.1038/srep42717
18. Baell JB, Holloway GA. New substructure filters for removal of pan assay interference compounds (PAINS). J Med Chem. 2010. DOI: 10.1021/jm901137j

Activity Data

19. Gaulton A et al. ChEMBL: a large-scale bioactivity database for drug discovery. Nucleic Acids Res. 2012. DOI: 10.1093/nar/gkr777
20.Gilson MK et al. BindingDB in 2015: A public database for medicinal chemistry. Nucleic Acids Res. 2016. DOI: 10.1093/nar/gkv1072
21. Irwin JJ et al. ZINC20 — A free ultralarge-scale chemical database for ligand discovery. J Chem Inf Model. 2020. DOI: 10.1021/acs.jcim.0c00675

Binding Site & Structure Prediction

22. Le Guilloux V et al. Fpocket: An open source platform for ligand pocket detection. BMC Bioinformatics. 2009. DOI: 10.1186/1471-2105-10-168
23. Mirdita M et al. ColabFold: making protein folding accessible to all. Nat Methods. 2022. DOI: 10.1038/s41592-022-01488-1
24. Evans R et al. Protein complex prediction with AlphaFold-Multimer. bioRxiv. 2022. DOI: 10.1101/2021.10.04.463034

Ligand Preparation

25. Eberhardt J et al. AutoDock Vina 1.2.0: New docking methods, expanded force field, and Python bindings. J Chem Inf Model. 2021. DOI: 10.1021/acs.jcim.1c00203

MD Simulations

26. Abraham MJ et al. GROMACS: High performance molecular simulations through multi-level parallelism from laptops to supercomputers. SoftwareX. 2015. DOI: 10.1016/j.softx.2015.06.001
27. Wang J et al. Development and testing of a general amber force field. J Comput Chem. 2004. DOI: 10.1002/jcc.20035
28. Lindorff-Larsen K et al. Improved side-chain torsion potentials for the Amber ff99SB protein force field. Proteins. 2010. DOI: 10.1002/prot.22711

Visualisation

29. The PyMOL Molecular Graphics System, Version 2.0, Schrödinger LLC.

Multi-objective Optimisation

30. Deb K et al. A fast and elitist multiobjective genetic algorithm: NSGA-II. IEEE Trans Evol Comput. 2002. DOI: 10.1109/4235.996017


---

## Licence
This project is licensed under the MIT Licence — see [LICENSE](LICENSE) for details.


## AI Disclosure

This project used Claude (Anthropic) as an assistance tool for code development and results analysis. The author has reviewed, edited, and assumes full responsibility for all content.
