# In Silico Drug Discovery Pipeline Targeting Mps1/TTK Kinase (PDB: 5LJJ)

> Structure-based virtual screening, ADME filtering, and machine learning-based affinity prediction for Mps1 kinase inhibitors — 1st year Master's internship project in Bioinformatics.

---

## About Me

I am a 1st year Master's student in Bioinformatics currently completing an internship in computational drug discovery. This project represents my first full end-to-end in silico drug discovery pipeline, combining structural bioinformatics, cheminformatics, and machine learning.

---

## Biological Context

The protein kinase **Mps1** (monopolar spindle 1), also known as **TTK**, is the most upstream regulator of the **Spindle Assembly Checkpoint (SAC)** — the evolutionarily conserved surveillance mechanism that ensures accurate chromosome segregation during cell division. Mps1 is overexpressed in a wide range of aggressive tumour types including breast, bladder, thyroid, pancreatic, and brain cancers, making it a promising drug target for cancer therapy.

Several Mps1 inhibitors have reached clinical trials, including **BAY-1217389**, **BAY-1161909**, **BOS-172722**, and **CFI-402257**, though none have yet received clinical approval. All known small-molecule inhibitors are ATP-competitive and bind the kinase hinge region, establishing key hydrogen bond interactions with residues **Gly605** and **Glu603**.

This project uses the crystal structure of Mps1 in complex with reversine (PDB: **5LJJ**, resolution 3.00 Å) as the receptor for molecular docking.

---

## Project Overview

This project is structured in two phases:

### Phase 1 — Virtual Screening of Known Inhibitors
1. Download and prepare 45 co-crystallised Mps1 inhibitors (RCSB PDB)
2. Receptor preparation (Meeko/AutoDock Vina)
3. Molecular docking (AutoDock Vina, exhaustiveness=16)
4. **PLIP interaction analysis** — proper H-bond detection
   (replaces distance-only contact counting)
5. ADME filtering (RDKit, Lipinski Ro5 + TPSA + PAINS)
6. Two complementary ML models:
   - **Model 1** — Interpretable Ridge Regression identifying
     physicochemical drivers of hinge binding (n=45)
   - **Model 2** — SVR predicting experimental pIC50 from ChEMBL data (n=2,352, R²=0.729)

### Phase 2 — Novel Candidate Discovery *(complete)*
The best inhibitors identified in Phase 1 are used as seeds for:
1. PubChem structural similarity search (Tanimoto index ≥ 0.9,
   Pugh et al. 2022 physicochemical filters)
2. Ligand preparation and molecular docking (AutoDock Vina,
   exhaustiveness=8, validated at exhaustiveness=16)
3. PLIP hinge interaction analysis (Gly605/Glu603 H-bond
   detection with distance + angle criteria)
4. ADME filtering (RDKit, Lipinski Ro5 + TPSA + PAINS)
5. Affinity prediction using the ChEMBL-trained SVR model
6. Combined ranking (Vina score + predicted pIC50)

---
## Pipeline Architecture

```
RCSB PDB (45 Mps1 inhibitor complexes)
         │
         ▼
Ligand Download (PubChem API / RCSB fallback)
         │
         ▼
Receptor & Ligand Preparation (Meeko)
         │
         ▼
Molecular Docking (AutoDock Vina 1.2)
    center=[-34.48, -15.66, -10.38]
    box_size=[20, 33, 21] Å
    exhaustiveness=16
         │
         ▼
Interaction Analysis (Gly605 / Glu603 H-bonds)
         │
         ▼
ADME Filtering (RDKit)
    Lipinski Ro5 + TPSA + PAINS
         │
         ├─────────────────────────────────────┐
         ▼                                     ▼
Model 1 — Vina scores (n=45)       Model 2 — Experimental (n=2352)
Ridge Regression                   SVR (PhysChem + Fingerprints)
PhysChem descriptors               ChEMBL Mps1 IC50 data
R²=0.706 (LOOCV)                   R²=0.729 (80/20 split)
         │                                     │
         └─────────────────┬───────────────────┘
                           ▼
              Phase 2: Similarity Search
              Novel candidate prediction
```

---

## Key Results

### Phase 1 — Docking & Interaction Analysis

Top 5 predicted binders by AutoDock Vina score:

| Ligand | PDB ID | Affinity (kcal/mol) | Gly605 H-bonds | Glu603 H-bonds | ADME |
|--------|--------|---------------------|----------------|----------------|------|
| 7CHN_LIG | 7CHN | -10.324 | 2 | 0 | ✓ |
| 7CHM_LIG | 7CHM | -10.288 | 4 | 1 | ✓ |
| 7CJA_LIG | 7CJA | -10.281 | 4 | 1 | ✓ |
| 7CHT_LIG | 7CHT | -9.970  | 2 | 1 | ✓ |
| 7CIL_LIG | 7CIL | -9.482  | — | — | ✓ |

Reversine (5LJJ_AD5) scored **-9.645 kcal/mol**, consistent with 
published Glide XP scores (-10.97 kcal/mol, Pugh et al. 2022), 
validating the docking setup. Clinically advanced compounds 
BAY-1217389 and BAY-1161909 scored -9.91 and -9.644 kcal/mol 
respectively, in line with their known potency. 39/45 inhibitors
confirmed genuine hinge binding by PLIP H-bond analysis 
(Gly605/Glu603, distance + angle criteria). Notably,
BAY-1161909 (5N9S_LIG, Gly605 : 1 H-bond) was correctly
reclassified as a hinge binder, previously missed by 
distance-only analysis. 28/45 passed both ADME filtering and 
hinge binding criteria. 

> Note: AutoDock Vina scores systematically overestimate binding 
> affinities relative to experimental values. Results are interpreted 
> in terms of relative ranking rather than absolute affinity prediction, 
> consistent with established practice in the field 
> (Pugh et al., 2022; Bolanos-Garcia, 2025).

**Docking Validation (three-tier hierarchy):**

| Test | Result | Status |
|------|--------|--------|
| Self-docking RMSD (5LJJ) | 0.666 Å | ✓ Exceptional |
| Cross-docking 5N7V (purine scaffold) | 1.332 Å | ✓ Pass |
| Cross-docking 4JS8 (indazole scaffold) | 1.590 Å | ✓ Pass |
| Cross-docking 5NAD (methylbenzamide) | 2.786 Å | ⚠ Borderline |
| Cross-docking 7LQD (covalent RMS-07) | 3.193 Å | ✗ Expected |
| Homologous docking Spearman ρ | 0.833 (p=0.005) | ✓ Validated |
| Mean score difference | 0.094 ± 0.524 kcal/mol | ✓ No bias |

All cross-docking structures aligned using DSSP rigid
secondary structure core (helices + sheets, excluding
activation loop 672-690 and hinge loop 603-605).
7LQD failure attributed to irreversible Cys604 covalent
bond — expected and scientifically appropriate.

---

### Machine Learning Models

**Model 1 — Interpretable ML model to identify physicochemical drivers of Mps1 hinge binding (Ridge Regression, n=45, LOOCV)**

| Model | Features | R² | RMSE | MAE |
|---|---|---|---|---|
| **Ridge Regression** | PhysChem only | **0.706** | 0.467 | 0.370 |
| Random Forest | PhysChem only | 0.471 | 0.626 | 0.462 |
| SVR | PhysChem only | 0.409 | 0.662 | 0.467 |

Ridge coefficients reveal key molecular drivers of hinge binding:
aromatic rings (π-stacking), TPSA (polar hinge contacts),
Gly605/Glu603 contacts (direct anchoring), and low HBA
(reduced desolvation penalty). Model used for mechanistic
interpretation, not affinity prediction.

**Model 2 — Experimental pIC50 prediction (n=2,352, ChEMBL, 80/20 split)**
| Model | Features | R² | RMSE | MAE |
|---|---|---|---|---|
| **SVR** | PhysChem + Fingerprints | **0.729** | 0.623 | 0.445 |
| Random Forest | PhysChem + Fingerprints | 0.701 | 0.653 | 0.492 |
| Ridge Regression | PhysChem + Fingerprints | 0.601 | 0.755 | 0.596 |

**10-fold cross validation (SVR, best model):**
| Metric | Mean | Std | 95% CI |
|---|---|---|---|
| R² | 0.702 | 0.049 | [0.605, 0.799] |
| RMSE | 0.640 | 0.043 | [0.555, 0.725] |
| MAE | 0.453 | 0.028 | [0.399, 0.507] |

Single split R²=0.729 vs CV R²=0.702 (ΔR²=0.027) — no overfitting confirmed.

Trained on 2,352 Mps1/TTK inhibitors from ChEMBL with experimental 
IC50 values. Hyperparameter tuning (GridSearchCV, 120 fits) confirmed 
near-optimal default parameters.
> The transition from Ridge Regression (n=45) to SVR (n=2,352) 
> directly illustrates the bias-variance tradeoff: complex models 
> require sufficient data to outperform simple linear baselines. 
> Morgan fingerprints were detrimental with n=45 but essential 
> with n=2,352, further illustrating this principle.

### Phase 2 — Novel Candidate Discovery

**235 novel candidates** retrieved from PubChem
(Tanimoto ≥ 0.9, 5 seed compounds).

| Stage | Count | Rate |
|-------|-------|------|
| Retrieved from PubChem | 235 | — |
| Successfully prepared | 232 | 99% |
| PLIP hinge binders | 216 | 93% |
| Pass ADME + hinge | 210 | 91% |

| Rank | CID | Seed | Affinity (kcal/mol) | Gly605 | Glu603 | IC50 pred |
|------|-----|------|---------------------|--------|--------|-----------|
| 1 | 155119000 | 7CHN | -10.860 | 5 | 1 | — |
| 2 | 155118997 | 7CHN | -10.735 | 1 | 1 | — |
| 3 | 155118910 | 7CHN | -10.439 | 3 | 0 | — |
| 4 | 155109147 | 7CHN | -10.341 | 3 | 1 | — |
| 5 | 155109167 | 7CHN | -10.333 | 4 | 1 | — |

**Top candidate by combined score: phase2_7CHM_LIG_146393701**
- Rank 1 under equal weighting (0.5/0.5) and Vina-heavy (0.7/0.3)
- Identified only by PLIP analysis — missed by distance-based
  interaction analysis, demonstrating the value of proper H-bond
  detection

**Priority candidate for experimental validation: CID 142416385**
- IUPAC: 4-[(1-methylcyclopropyl)amino]-2-[(5-methyl-1-
  propan-2-ylpyrazol-4-yl)amino]-7H-pyrrolo[2,3-d]
  pyrimidine-5-carbonitrile
- Formula: C₁₈H₂₂N₈ | MW: 350.4 Da
- Vina: −8.724 kcal/mol (confirmed exhaustiveness=16)
- Predicted pIC50: 7.831 → **IC50 = 14.74 nM**
- PLIP H-bonds: 6 (Gly605 + Glu603) ✓
- ADME: all filters passed ✓
- Scaffold: 7H-pyrrolo[2,3-d]pyrimidine
- Consistent across all weighting schemes (ranks 2-18) —
  robust to combined score weighting choice

**Combined score sensitivity analysis:**
| Weighting | Vina-heavy (0.7/0.3) | Equal (0.5/0.5) | pIC50-heavy (0.3/0.7) |
|---|---|---|---|
| Spearman ρ vs equal | 0.889 | — | 0.933 |

Equal weighting validated — top candidates robust across all schemes.

> Re-docking of top 5 candidates at exhaustiveness=16 yielded
> scores within 0.08 kcal/mol of screening values, confirming
> pose convergence and validating the exhaustiveness=8 screening
> protocol.

---

## Repository Structure

```
docking-5LJJ/
├── README.md                                # Project description and results
├── environment.yml                          # Conda reproducible environment
│
├── data/
│   ├── raw/
│   │   ├── 5LJJ.pdb                        # Original unmodified PDB from RCSB
│   │   ├── 5N7V.pdb                        # Cross-docking validation structure
│   │   ├── 4JS8.pdb                        # Cross-docking validation structure
│   │   ├── 5NAD.pdb                        # Cross-docking validation structure
│   │   └── 7LQD.pdb                        # Cross-docking validation structure
│   ├── receptor/
│   │   ├── receptor_clean.pdb              # Waters/ligands/artefacts removed
│   │   ├── receptor.pdbqt                  # Vina-ready receptor (Meeko)
│   │   ├── reversine_crystal_coords.npy    # Crystal coordinates for RMSD validation
│   │   ├── 5N7V_aligned_dssp.pdb           # 5N7V aligned to 5LJJ (DSSP rigid core)
│   │   ├── 4JS8_aligned.pdb                # 4JS8 aligned to 5LJJ
│   │   ├── 5NAD_aligned.pdb                # 5NAD aligned to 5LJJ
│   │   ├── 7LQD_aligned_trimmed.pdb        # 7LQD aligned to 5LJJ (trimmed)
│   │   └── *_homo_aligned.pdb              # Top 10 receptors for homologous docking
│   ├── ligands/
│   │   ├── compounds.csv                   # PDB ID input list
│   │   ├── download_log.csv                # Download status per ligand
│   │   ├── preparation_log.txt             # Meeko preparation log
│   │   ├── native_ligand.pdb               # Reversine (AD5) — reference ligand
│   │   ├── raw/                            # 45 SDF files (PubChem/RCSB)
│   │   └── pdbqt/                          # 45 prepared PDBQT files (Meeko)
│   └── phase2/
│       ├── raw/                            # 235 novel candidate SDF files
│       ├── pdbqt/                          # 232 prepared PDBQT files
│       └── preparation_log.txt             # Meeko preparation log (Phase 2)
│
├── docking/
│   ├── results/
│   │   ├── phase1_docking_scores.csv       # Vina scores for all 45 ligands
│   │   ├── *_out.pdbqt                     # Docked poses (one per ligand)
│   │   ├── 5LJJ_AD5_redock_out.pdbqt       # Reversine re-docked (RMSD validation)
│   │   └── reversine_crossdock_*_out.pdbqt # Cross-docking poses
│   └── phase2_results/
│       ├── docking_scores.csv              # Vina scores for 232 candidates
│       ├── *_out.pdbqt                     # Docked poses (one per candidate)
│       └── *_refined_out.pdbqt             # Top 5 re-docked at exhaustiveness=16
│
├── analysis/
│   ├── interactions/
│   │   ├── interaction_analysis.csv        # DEPRECATED — distance-based contacts
│   │   └── figures/
│   │       ├── 5LJJ-AD5.png               # Reversine in binding site
│   │       ├── 5LJJ-AD5_zoom.png          # Reversine — zoom on hinge contacts
│   │       ├── 5LJJ-7CHN.png              # Best Phase 1 binder (7CHN)
│   │       ├── 5LJJ-7CHN_zoom.png         # 7CHN — zoom on hinge contacts
│   │       ├── 5LJJ-5N9S.png              # BAY-1161909 binding mode
│   │       ├── 5LJJ-5N9S_zoom.png         # BAY-1161909 — zoom on hinge
│   │       ├── 5LJJ_reversine_pocket.png  # Reversine with pocket surface
│   │       ├── comparison_reversine_candidate.png  # Reversine vs CID 142416385
│   │       ├── phase2_top_candidate.png   # CID 142416385 in binding site
│   │       ├── phase2_top_candidate_pocket.png     # CID 142416385 with pocket
│   │       └── Figure_prep_PyMOL.ipynb    # PyMOL figure preparation notebook
│   ├── interactions_plip/
│   │   ├── phase1/
│   │   │   └── plip_interactions.csv      # PLIP H-bond analysis — 39/45 binders
│   │   └── phase2/
│   │       └── plip_interactions.csv      # PLIP H-bond analysis — 216/232 binders
│   ├── adme/
│   │   ├── adme_full.csv                   # All 45 ligands + descriptors + flags
│   │   └── adme_candidates.csv             # 26 candidates (distance-based, deprecated)
│   ├── adme_plip/
│   │   ├── phase1/
│   │   │   ├── adme_full.csv               # All 45 ligands + descriptors + flags
│   │   │   └── adme_candidates.csv         # 28 candidates passing ADME + hinge
│   │   └── phase2/
│   │       ├── adme_full.csv               # All 232 candidates + descriptors
│   │       └── adme_candidates.csv         # 210 candidates passing ADME + hinge
│   ├── ic50/
│   │   ├── chembl_mps1_ic50.csv            # 2,352 Mps1 IC50 values from ChEMBL
│   │   └── ic50_matches.csv                # Matches between our 45 and ChEMBL
│   ├── ml_model/
│   │   ├── best_model.pkl                  # Model 1 — Ridge Regression (Vina)
│   │   ├── scaler.pkl                      # Scaler for Model 1
│   │   ├── variance_selector.pkl           # Variance filter for Model 1
│   │   ├── predictions.csv                 # Model 1 predicted vs actual
│   │   ├── ridge_coefficients.csv          # Ridge feature coefficients
│   │   └── chembl/
│   │       ├── chembl_best_model.pkl       # Model 2 — SVR (ChEMBL pIC50)
│   │       ├── chembl_scaler.pkl           # Scaler for Model 2
│   │       ├── chembl_variance_selector.pkl
│   │       ├── chembl_predictions.csv      # Model 2 predicted vs actual
│   │       ├── chembl_model_comparison.csv
│   │       └── chembl_svr_tuned.pkl        # Hyperparameter-tuned SVR
│   ├── validation/
│   │   ├── rmsd_validation.csv             # Self-docking RMSD (0.666 Å)
│   │   ├── crossdocking_validation.csv     # Cross-docking across 4 conformations
│   │   ├── homologous_docking.csv          # Top 9 in native receptors (ρ=0.833)
│   │   └── top5_refined_scores.csv         # Top 5 Phase 2 re-docked (e=16)
│   └── phase2/
│       ├── phase2_candidates.csv           # 235 candidates from similarity search
│       ├── phase2_final_candidates.csv     # 210 ranked by combined score (PLIP)
│       ├── interactions/
│       │   └── interaction_analysis.csv    # DEPRECATED — distance-based contacts
│       └── adme/
│           ├── adme_full.csv               # All 232 candidates + descriptors
│           └── adme_candidates.csv         # 193 (distance-based, deprecated)
│
├── notebooks/
│   ├── 01_docking_analysis.ipynb           # Docking results + validation
│   ├── 02_adme_filter.ipynb                # ADME filtering visualisation
│   └── 03_ml_model.ipynb                   # ML model analysis and plots
│
└── scripts/
    ├── download_ligands.py                 # PubChem/RCSB automated download
    ├── cleanup_ligands.py                  # File renaming and deduplication
    ├── prepare_ligands.py                  # Meeko PDBQT preparation (Phase 1+2)
    ├── run_docking.py                      # AutoDock Vina batch docking (Phase 1+2)
    ├── parse_results.py                    # Score extraction and ranking
    ├── interaction_analysis.py             # DEPRECATED — replaced by plip_analysis.py
    ├── plip_analysis.py                    # PLIP H-bond detection (Phase 1+2)
    ├── adme_filter.py                      # RDKit ADME descriptors + PAINS (Phase 1+2)
    ├── ml_model.py                         # Model 1 — interpretable feature importance
    ├── get_ic50.py                         # ChEMBL IC50 retrieval
    ├── ml_chembl.py                        # Model 2 — pIC50 prediction (ChEMBL)
    ├── phase2_similarity_search.py         # PubChem similarity search (5 seeds)
    ├── phase2_predict.py                   # pIC50 prediction + combined ranking
    ├── rmsd_validation.py                  # Self-docking RMSD validation
    ├── crossdock_5N7V.py                   # Cross-docking 5N7V (DSSP alignment)
    ├── crossdock_validation.py             # Cross-docking across 4 conformations
    └── homologous_docking.py               # Top 10 in native receptors (Spearman ρ)
```

---

## Requirements

Create and activate the conda environment:

```bash
conda env create -f environment.yml
conda activate docking
```
Key dependencies: `vina`, `meeko`, `rdkit`, `pubchempy`, `pandas`, 
`scikit-learn`, `requests`, `numpy`, `scipy`, `joblib`, 
`chembl-webresource-client`, `mdanalysis`

---

## Reproducing the Pipeline

### Setup
```bash
conda env create -f environment.yml
conda activate docking
```

### Phase 1

#### 1. Download ligands
```bash
python scripts/download_ligands.py
```

#### 2. Clean and rename ligand files
```bash
python scripts/cleanup_ligands.py
```

#### 3. Prepare receptor
```bash
mk_prepare_receptor.py -i data/receptor/receptor_clean.pdb \
                       -o data/receptor/receptor -p
```

#### 4. Prepare ligands
```bash
python scripts/prepare_ligands.py
```

#### 5. Run docking
```bash
python scripts/run_docking.py
```

#### 6. PLIP interaction analysis (Gly605/Glu603)

Replaces distance-only contact counting with proper H-bond
detection (distance + angle criteria, PLIP algorithm).

```bash
python scripts/plip_analysis.py \
  --scores    docking/results/phase1_docking_scores.csv \
  --results   docking/results \
  --output    analysis/interactions_plip/phase1 \
  --old_interactions analysis/interactions/interaction_analysis.csv
```

> **Note:** `interaction_analysis.py` is retained for
> reference but deprecated. PLIP identified 39/45 hinge
> binders vs 42/45 by distance-only, correcting 9
> misclassifications including BAY-1161909.


#### 7. ADME filtering
```bash
python scripts/adme_filter.py
```

#### 8. ML Model 1 — Vina score prediction (n=45)
```bash
python scripts/ml_model.py
```

#### 9. Retrieve ChEMBL IC50 data
```bash
python scripts/get_ic50.py
```

#### 10. ML Model 2 — Experimental pIC50 prediction (n=2,352)
```bash
python scripts/ml_chembl.py
```

### Phase 2 — Novel Candidate Discovery

#### 11. Similarity search
```bash
python scripts/phase2_similarity_search.py
```

#### 12. Prepare Phase 2 ligands
```bash
python scripts/prepare_ligands.py \
  --raw   data/phase2/raw \
  --pdbqt data/phase2/pdbqt \
  --log   data/phase2/preparation_log.txt
```

#### 13. Dock Phase 2 candidates
```bash
python scripts/run_docking.py \
  --ligands        data/phase2/pdbqt \
  --results        docking/phase2_results \
  --exhaustiveness 8
```

#### 14. PLIP interaction analysis (Phase 2)

```bash
python scripts/plip_analysis.py \
  --scores    docking/phase2_results/docking_scores.csv \
  --results   docking/phase2_results \
  --output    analysis/interactions_plip/phase2 \
  --old_interactions analysis/phase2/interactions/interaction_analysis.csv
```

#### 15. ADME filtering
```bash
python scripts/adme_filter.py \
  --ligands      data/phase2/raw \
  --interactions analysis/phase2/interactions/interaction_analysis.csv \
  --scores       docking/phase2_results/docking_scores.csv \
  --output       analysis/phase2/adme
```

#### 16. Predict pIC50 and rank candidates
```bash
python scripts/phase2_predict.py
```

> Phase 2 docking uses exhaustiveness=8 for screening efficiency.
> Top candidates were validated at exhaustiveness=16 (score
> differences < 0.08 kcal/mol confirming convergence).

> Note: Steps 12–16 reuse the same scripts as Phase 1 via
> command-line arguments, ensuring methodological consistency
> across both phases. Phase 2 docking uses exhaustiveness=8
> (vs 16 in Phase 1) for computational efficiency during
> screening; top candidates can be re-docked at
> exhaustiveness=16 for validation.

---

## References

### Biological Context
- Bolanos-Garcia, V.M. (2025). Mps1 kinase functions in mitotic 
  spindle assembly and error correction. *Trends in Biochemical 
  Sciences*, 50(5), 438–453.
- Pugh, L. et al. (2022). Computational Biology Dynamics of Mps1 
  Kinase Molecular Interactions with Isoflavones Reveals a Chemical 
  Scaffold with Potential to Develop New Therapeutics for the 
  Treatment of Cancer. *Int. J. Mol. Sci.*, 23, 14228.
- Hiruma, Y. et al. (2016). Structural basis of reversine selectivity 
  in inhibiting Mps1 more potently than Aurora B kinase. *Proteins*, 
  84, 1761–1766. *(PDB: 5LJJ)*

### Molecular Docking
- Eberhardt, J. et al. (2021). AutoDock Vina 1.2.0: New Docking 
  Methods, Expanded Force Field, and Python Bindings. *J. Chem. Inf. 
  Model.*, 61, 3891–3898.
- Trott, O. & Olson, A.J. (2010). AutoDock Vina: improving the speed 
  and accuracy of docking. *J. Comput. Chem.*, 31, 455–461.
- Salentin S. et al. (2015) PLIP: fully automated 
  protein-ligand interaction profiler. *Nucleic Acids Res.*
  43:W443-W447. https://doi.org/10.1093/nar/gkv315

### Ligand & Receptor Preparation
- Forli, S. et al. (2016). Computational protein-ligand docking and 
  virtual drug screening with the AutoDock suite. *Nature Protocols*, 
  11, 905–919. *(Meeko/AutoDockTools)*

### Cheminformatics & ADME
- RDKit: Open-source cheminformatics software. 
  https://www.rdkit.org
- Daina, A. et al. (2017). SwissADME: A free web tool to evaluate 
  pharmacokinetics, drug-likeness and medicinal chemistry friendliness 
  of small molecules. *Sci. Rep.*, 7, 42717.
- Baell, J.B. & Holloway, G.A. (2010). New substructure filters for 
  removal of pan assay interference compounds (PAINS) from screening 
  libraries. *J. Med. Chem.*, 53, 2719–2740.
- Lipinski, C.A. et al. (2001). Experimental and computational 
  approaches to estimate solubility and permeability in drug discovery 
  and development settings. *Adv. Drug Deliv. Rev.*, 46, 3–26.

### Machine Learning
- Pedregosa, F. et al. (2011). Scikit-learn: Machine Learning in 
  Python. *JMLR*, 12, 2825–2830.
- Morgan, H.L. (1965). The generation of a unique machine description 
  for chemical structures. *J. Chem. Doc.*, 5, 107–113. 
  *(Morgan fingerprints)*

### Databases
- Berman, H.M. et al. (2000). The Protein Data Bank. 
  *Nucleic Acids Res.*, 28, 235–242.
- Mendez, D. et al. (2019). ChEMBL: towards direct deposition of 
  bioassay data. *Nucleic Acids Res.*, 47, D930–D940.
- Kim, S. et al. (2016). PubChem substance and compound databases. 
  *Nucleic Acids Res.*, 44, D1202–D1213.

---

## Author

**Enya Gavilan Esquitino** — Master's in Bioinformatics Internship, 2026
