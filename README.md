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
A dataset of **45 co-crystallised Mps1 inhibitors** (retrieved from 
the RCSB PDB) is subjected to:
1. Molecular docking against the Mps1 kinase domain (AutoDock Vina)
2. Binding interaction analysis (H-bond contacts with Gly605/Glu603)
3. ADME filtering using RDKit (Lipinski's Rule of Five, TPSA, PAINS)
4. Two complementary ML models:
   - **Model 1** — Ridge Regression on Vina scores (n=45, R²=0.706)
   - **Model 2** — SVR on experimental pIC50 from ChEMBL (n=2,352, R²=0.729)

### Phase 2 — Novel Candidate Discovery *(in progress)*
The best inhibitors identified in Phase 1 are used as seeds for:
1. PubChem structural similarity search (Tanimoto index ≥ 0.9)
2. Docking and ADME filtering of retrieved candidates
3. Hinge interaction analysis (Gly605/Glu603)
4. Affinity prediction using the ChEMBL-trained SVR model (Model 2)
5. Ranking and characterisation of novel predicted Mps1 inhibitors

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

| Ligand | PDB ID | Affinity (kcal/mol) | Gly605 | Glu603 | ADME |
|--------|--------|---------------------|--------|--------|------|
| 7CHN_LIG | 7CHN | -10.324 | 4 | 2 | ✓ |
| 7CHM_LIG | 7CHM | -10.288 | 4 | 2 | ✓ |
| 7CJA_LIG | 7CJA | -10.281 | 3 | 2 | ✓ |
| 7CHT_LIG | 7CHT | -9.970  | 4 | 2 | ✓ |
| 7CIL_LIG | 7CIL | -9.482  | 5 | 2 | ✓ |

Reversine (5LJJ_AD5) scored **-9.645 kcal/mol**, consistent with 
published Glide XP scores (-10.97 kcal/mol, Pugh et al. 2022), 
validating the docking setup. Clinically advanced compounds 
BAY-1217389 and BAY-1161909 scored -9.91 and -9.644 kcal/mol 
respectively, in line with their known potency. 42/45 inhibitors 
confirmed correct hinge binding (Gly605/Glu603 contacts). 
26/45 passed both ADME filtering and hinge binding criteria.

> Note: AutoDock Vina scores systematically overestimate binding 
> affinities relative to experimental values. Results are interpreted 
> in terms of relative ranking rather than absolute affinity prediction, 
> consistent with established practice in the field 
> (Pugh et al., 2022; Bolanos-Garcia, 2025).

---

### Phase 2 — Machine Learning Models

**Model 1 — Vina score prediction (n=45, LOOCV)**

| Model | Features | R² | RMSE | MAE |
|---|---|---|---|---|
| **Ridge Regression** | PhysChem only | **0.706** | 0.467 | 0.370 |
| Random Forest | PhysChem only | 0.471 | 0.626 | 0.462 |
| SVR | PhysChem only | 0.409 | 0.662 | 0.467 |

Key drivers of affinity: aromatic rings, LogP, TPSA, 
Gly605/Glu603 contacts (Ridge coefficients).

**Model 2 — Experimental pIC50 prediction (n=2,352, ChEMBL, 80/20 split)**

| Model | Features | R² | RMSE | MAE |
|---|---|---|---|---|
| **SVR** | PhysChem + Fingerprints | **0.729** | 0.623 | 0.445 |
| Random Forest | PhysChem + Fingerprints | 0.701 | 0.653 | 0.492 |
| Ridge Regression | PhysChem + Fingerprints | 0.601 | 0.755 | 0.596 |

Trained on 2,352 Mps1/TTK inhibitors from ChEMBL with experimental 
IC50 values. Hyperparameter tuning (GridSearchCV, 120 fits) confirmed 
near-optimal default parameters.

> The transition from Ridge Regression (n=45) to SVR (n=2,352) 
> directly illustrates the bias-variance tradeoff: complex models 
> require sufficient data to outperform simple linear baselines. 
> Morgan fingerprints were detrimental with n=45 but essential 
> with n=2,352, further illustrating this principle.

---

## Repository Structure

```
docking-5LJJ/
├── README.md                            # Project description and results
├── environment.yml                      # Conda reproducible environment
│
├── data/
│   ├── raw/
│   │   └── 5LJJ.pdb                    # Original unmodified PDB from RCSB
│   ├── receptor/
│   │   ├── receptor_clean.pdb          # Waters/ligands/artefacts removed
│   │   └── receptor.pdbqt              # Vina-ready receptor (Meeko)
│   └── ligands/
│       ├── compounds.csv               # PDB ID input list
│       ├── download_log.csv            # Download status per ligand
│       ├── preparation_log.txt         # Meeko preparation log
│       ├── native_ligand.pdb           # Reversine (AD5) — reference ligand
│       ├── raw/                        # 45 SDF files (PubChem/RCSB)
│       └── pdbqt/                      # 45 prepared PDBQT files (Meeko)
│
├── docking/
│   └── results/
│       ├── docking_scores.csv          # Vina scores for all 45 ligands
│       └── *_out.pdbqt                 # Docked poses (one per ligand)
│
├── analysis/
│   ├── interactions/
│   │   ├── interaction_analysis.csv    # Gly605/Glu603 contact counts
│   │   └── figures/
│   │       ├── 5LJJ-AD5.png            # Reversine in binding site (reference)
│   │       ├── 5LJJ-7CHN.png           # Best binder (7CHN) — dual hinge contact
│   │       └── 5LJJ-5N9S.png           # Non-binder (5N9S/5EHL) — no hinge contact
│   ├── adme/
│   │   ├── adme_full.csv               # All 45 ligands + descriptors + flags
│   │   └── adme_candidates.csv         # 26 candidates passing ADME + hinge
│   ├── ic50/
│   │   ├── chembl_mps1_ic50.csv        # 2,352 Mps1 IC50 values from ChEMBL
│   │   └── ic50_matches.csv            # Matches between our 45 and ChEMBL
│   └── ml_model/
│       ├── best_model.pkl              # Model 1 — Ridge Regression (Vina)
│       ├── scaler.pkl                  # Scaler for Model 1
│       ├── variance_selector.pkl       # Variance filter for Model 1
│       ├── predictions.csv             # Model 1 predicted vs actual
│       ├── ridge_coefficients.csv      # Ridge feature coefficients
│       └── chembl/
│           ├── chembl_best_model.pkl   # Model 2 — SVR (ChEMBL pIC50)
│           ├── chembl_scaler.pkl       # Scaler for Model 2
│           ├── chembl_variance_selector.pkl
│           ├── chembl_predictions.csv  # Model 2 predicted vs actual
│           ├── chembl_model_comparison.csv
│           └── chembl_svr_tuned.pkl    # Hyperparameter-tuned SVR
│
├── notebooks/
│   ├── 01_docking_analysis.ipynb       # Docking results exploration
│   ├── 02_adme_filter.ipynb            # ADME filtering visualisation
│   └── 03_ml_model.ipynb               # ML model analysis and plots
│
└── scripts/
    ├── download_ligands.py             # PubChem/RCSB automated download
    ├── cleanup_ligands.py              # File renaming and deduplication
    ├── prepare_ligands.py              # Meeko PDBQT preparation
    ├── run_docking.py                  # AutoDock Vina batch docking
    ├── parse_results.py                # Score extraction and ranking
    ├── interaction_analysis.py         # Gly605/Glu603 contact analysis
    ├── adme_filter.py                  # RDKit ADME descriptors + PAINS
    ├── ml_model.py                     # Model 1 — Vina score prediction
    ├── get_ic50.py                     # ChEMBL IC50 retrieval
    └── ml_chembl.py                    # Model 2 — pIC50 prediction (ChEMBL)
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

#### 6. Interaction analysis (Gly605/Glu603)
```bash
python scripts/interaction_analysis.py
```

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

### Phase 2 *(in progress)*

#### 11. Similarity search and novel candidate docking
```bash
python scripts/phase2_similarity_search.py  # coming soon
```
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
