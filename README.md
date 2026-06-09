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
A dataset of **45 co-crystallised Mps1 inhibitors** (retrieved from the RCSB PDB) is subjected to:
1. Molecular docking against the Mps1 kinase domain (AutoDock Vina)
2. ADME filtering using RDKit (Lipinski's Rule of Five + drug-likeness descriptors)
3. Binding interaction analysis (H-bond contacts with Gly605/Glu603)
4. Machine learning model training to predict binding affinity from structural features

### Phase 2 — Novel Candidate Discovery
The best inhibitors identified in Phase 1 are used as seeds for:
1. PubChem structural similarity search (Tanimoto index ≥ 0.9)
2. Docking and ADME filtering of retrieved candidates
3. Affinity prediction using the Phase 1 trained ML model
4. Ranking and characterisation of novel predicted Mps1 inhibitors

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
ADME Filtering (RDKit)
    Lipinski Ro5 + solubility + PAINS
         │
         ▼
Interaction Analysis (Gly605 / Glu603 H-bonds)
         │
         ▼
ML Model (scikit-learn)
    Features: structural descriptors, fingerprints
    Target: predicted binding affinity (kcal/mol)
         │
         ▼
Phase 2: Similarity Search → Novel Candidates → Prediction
```

---

## Key Results (Phase 1 — Docking)

Top 5 predicted binders by AutoDock Vina score:

| Ligand | PDB ID | Affinity (kcal/mol) |
|--------|--------|---------------------|
| 6B4W_LIG | 6B4W | -10.844 |
| 7CHN_LIG | 7CHN | -10.324 |
| 4JS8_LIG | 4JS8 | -10.316 |
| 7CHM_LIG | 7CHM | -10.288 |
| 7CJA_LIG | 7CJA | -10.281 |

Reversine (5LJJ_AD5) scored **-9.645 kcal/mol**, consistent with published Glide XP scores (-10.97 kcal/mol, Pugh et al. 2022), validating the docking setup. Clinically advanced compounds BAY-1217389 and BAY-1161909 scored -9.91 and -9.644 kcal/mol respectively, in line with their known potency.

> Note: AutoDock Vina scores systematically overestimate binding affinities relative to experimental values. Results are interpreted in terms of relative ranking rather than absolute affinity prediction, consistent with established practice in the field (Pugh et al., 2022; Bolanos-Garcia, 2025).

---

## Repository Structure

```
docking-5LJJ/
├── README.md
├── environment.yml                  # Conda environment
├── data/
│   ├── raw/
│   │   └── 5LJJ.pdb                # Original PDB structure
│   ├── receptor/
│   │   ├── receptor_clean.pdb      # Waters/ligands removed
│   │   └── receptor.pdbqt          # Vina-ready receptor
│   └── ligands/
│       ├── compounds.csv           # PDB ID list
│       ├── download_log.csv        # Download status log
│       ├── preparation_log.txt     # Meeko preparation log
│       ├── native_ligand.pdb       # Reversine (AD5) reference
│       ├── raw/                    # 45 SDF files
│       └── pdbqt/                  # 45 prepared PDBQT files
├── docking/
│   └── results/                    # Vina output + docking_scores.csv
├── analysis/
│   ├── adme/                       # ADME filtering results
│   └── ml_model/                   # Model outputs and evaluation
├── notebooks/
│   ├── 01_docking_analysis.ipynb
│   ├── 02_adme_filter.ipynb
│   └── 03_ml_model.ipynb
└── scripts/
    ├── download_ligands.py         # PubChem/RCSB automated download
    ├── cleanup_ligands.py          # File renaming and deduplication
    ├── prepare_ligands.py          # Meeko PDBQT preparation
    ├── run_docking.py              # AutoDock Vina batch docking
    ├── parse_results.py            # Score extraction and ranking
    ├── adme_filter.py              # RDKit ADME descriptors
    └── ml_model.py                 # Affinity prediction model
```

---

## Requirements

Create and activate the conda environment:

```bash
conda env create -f environment.yml
conda activate docking
```

Key dependencies: `vina`, `meeko`, `rdkit`, `pubchempy`, `pandas`, `scikit-learn`, `requests`

---

## Reproducing the Pipeline

### 1. Download ligands
```bash
python scripts/download_ligands.py
```

### 2. Clean and rename ligand files
```bash
python scripts/cleanup_ligands.py
```

### 3. Prepare receptor
```bash
mk_prepare_receptor.py -i data/receptor/receptor_clean.pdb -o data/receptor/receptor -p
```

### 4. Prepare ligands
```bash
python scripts/prepare_ligands.py
```

### 5. Run docking
```bash
python scripts/run_docking.py
```

### 6. ADME filtering
```bash
python scripts/adme_filter.py
```

### 7. Train ML model
```bash
python scripts/ml_model.py
```

---

## References

- Bolanos-Garcia, V.M. (2025). Mps1 kinase functions in mitotic spindle assembly and error correction. *Trends in Biochemical Sciences*, 50(5), 438–453.
- Pugh, L. et al. (2022). Computational Biology Dynamics of Mps1 Kinase Molecular Interactions with Isoflavones Reveals a Chemical Scaffold with Potential to Develop New Therapeutics for the Treatment of Cancer. *Int. J. Mol. Sci.*, 23, 14228.
- Trott, O. & Olson, A.J. (2010). AutoDock Vina: improving the speed and accuracy of docking. *J. Comput. Chem.*, 31, 455–461.
- Daina, A. et al. (2017). SwissADME: A free web tool to evaluate pharmacokinetics, drug-likeness and medicinal chemistry friendliness of small molecules. *Sci. Rep.*, 7, 42717.

---

## Author

**Enya Gavilan Esquitino** — Master's in Bioinformatics Internship, 2025/26
