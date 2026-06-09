# Molecular Docking and Affinity Prediction of 5LJJ Inhibitors

## Project Overview
This project is part of a 1st year Master's internship in Bioinformatics.
It involves structure-based virtual screening of 51 inhibitors against 
protein 5LJJ using AutoDock Vina, followed by ADME filtering and the 
development of a machine learning model to predict binding affinity from 
structural features.

## Protein Target
- **PDB ID:** 5LJJ
- **Native ligand:** Reversine (AD5)
- **Binding site:** defined from co-crystallised reversine coordinates

## Pipeline Overview
1. Receptor preparation (Meeko)
2. Ligand preparation (Meeko + RDKit)
3. Molecular docking (AutoDock Vina 1.2)
4. ADME filtering (RDKit)
5. Affinity prediction model (scikit-learn)

## Requirements
Create and activate the conda environment:
```bash
conda env create -f environment.yml
conda activate docking
```

## Usage

### 1. Prepare receptor
```bash
mk_prepare_receptor.py -i data/receptor/receptor_clean.pdb -o data/receptor/receptor -p
```

### 2. Prepare ligands
```bash
python scripts/prepare_ligands.py
```

### 3. Run docking
```bash
python scripts/run_docking.py
```

### 4. Parse results
```bash
python scripts/parse_results.py
```

## Results
*(to be filled in as the project progresses)*

## Notebooks
- [Docking Analysis](notebooks/01_docking_analysis.ipynb)
- [ADME Filter](notebooks/02_adme_filter.ipynb)

## Author
Enya Gavilan Esquitino — Master's in Bioinformatics, [Université Clermont Auvergne (France)], [Oxford Brookes University (UK)], 2026
