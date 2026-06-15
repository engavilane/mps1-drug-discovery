import pubchempy as pcp
import pandas as pd
import requests
import os
import time

# Paths
INPUT_CSV  = "kinase_domain/data/ligands/compounds.csv"
OUTPUT_DIR = "kinase_domain/data/ligands/raw"
LOG_FILE   = "kinase_domain/data/ligands/download_log.csv"

os.makedirs(OUTPUT_DIR, exist_ok=True)


# Load PDB ID list
df = pd.read_csv(INPUT_CSV)
pdb_ids = df["pdb_id"].tolist()

log = []

for pdb_id in pdb_ids:
    pdb_id = pdb_id.strip().upper()
    print(pdb_id)
