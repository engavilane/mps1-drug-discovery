import pubchempy as pcp
import pandas as pd
import requests
import os
import time

# Paths
INPUT_CSV  = "data/ligands/compounds.csv"
OUTPUT_DIR = "data/ligands/raw"
LOG_FILE   = "data/ligands/download_log.csv"

os.makedirs(OUTPUT_DIR, exist_ok=True)


# Load PDB ID list
df = pd.read_csv(INPUT_CSV)
pdb_ids = df["pdb_id"].tolist()

log = []

for pdb_id in pdb_ids:
    pdb_id = pdb_id.strip().upper()
    print(f"\nProcessing: {pdb_id}")

    # Step 1: Get ligand info from RCSB API
    try:
        url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
        r = requests.get(url)
        data = r.json()

        # Get the list of non-polymer entity IDs (ligands)
        ligand_ids = data.get("rcsb_entry_container_identifiers", {}) \
                         .get("non_polymer_entity_ids", [])

        if not ligand_ids:
            print(f"   No ligand found in {pdb_id}")
            log.append({"pdb_id": pdb_id, "ligand_id": None,
                        "ligand_name": None, "cid": None, "status": "no ligand"})
            continue

    except Exception as e:
        print(f"   RCSB API error for {pdb_id}: {e}")
        log.append({"pdb_id": pdb_id, "ligand_id": None,
                    "ligand_name": None, "cid": None, "status": f"RCSB error: {e}"})
        continue

    # Step 2: Get ligand name from its entity page 
    try:
        entity_id = ligand_ids[0]  # take the first ligand
        url2 = f"https://data.rcsb.org/rest/v1/core/nonpolymer_entity/{pdb_id}/{entity_id}"
        r2 = requests.get(url2)
        entity_data = r2.json()

        ligand_comp_id = entity_data["pdbx_entity_nonpoly"]["comp_id"]  
        ligand_name    = entity_data["pdbx_entity_nonpoly"]["name"]
        ligand_name = ligand_name.replace("~{N}-", "N-") \
                         .replace("~{O}-", "O-") \
                         .replace("~{S}-", "S-")
        print(f"  → Ligand found: {ligand_name} ({ligand_comp_id})")

    except Exception as e:
        print(f"   Could not extract ligand name for {pdb_id}: {e}")
        log.append({"pdb_id": pdb_id, "ligand_id": ligand_ids[0],
                    "ligand_name": None, "cid": None, "status": f"name error: {e}"})
        continue

    # Step 3: Search PubChem by ligand name 
    try:
        results = pcp.get_compounds(ligand_name, "name", record_type="3d")

        if not results:
            # fallback: try with the 3-letter comp_id
            results = pcp.get_compounds(ligand_comp_id, "name", record_type="3d")

        if not results:
            print(f"  → Trying RCSB direct download for {ligand_comp_id}...")
            rcsb_url = f"https://files.rcsb.org/ligands/download/{ligand_comp_id}_ideal.sdf"
            r = requests.get(rcsb_url)
            if r.status_code == 200:
                filepath = os.path.join(OUTPUT_DIR, f"{pdb_id}_{ligand_comp_id}.sdf")
                with open(filepath, "w") as f:
                    f.write(r.text)
                print(f"   Downloaded from RCSB: {ligand_comp_id}")
                log.append({"pdb_id": pdb_id, "ligand_id": ligand_comp_id,
                            "ligand_name": ligand_name, "cid": None, "status": "success (RCSB)"})
            else:
                print(f"   Not found anywhere: {ligand_comp_id}")
                log.append({"pdb_id": pdb_id, "ligand_id": ligand_comp_id,
                            "ligand_name": ligand_name, "cid": None, "status": "not found anywhere"})
            continue

        cid = results[0].cid

    except Exception as e:
        print(f"   PubChem search error for {ligand_name}: {e}")
        log.append({"pdb_id": pdb_id, "ligand_id": ligand_comp_id,
                    "ligand_name": ligand_name, "cid": None, "status": f"PubChem error: {e}"})
        continue

    # Step 4: Download 3D SDF 
    try:
        sdf_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=3d"
        response = requests.get(sdf_url)

        if response.status_code == 200:
            filename = ligand_name.replace(" ", "_").replace("/", "-")
            filepath = os.path.join(OUTPUT_DIR, f"{pdb_id}_{ligand_comp_id}.sdf")
            with open(filepath, "w") as f:
                f.write(response.text)
            print(f"   Downloaded: {ligand_name} (CID: {cid})")
            log.append({"pdb_id": pdb_id, "ligand_id": ligand_comp_id,
                        "ligand_name": ligand_name, "cid": cid, "status": "success"})
        else:
            print(f"   SDF download failed for {ligand_name}")
            log.append({"pdb_id": pdb_id, "ligand_id": ligand_comp_id,
                        "ligand_name": ligand_name, "cid": cid, "status": "SDF download failed"})

    except Exception as e:
        log.append({"pdb_id": pdb_id, "ligand_id": ligand_comp_id,
                    "ligand_name": ligand_name, "cid": cid, "status": f"download error: {e}"})

    time.sleep(0.5)


# Save log 
log_df = pd.DataFrame(log)
log_df.to_csv(LOG_FILE, index=False)

success = log_df[log_df.status.str.startswith("success")].shape[0]
print(f"\n Done. {success}/{len(pdb_ids)} ligands downloaded.")
print(f"Log saved to {LOG_FILE}")

