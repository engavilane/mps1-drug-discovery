# Usage:
#    Phase 1 : python scripts/prepare_ligands.py
#    Phase 2 : python scripts/prepare_ligands \
#              --ligands data/phase2/pdbqt \
#              --results docking/phase2_results \
#              --exhaustiveness 8
    
import os
import subprocess

RAW_DIR    = "data/ligands/raw"
PDBQT_DIR  = "data/ligands/pdbqt"
LOG_FILE   = "data/ligands/preparation_log.txt"

os.makedirs(PDBQT_DIR, exist_ok=True)

sdf_files = sorted([f for f in os.listdir(RAW_DIR) if f.endswith(".sdf")])

print(f"Found {len(sdf_files)} SDF files to prepare\n")

success = []
failed  = []

with open(LOG_FILE, "w") as log:
    for sdf in sdf_files:
        name     = sdf.replace(".sdf", "")
        input_path  = os.path.join(RAW_DIR, sdf)
        output_path = os.path.join(PDBQT_DIR, f"{name}.pdbqt")

        cmd = [
            "mk_prepare_ligand.py",
            "-i", input_path,
            "-o", output_path
        ]

        print(f"Preparing: {sdf}")
        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode == 0 and os.path.exists(output_path):
            print(f"   Success: {name}.pdbqt")
            success.append(name)
            log.write(f"SUCCESS: {name}\n")
        else:
            print(f"   Failed: {name}")
            print(f"    {result.stderr.strip()}")
            failed.append(name)
            log.write(f"FAILED: {name}\n")
            log.write(f"  {result.stderr.strip()}\n")

# Summary 
print(f"\n{'='*50}")
print(f" Success: {len(success)}/{len(sdf_files)}")
print(f" Failed:  {len(failed)}/{len(sdf_files)}")

if failed:
    print("\nFailed ligands:")
    for f in failed:
        print(f"  - {f}")

print(f"\nLog saved to {LOG_FILE}")
