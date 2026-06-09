import os

RAW_DIR = "data/ligands/raw"

# 1. Files to DELETE (artefacts + duplicates to remove)
to_delete = [
    # Artefacts — not real inhibitors
    "3WZK_CHLORIDE_ION.sdf",
    "4BI2_1,2-ETHANEDIOL.sdf",
    "4JT3_PHOSPHATE_ION.sdf",

    # Duplicates — keep the ~{N} cleaned version, delete the dirty one
    "2X9E.sdf",
    "2ZMD.sdf",
    "3GFW.sdf",
    "3HMP.sdf",
    "3W1F.sdf",
    "4C4E.sdf",
    "4C4G.sdf",
    "4C4H.sdf",
    "4C4I.sdf",
    "4C4J.sdf",

    # ~{N} dirty versions — superseded by clean renames below
    "5EHO_~{N}8-cyclohexyl-~{N}2-[2-methoxy-4-(1-methylpyrazol-4-yl)phenyl]pyrido[3,4-d]pyrimidine-2,8-diamine.sdf",
    "5EI2_~{N}-(2,4-dimethoxyphenyl)-8-(1-methylpyrazol-4-yl)pyrido[3,4-d]pyrimidin-2-amine.sdf",
    "5EI6_~{N}-(2,4-dimethoxyphenyl)-5-(1-methylpyrazol-4-yl)isoquinolin-3-amine.sdf",
    "5EI8_~{N}-[2-methoxy-4-(1-methylpyrazol-4-yl)phenyl]-8-(1-methylpyrazol-4-yl)pyrido[3,4-d]pyrimidin-2-amine.sdf",
    "5MRB_~{N}-(2,6-diethylphenyl)-8-[[2-methoxy-4-(4-methylpiperazin-1-yl)phenyl]amino]-1-methyl-4,5-dihydropyrazolo[4,3-h]quinazoline-3-carboxamide.sdf",
    "5N84_~{N}-cyclopropyl-4-[8-(2-methylpropylamino)imidazo[1,2-a]pyrazin-3-yl]benzamide.sdf",
    "5N87_~{N}-(2,6-diethylphenyl)-2-[[2-methoxy-4-(4-methylpiperazin-1-yl)phenyl]amino]-5,6-dihydropyrimido[4,5-e]indolizine-7-carboxamide.sdf",
    "5N9S_(2~{R})-2-(4-fluorophenyl)-~{N}-[4-[2-[(2-methoxy-4-methylsulfonyl-phenyl)amino]-[1,2,4]triazolo[1,5-a]pyridin-6-yl]phenyl]propanamide.sdf",
    "5N7V_~{N}6-cyclohexyl-~{N}2-(2-methyl-4-morpholin-4-yl-phenyl)-7~{H}-purine-2,6-diamine.sdf",
    "6H3K_~{N}8-(2,2-dimethylpropyl)-~{N}2-[2-ethoxy-4-(4-methyl-1,2,4-triazol-3-yl)phenyl]-6-methyl-pyrido[3,4-d]pyrimidine-2,8-diamine.sdf",
    "6TNB_(2~{R})-2-(4-fluorophenyl)-~{N}-[4-[2-[(2-methoxy-4-methylsulfonyl-phenyl)amino]-[1,2,4]triazolo[1,5-a]pyridin-6-yl]phenyl]propanamide.sdf",
    "7LQD_~{N}-(2,6-diethylphenyl)-2-[[4-(4-methylpiperazin-1-yl)-2-(propanoylamino)phenyl]amino]-5,6-dihydropyrimido[4,5-e]indolizine-7-carboxamide.sdf",
]

# 2. Files to RENAME > clean PDBID_COMPID.sdf format 
to_rename = {
    # long name → clean name
    "2X9E_N-(2,6-DIETHYLPHENYL)-1-METHYL-8-({4-[(1-METHYLPIPERIDIN-4-YL)CARBAMOYL]-2-(TRIFLUOROMETHOXY)PHENYL}AMINO)-4,5-DIHYDRO-1H-PYRAZOLO[4,3-H]QUINAZOLINE-3-CARBOXAMIDE.sdf": "2X9E_NMS.sdf",
    "2ZMD_537.sdf":                                          "2ZMD_537.sdf",  # already clean
    "3GFW_1-(4-(4-(2-(isopropylsulfonyl)phenylamino)-1H-pyrrolo[2,3-b]pyridin-6-ylamino)-3-methoxyphenyl)piperidin-4-ol.sdf": "3GFW_S22.sdf",
    "3HMP_7-chloro-N-(cyclopropylmethyl)quinazolin-4-amine.sdf": "3HMP_Q36.sdf",
    "3W1F_5-[5-ethoxy-6-(1-methyl-1H-pyrazol-4-yl)-1H-indazol-3-yl]-2-methylbenzenesulfonamide.sdf": "3W1F_LIG.sdf",
    "3WYY_(2E)-3-[4-({4-amino-5-cyano-6-[(3s,5s,7s)-tricyclo[3.3.1.1~3,7~]dec-1-ylamino]pyridin-2-yl}amino)-2-(cyanomethoxy)phenyl]-N-(2-methoxyethyl)prop-2-enamide.sdf": "3WYY_LIG.sdf",
    "3WZJ_4-{6-(cyclohexylamino)-8-[(tetrahydro-2H-pyran-4-ylmethyl)amino]imidazo[1,2-b]pyridazin-3-yl}-N-cyclopropylbenzamide.sdf": "3WZJ_LIG.sdf",
    "4C4E_N-(3,4-dimethoxyphenyl)-2-(1H-pyrazol-4-yl)-1H-pyrrolo[3,2-c]pyridin-6-amine.sdf": "4C4E_LIG.sdf",
    "4C4G_tert-butyl_6-((2-chloro-4-(dimethylcarbamoyl)phenyl)amino)-2-(1-methyl-1H-pyrazol-4-yl)-1H-pyrrolo[3,2-c]pyridine-1-carboxylate.sdf": "4C4G_LIG.sdf",
    "4C4H_tert-butyl_6-((2-chloro-4-(dimethylcarbamoyl)phenyl)amino)-2-(1-methyl-1H-pyrazol-4-yl)-1H-pyrrolo[3,2-c]pyridine-1-carboxylate.sdf": "4C4H_LIG.sdf",
    "4C4I_tert-butyl_6-{[2-chloro-4-(dimethylcarbamoyl)phenyl]amino}-2-(1,3-oxazol-5-yl)-1H-pyrrolo[3,2-c]pyridine-1-carboxylate.sdf": "4C4I_LIG.sdf",
    "4C4J_tert-butyl_6-{[2-methoxy-4-(4-methylpiperazin-1-yl)phenyl]amino}-2-(1-methyl-1H-pyrazol-4-yl)-1H-pyrrolo[3,2-c]pyridine-1-carboxylate.sdf": "4C4J_LIG.sdf",
    "5EHO_~{N}8-cyclohexyl-~{N}2-[2-methoxy-4-(1-methylpyrazol-4-yl)phenyl]pyrido[3,4-d]pyrimidine-2,8-diamine.sdf": "5EHO_LIG.sdf",
    "5EI2_N-(2,4-dimethoxyphenyl)-8-(1-methylpyrazol-4-yl)pyrido[3,4-d]pyrimidin-2-amine.sdf": "5EI2_LIG.sdf",
    "5EI6_N-(2,4-dimethoxyphenyl)-5-(1-methylpyrazol-4-yl)isoquinolin-3-amine.sdf": "5EI6_LIG.sdf",
    "5EI8_N-[2-methoxy-4-(1-methylpyrazol-4-yl)phenyl]-8-(1-methylpyrazol-4-yl)pyrido[3,4-d]pyrimidin-2-amine.sdf": "5EI8_LIG.sdf",
    "5LJJ_N~6~-cyclohexyl-N~2~-(4-morpholin-4-ylphenyl)-9H-purine-2,6-diamine.sdf": "5LJJ_AD5.sdf",
    "5MRB_N-(2,6-diethylphenyl)-8-[[2-methoxy-4-(4-methylpiperazin-1-yl)phenyl]amino]-1-methyl-4,5-dihydropyrazolo[4,3-h]quinazoline-3-carboxamide.sdf": "5MRB_LIG.sdf",
    "5N84_N-cyclopropyl-4-[8-(2-methylpropylamino)imidazo[1,2-a]pyrazin-3-yl]benzamide.sdf": "5N84_LIG.sdf",
    "5N87_N-(2,6-diethylphenyl)-2-[[2-methoxy-4-(4-methylpiperazin-1-yl)phenyl]amino]-5,6-dihydropyrimido[4,5-e]indolizine-7-carboxamide.sdf": "5N87_LIG.sdf",
    "5N9S_(2~{R})-2-(4-fluorophenyl)-N-[4-[2-[(2-methoxy-4-methylsulfonyl-phenyl)amino]-[1,2,4]triazolo[1,5-a]pyridin-6-yl]phenyl]propanamide.sdf": "5N9S_LIG.sdf",
    "5N7V_~{N}6-cyclohexyl-~{N}2-(2-methyl-4-morpholin-4-yl-phenyl)-7~{H}-purine-2,6-diamine.sdf": "5N7V_LIG.sdf",
    "5NA0_8QZ.sdf":                                          "5NA0_8QZ.sdf",  # already clean
    "5NAD_BAY_1217389.sdf":                                  "5NAD_BAY1217389.sdf",
    "5NTT_N-(2,6-DIETHYLPHENYL)-1-METHYL-8-({4-[(1-METHYLPIPERIDIN-4-YL)CARBAMOYL]-2-(TRIFLUOROMETHOXY)PHENYL}AMINO)-4,5-DIHYDRO-1H-PYRAZOLO[4,3-H]QUINAZOLINE-3-CARBOXAMIDE.sdf": "5NTT_LIG.sdf",
    "6H3K_~{N}8-(2,2-dimethylpropyl)-~{N}2-[2-ethoxy-4-(4-methyl-1,2,4-triazol-3-yl)phenyl]-6-methyl-pyrido[3,4-d]pyrimidine-2,8-diamine.sdf": "6H3K_LIG.sdf",
    "6TNB_(2~{R})-2-(4-fluorophenyl)-N-[4-[2-[(2-methoxy-4-methylsulfonyl-phenyl)amino]-[1,2,4]triazolo[1,5-a]pyridin-6-yl]phenyl]propanamide.sdf": "6TNB_LIG.sdf",
    "7LQD_N-(2,6-diethylphenyl)-2-[[4-(4-methylpiperazin-1-yl)-2-(propanoylamino)phenyl]amino]-5,6-dihydropyrimido[4,5-e]indolizine-7-carboxamide.sdf": "7LQD_LIG.sdf",
}

# Run deletions
print("=== Deleting artefacts and duplicates ===")
for fname in to_delete:
    fpath = os.path.join(RAW_DIR, fname)
    if os.path.exists(fpath):
        os.remove(fpath)
        print(f"   Deleted: {fname}")
    else:
        print(f"   Not found: {fname}")

# Run renames 
print("\n=== Renaming to clean format ===")
for old_name, new_name in to_rename.items():
    old_path = os.path.join(RAW_DIR, old_name)
    new_path = os.path.join(RAW_DIR, new_name)
    if old_name == new_name:
        print(f"  – Already clean: {new_name}")
        continue
    if os.path.exists(old_path):
        os.rename(old_path, new_path)
        print(f"   {old_name[:40]}... → {new_name}")
    else:
        print(f"   Not found: {old_name}")

#  Final count 
remaining = [f for f in os.listdir(RAW_DIR) if f.endswith(".sdf")]
print(f"\n Done. {len(remaining)} clean SDF files remaining:")
for f in sorted(remaining):
    print(f"  {f}")
