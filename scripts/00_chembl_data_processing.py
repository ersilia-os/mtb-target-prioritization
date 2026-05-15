"""
Extracts and processes ChEMBL bioactivity data for M. tuberculosis targets, mapping
ChEMBL target IDs to UniProt accessions (including manual curation of UNCHECKED assays),
building molecule-to-target activity matrices, and standardising SMILES strings.

Inputs:
    data/raw/mtuberculosis_ChEMBL_data.csv       — raw ChEMBL bioactivity records for
                                                    M. tuberculosis targets
    data/raw/chembl_uniprot_mapping.txt          — tab-separated ChEMBL target ID to
                                                    UniProt accession mapping file
    data/raw/chembl_unchecked_manual.csv         — MANUALLY GENERATED: target assignments
                                                    for UNCHECKED assays, curated by hand
                                                    from chembl_unchecked.csv; never
                                                    overwritten by this script

Outputs:
    data/processed/00_chembl/chembl_unchecked.csv             — activity counts per UNCHECKED
                                                                 assay, used as basis for
                                                                 manual curation
    data/processed/00_chembl/mtuberculosis_chembl_edited.csv  — ChEMBL data with UNCHECKED
                                                                 targets resolved
    data/processed/00_chembl/uniprot_ac_in_chembl.csv         — list of UniProt accessions
                                                                 present in the MTB ChEMBL
                                                                 dataset
    data/processed/00_chembl/chembl2uniprot.csv               — mapping of ChEMBL target IDs
                                                                 to UniProt accessions (or
                                                                 organism names for non-protein
                                                                 targets)
    data/processed/00_chembl/mols2target.csv                  — molecule × ChEMBL target
                                                                 activity count pivot table
    data/processed/00_chembl/mols2uniprot.csv                 — molecule × UniProt accession
                                                                 activity count pivot table
    data/processed/00_chembl/mols2uniprot_st.csv              — mols2uniprot with standardised
                                                                 SMILES and duplicate molecules
                                                                 collapsed
"""

import os
import sys

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))

import pandas as pd
from rdkit import Chem
from standardiser import standardise
from tqdm import tqdm

raw_dir = os.path.join(root, "..", "data", "raw")
chembl_dir = os.path.join(root, "..", "data", "processed", "00_chembl")
os.makedirs(chembl_dir, exist_ok=True)

# Non-protein or organism-level ChEMBL targets with no UniProt accession
non_protein_targets = {
    "CHEMBL2111188": "Mycobacterium tuberculosis H37Rv",
    "CHEMBL2364031": "Cell membrane",
    "CHEMBL360": "Mycobacterium tuberculosis",
    "CHEMBL3879801": "NON-PROTEIN TARGET",
    "CHEMBL612545": "UNCHECKED",
    "CHEMBL612558": "ADMET",
    "CHEMBL612960": "Mycobacterium tuberculosis variant microti",
    "CHEMBL613086": "Mycobacterium tuberculosis variant bovis",
    "CHEMBL614987": "Mycobacterium marinum",
    "CHEMBL615052": "Mycobacterium tuberculosis variant bovis BCG",
    "UNKNOWN": "UNCHECKED",
}

# Load raw ChEMBL data
df_raw = pd.read_csv(os.path.join(raw_dir, "mtuberculosis_ChEMBL_data.csv"))

# Export UNCHECKED assay counts for manual curation review
df_unchecked = df_raw[df_raw["target_type"] == "UNCHECKED"]
df_unchecked_counts = df_unchecked["assay_chembl_id"].value_counts().reset_index()
df_unchecked_counts.columns = ["assay_chembl_id", "counts"]
df_unchecked_counts.to_csv(os.path.join(chembl_dir, "chembl_unchecked.csv"), index=False)

# chembl_unchecked_manual.csv is a manually curated file and must never be overwritten.
# If it does not exist yet, review chembl_unchecked.csv and create it by hand before re-running.
manual_csv = os.path.join(raw_dir, "chembl_unchecked_manual.csv")
if not os.path.exists(manual_csv):
    raise FileNotFoundError(
        f"{manual_csv} not found. "
        "Manually curate chembl_unchecked.csv to produce this file before re-running."
    )

# Apply manual curation: map UNCHECKED assays to their correct ChEMBL target IDs
df_manual = pd.read_csv(manual_csv)
df_manual = df_manual[~df_manual["target_manual"].isna()]
df_manual = df_manual[["assay_chembl_id", "target_chembl_id"]].drop_duplicates(keep="first")

lookup = df_manual.set_index("assay_chembl_id")["target_chembl_id"]
mapped = df_raw["assay_chembl_id"].map(lookup)
df_raw["target_chembl_id"] = mapped.fillna(df_raw["target_chembl_id"])
df_raw.to_csv(os.path.join(chembl_dir, "mtuberculosis_chembl_edited.csv"), index=False)

# Build ChEMBL target ID → UniProt accession mapping
df_chembl_uniprot = pd.read_csv(os.path.join(raw_dir, "chembl_uniprot_mapping.txt"), sep="\t")
targets = list(set(df_raw["target_chembl_id"].tolist()))

df_mapped = (
    df_chembl_uniprot[df_chembl_uniprot["chembl_id"].isin(targets)]
    .dropna(subset=["uniprot_ac"])
    .rename(columns={"chembl_id": "target_chembl_id"})
)

# Add manually curated UniProt accessions for UNCHECKED assays
df_manual_uniprots = pd.read_csv(manual_csv)
df_manual_uniprots = df_manual_uniprots[~df_manual_uniprots["uniprot_ac"].isna()]
df_manual_uniprots = df_manual_uniprots[df_manual_uniprots["uniprot_ac"] != "UNKNOWN"]
df_manual_uniprots = df_manual_uniprots[["target_chembl_id", "uniprot_ac"]]

uniprot_ac_in_chembl = pd.DataFrame(
    {"uniprot_ac": list(set(df_manual_uniprots["uniprot_ac"].tolist() + df_mapped["uniprot_ac"].tolist()))}
)
uniprot_ac_in_chembl.to_csv(os.path.join(chembl_dir, "uniprot_ac_in_chembl.csv"), index=False)

# Build complete chembl_id → uniprot_ac table, including non-protein targets
df_chembl2uniprot = df_chembl_uniprot[df_chembl_uniprot["chembl_id"].isin(targets)][["chembl_id", "uniprot_ac"]]
df_non_protein = pd.DataFrame(
    [(k, v) for k, v in non_protein_targets.items()],
    columns=["chembl_id", "uniprot_ac"],
)
df_chembl2uniprot = pd.concat([df_chembl2uniprot, df_non_protein], ignore_index=True)
df_chembl2uniprot.to_csv(os.path.join(chembl_dir, "chembl2uniprot.csv"), index=False)

# Build molecule × ChEMBL target pivot table
df_edited = pd.read_csv(os.path.join(chembl_dir, "mtuberculosis_chembl_edited.csv"))
df_sub = df_edited[["canonical_smiles", "compound_chembl_id", "target_chembl_id"]].copy()

counts = (
    df_sub.groupby(["canonical_smiles", "compound_chembl_id", "target_chembl_id"])
    .size()
    .reset_index(name="count")
)

pivot = counts.pivot_table(
    index=["canonical_smiles", "compound_chembl_id"],
    columns="target_chembl_id",
    values="count",
    aggfunc="sum",
    fill_value=0,
)
pivot.columns = pivot.columns.astype(str)
pivot = pivot.reset_index().rename_axis(None, axis="index")
pivot.to_csv(os.path.join(chembl_dir, "mols2target.csv"), index=False)

# Map ChEMBL target columns to UniProt accessions
mapping = (
    df_chembl2uniprot.dropna(subset=["uniprot_ac"])
    .drop_duplicates(subset=["chembl_id", "uniprot_ac"])
    .groupby("chembl_id")["uniprot_ac"]
    .apply(list)
    .to_dict()
)


def merge_targets_to_uniprot_sum(pivot_df, mapping):
    """Remap pivot columns from ChEMBL IDs to UniProt ACs, summing where multiple
    ChEMBL IDs share a UniProt accession."""
    new_cols = {}
    for col in pivot_df.columns:
        if col not in mapping:
            new_cols[col] = new_cols.get(col, pivot_df[col]) if col in new_cols else pivot_df[col]
            continue
        uniprots = mapping[col]
        if not uniprots:
            new_cols[col] = new_cols.get(col, pivot_df[col]) if col in new_cols else pivot_df[col]
        else:
            for up in uniprots:
                if up in new_cols:
                    new_cols[up] = new_cols[up] + pivot_df[col]
                else:
                    new_cols[up] = pivot_df[col]
    return pd.DataFrame(new_cols, index=pivot_df.index)


df_uniprot = merge_targets_to_uniprot_sum(pivot, mapping).reset_index(drop=True)
df_uniprot.to_csv(os.path.join(chembl_dir, "mols2uniprot.csv"), index=False)

# Standardise SMILES and collapse duplicate molecules


def standardise_smiles(df, smi_col):
    std_smiles = []
    for smiles in tqdm(list(df[smi_col])):
        try:
            mol = Chem.MolFromSmiles(smiles)
            mol = standardise.run(mol)
            std_smiles.append(Chem.MolToSmiles(mol))
        except Exception:
            std_smiles.append(None)
    df["std_smiles"] = std_smiles
    return df[df["std_smiles"].notnull()]


df_uniprot = pd.read_csv(os.path.join(chembl_dir, "mols2uniprot.csv"))
df_std = standardise_smiles(df_uniprot, "canonical_smiles")
df_std = df_std.drop(columns=["canonical_smiles"])

value_cols = df_std.columns.difference(["std_smiles", "compound_chembl_id"])

compound_sums = df_std.groupby("compound_chembl_id")[value_cols].sum().sum(axis=1)


def pick_best_compound(series):
    return compound_sums.loc[series].idxmax()


collapsed = df_uniprot.groupby("std_smiles").agg(
    {**{col: "sum" for col in value_cols}, "compound_chembl_id": pick_best_compound}
).reset_index()

collapsed.to_csv(os.path.join(chembl_dir, "mols2uniprot_st.csv"), index=False)
