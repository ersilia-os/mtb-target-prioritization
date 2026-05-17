"""
Merges annotated protein features with the master UniProt list, derives target
prioritisation flags, and produces a ranked selection of 417 MTB drug targets.

Inputs:
    data/processed/02_master_uniprots.csv    — cleaned master UniProt table with
                                               essentiality, VI scores, and flags
                                               (produced by 02_data_annotation.py)
    data/processed/03_annotated_proteins.csv — per-protein annotations including
                                               ChEMBL count, PDB count, AlphaFold
                                               confidence, PANTHER classification,
                                               and protein evidence flag
                                               (produced by 03_annotate_proteins.py)

Outputs:
    data/processed/04_first_selection_417.csv  — ranked table of 417 prioritised MTB
                                                  targets with all annotation columns
    output/04_data_merging/04_upset_plot.png             — upset plot showing overlap between
                                                  target selection criteria
    output/04_data_merging/04_panther_annotation_counts.png — bar chart of PANTHER annotation
                                                  frequencies across the 417 targets
"""

import os
import sys

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))

import collections

import matplotlib.pyplot as plt
import pandas as pd
from upsetplot import UpSet, from_contents

data_dir = os.path.join(root, "..", "data", "processed")
plots_dir = os.path.join(root, "..", "output", "04_data_merging")
os.makedirs(data_dir, exist_ok=True)
os.makedirs(plots_dir, exist_ok=True)

df_0 = pd.read_csv(os.path.join(data_dir, "02_master_uniprots.csv"))
df_1 = pd.read_csv(os.path.join(data_dir, "03_annotated_proteins.csv"))

df = pd.merge(df_0, df_1, on="uniprot_ac", how="left")
df["vi"] = pd.to_numeric(df["vi"], errors="coerce")
df["high_confidence"] = df["high_confidence"].map(lambda x: str(x).strip().lower() == "true")

in_chembl = []
for v in df[["in_chembl", "chembl_count"]].values:
    if float(v[1]) > 0:
        in_chembl += [1]
    else:
        in_chembl += [int(v[0])]
df["in_chembl"] = in_chembl

is_vulnerable = []
for _, row in df.iterrows():
    if row["high_confidence"] == False:
        is_vulnerable += [0]
    else:
        if row["vi"] < -9.17575:
            is_vulnerable += [1]
        else:
            is_vulnerable += [0]
df["is_vulnerable"] = is_vulnerable

is_essential = []
for _, row in df.iterrows():
    if row["essentiality"] == "Essential" or row["essentiality"] == "Essential Domain":
        is_essential += [1]
    else:
        is_essential += [0]

is_growth_defect = []
for _, row in df.iterrows():
    if row["essentiality"] == "Growth Defect":
        is_growth_defect += [1]
    else:
        is_growth_defect += [0]

df["is_essential"] = is_essential
df["is_growth_defect"] = is_growth_defect

df = df.rename(columns={"high_confidence": "vi_hc"})
df = df[df["protein_evidence"] == True]

vi_hc = []
for v in df["vi_hc"].tolist():
    if v:
        vi_hc += [1]
    else:
        vi_hc += [0]
df["vi_hc"] = vi_hc

columns = [
    "uniprot_ac", "gene_name", "orf_id", "full_name", "protein_length",
    "is_known", "in_chembl", "is_essential", "is_growth_defect", "is_vulnerable",
    "chembl_count", "essentiality", "vi", "vi_hc", "pdb_count", "alphafold_conf",
    "panther_sf_name", "panther_annotation",
]
df = df[columns]

# Upset plot of target selection criteria overlap
data = collections.defaultdict(list)
for v in df[["uniprot_ac", "is_known", "in_chembl", "is_essential", "is_growth_defect", "is_vulnerable"]].values:
    if v[1] == 1:
        data["Known"] += [v[0]]
    if v[2] == 1:
        data["In ChEMBL"] += [v[0]]
    if v[3] == 1:
        data["Essential"] += [v[0]]
    if v[4] == 1:
        data["Growth Defect"] += [v[0]]
    if v[5] == 1:
        data["Vulnerable"] += [v[0]]

upset_data = from_contents(data)
UpSet(upset_data, subset_size="count", show_counts=True).plot()
plt.savefig(os.path.join(plots_dir, "04_upset_plot.png"), dpi=300)
plt.close()

# Target selection: known + ChEMBL + (vulnerable & essential) + (vulnerable & growth-defect)
# topped up with highest-VI essential proteins to reach 417 total
all_selected = df[df["is_known"] == 1]["uniprot_ac"].tolist()
all_selected += df[df["in_chembl"] == 1]["uniprot_ac"].tolist()
all_selected += df[(df["is_vulnerable"] == 1) & (df["is_essential"] == 1)]["uniprot_ac"].tolist()
all_selected += df[(df["is_vulnerable"] == 1) & (df["is_growth_defect"] == 1)]["uniprot_ac"].tolist()
all_selected = sorted(set(all_selected))

n = 417 - len(all_selected)
df_top = df[~df["uniprot_ac"].isin(all_selected)]
df_top = df_top.sort_values("vi", ascending=True)
df_top = df_top[df_top["is_essential"] == 1]
df_top = df_top.head(n)
all_selected += df_top["uniprot_ac"].tolist()

df_ = df[df["uniprot_ac"].isin(all_selected)]
df_ = df_.sort_values(
    by=["is_known", "in_chembl", "is_vulnerable", "is_essential", "is_growth_defect", "vi"],
    ascending=[False, False, False, False, False, True],
).reset_index()
columns = list(df_.columns)
df_["rank"] = [i + 1 for i in range(df_.shape[0])]
df_ = df_[["rank"] + columns]

df_.to_csv(os.path.join(data_dir, "04_first_selection_417.csv"), index=False)

# PANTHER annotation frequency bar chart for the 417 selected targets
df_["panther_annotation"].value_counts().head(25).plot(kind="bar", color="black", figsize=(8, 6))
plt.title("PANTHER Annotation Counts")
plt.xlabel("")
plt.ylabel("Frequency")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.savefig(os.path.join(plots_dir, "04_panther_annotation_counts.png"), dpi=300)
plt.close()
